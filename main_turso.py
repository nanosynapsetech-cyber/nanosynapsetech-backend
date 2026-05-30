# -*- coding: utf-8 -*-
# backend/main_turso.py -- NanoSynapse API v4.0  (R2 + Turso Hybrid)
#
# Architecture (3-layer):
#   [1] R2 nanosynapse-kmer/{seed_rc}.json  → gene_id list   (O(1) seed filter)
#   [2] Turso genes_meta HTTP API            → description, biotype
#   [3] R2 nanosynapse-seq/{gene_id}.txt     → full sequence  (thermodynamic scan)
#
# Required environment variables:
#   TURSO_DATABASE_URL   = libsql://nanosynapse-meta-xxxx.turso.io
#   TURSO_AUTH_TOKEN     = eyJh...
#   R2_ACCOUNT_ID        = <Cloudflare Account ID>
#   R2_ACCESS_KEY        = <R2 Access Key ID>
#   R2_SECRET_KEY        = <R2 Secret Access Key>
#
# Optional:
#   R2_KMER_BUCKET       = nanosynapse-kmer  (default)
#   R2_SEQ_BUCKET        = nanosynapse-seq   (default)
#   ALLOWED_ORIGINS      = https://nanosynapse.tech,http://localhost:3000

import os
import asyncio
import hashlib
import json
import logging
import time
import urllib.request
import urllib.parse
import urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional

try:
    import boto3
    from botocore.config import Config as BotocoreConfig
    _BOTO3_OK = True
except ImportError:
    _BOTO3_OK = False

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from mirna_engine import find_targets

# ─── LOGGING ──────────────────────────────────────────────────────────────────
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# ─── APP & CORS ───────────────────────────────────────────────────────────────
app = FastAPI(title="NanoSynapseTech API", version="4.0")

_raw_origins = os.environ.get("ALLOWED_ORIGINS", "http://localhost:3000")
ALLOWED_ORIGINS = [o.strip() for o in _raw_origins.split(",")]

app.add_middleware(
    CORSMiddleware,
    allow_origins=ALLOWED_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ─── CONFIG ───────────────────────────────────────────────────────────────────
TURSO_URL    = os.environ.get("TURSO_DATABASE_URL", "").strip()
TURSO_TOKEN  = os.environ.get("TURSO_AUTH_TOKEN",   "").strip()

R2_ACCOUNT_ID  = os.environ.get("R2_ACCOUNT_ID",  "").strip()
R2_ACCESS_KEY  = os.environ.get("R2_ACCESS_KEY",  "").strip()
R2_SECRET_KEY  = os.environ.get("R2_SECRET_KEY",  "").strip()
R2_KMER_BUCKET = os.environ.get("R2_KMER_BUCKET", "nanosynapse-kmer")
R2_SEQ_BUCKET  = os.environ.get("R2_SEQ_BUCKET",  "nanosynapse-seq")

# Derive Turso HTTP pipeline endpoint
_TURSO_HTTP_URL     = TURSO_URL.replace("libsql://", "https://") if TURSO_URL else ""
_TURSO_PIPELINE_URL = f"{_TURSO_HTTP_URL}/v2/pipeline" if _TURSO_HTTP_URL else ""

USE_TURSO = bool(TURSO_URL and TURSO_TOKEN)
USE_R2    = bool(R2_ACCOUNT_ID and R2_ACCESS_KEY and R2_SECRET_KEY and _BOTO3_OK)

# ─── S3 CLIENT (Cloudflare R2 is S3-compatible) ───────────────────────────────
_s3 = None
if USE_R2:
    _s3 = boto3.client(
        "s3",
        endpoint_url=f"https://{R2_ACCOUNT_ID}.r2.cloudflarestorage.com",
        aws_access_key_id=R2_ACCESS_KEY,
        aws_secret_access_key=R2_SECRET_KEY,
        config=BotocoreConfig(
            retries={"max_attempts": 3, "mode": "adaptive"},
            max_pool_connections=20,
        ),
        region_name="auto",
    )

if USE_TURSO:
    logger.info("✅ Turso HTTP API active — %s", TURSO_URL)
else:
    logger.warning("⚠️  TURSO_DATABASE_URL / TURSO_AUTH_TOKEN not set")

if USE_R2:
    logger.info("✅ R2 S3 API active — account: %s...", R2_ACCOUNT_ID[:6])
else:
    logger.warning("⚠️  R2_ACCOUNT_ID / R2_ACCESS_KEY / R2_SECRET_KEY not set (or boto3 missing)")

# ─── ORGANISM REGISTRY ────────────────────────────────────────────────────────
ORGANISM_DATABASES = {
    "human": {
        "display_name": "Homo sapiens (Human)",
        "common_name":  "Human",
        "taxon_id":     9606,
        "genome_build": "GRCh38/hg38",
        "notes":        "Human transcriptome — NCBI RefSeq mRNA sequences",
    },
}
DEFAULT_ORGANISM = "human"

logger.info("--- STARTING BIOLOGICAL ENGINE v4.0 (R2 + Turso Hybrid) ---")

# ─── IN-MEMORY LRU CACHES ─────────────────────────────────────────────────────
# Simple dict-based LRU (avoids functools.lru_cache limitations with lists)
_MAX_CACHE        = 200   # analysis result cache
_MAX_KMER_CACHE   = 500   # kmer JSON lists (each ~a few KB)
_MAX_SEQ_CACHE    = 300   # sequence strings (each ~2–10 KB)
_MAX_META_CACHE   = 1000  # metadata dicts

_search_cache: dict = {}
_kmer_cache:   dict = {}
_seq_cache:    dict = {}
_meta_cache:   dict = {}


def _lru_get(cache: dict, key: str):
    return cache.get(key)


def _lru_set(cache: dict, key: str, value, max_size: int):
    if len(cache) >= max_size:
        # Evict oldest entry (insertion-ordered dict)
        oldest = next(iter(cache))
        del cache[oldest]
    cache[key] = value


def _cache_key(mirna_id: str, mirna_seq: str,
               mfe_threshold: float, max_mismatches: int,
               rule_set: str, organism: str) -> str:
    raw = f"{mirna_id}:{mirna_seq}:{mfe_threshold}:{max_mismatches}:{rule_set}:{organism}"
    return hashlib.md5(raw.encode()).hexdigest()


# ─── R2 FETCH HELPERS (S3 API — no public access needed) ─────────────────────

def _s3_get(bucket: str, key: str) -> Optional[bytes]:
    """Fetch an object from R2 via S3 API. Returns bytes or None on miss/error."""
    if _s3 is None:
        return None
    try:
        resp = _s3.get_object(Bucket=bucket, Key=key)
        return resp["Body"].read()
    except _s3.exceptions.NoSuchKey:
        return None
    except Exception as e:
        # ClientError code 'NoSuchKey' also arrives here
        err_code = getattr(getattr(e, "response", {}), "get", lambda *a: None)("Error", {}).get("Code", "")
        if err_code not in ("NoSuchKey", "404"):
            logger.debug("R2 S3 error [%s/%s]: %s", bucket, key, e)
        return None


def fetch_kmer_list(seed_rc: str) -> list:
    """
    Fetch gene_id list for a k-mer seed from R2.
    Returns [] on miss or error.
    Cached in _kmer_cache.
    """
    cached = _lru_get(_kmer_cache, seed_rc)
    if cached is not None:
        return cached

    if not USE_R2:
        return []

    body = _s3_get(R2_KMER_BUCKET, f"{seed_rc}.json")
    if body is None:
        result = []
    else:
        try:
            result = json.loads(body.decode("utf-8"))
            if not isinstance(result, list):
                result = []
        except Exception:
            result = []

    _lru_set(_kmer_cache, seed_rc, result, _MAX_KMER_CACHE)
    return result


def fetch_sequence(gene_id: str) -> Optional[str]:
    """
    Fetch raw sequence string for a gene from R2.
    Returns None on miss/error.
    Cached in _seq_cache.
    """
    cached = _lru_get(_seq_cache, gene_id)
    if cached is not None:
        return cached

    if not USE_R2:
        return None

    body = _s3_get(R2_SEQ_BUCKET, f"{gene_id}.txt")
    if body is None:
        return None

    seq = body.decode("utf-8").strip().upper()
    _lru_set(_seq_cache, gene_id, seq, _MAX_SEQ_CACHE)
    return seq


# ─── TURSO METADATA HELPERS ───────────────────────────────────────────────────

def _turso_request(statements: list, retries: int = 3) -> dict:
    """POST to Turso HTTP pipeline API."""
    payload = json.dumps({"requests": statements}).encode("utf-8")
    req = urllib.request.Request(
        _TURSO_PIPELINE_URL,
        data=payload,
        headers={
            "Authorization": f"Bearer {TURSO_TOKEN}",
            "Content-Type":  "application/json",
        },
        method="POST",
    )
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=20) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            body = e.read().decode("utf-8")
            logger.warning("Turso HTTP %s: %s", e.code, body[:200])
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
        except Exception as e:
            logger.warning("Turso request error: %s", e)
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
    raise RuntimeError("Turso request failed after retries")


def fetch_metadata_batch(gene_ids: list) -> dict:
    """
    Fetch metadata for a list of gene_ids from Turso.
    Returns dict: {gene_id: {"description": ..., "biotype": ...}}
    Missing ids get default values. Results cached in _meta_cache.
    """
    if not gene_ids:
        return {}

    result      = {}
    missing_ids = []

    for gid in gene_ids:
        cached = _lru_get(_meta_cache, gid)
        if cached is not None:
            result[gid] = cached
        else:
            missing_ids.append(gid)

    if not missing_ids:
        return result

    if not USE_TURSO:
        # Return defaults when Turso not configured
        for gid in missing_ids:
            default = {"description": gid, "biotype": "transcript"}
            result[gid] = default
            _lru_set(_meta_cache, gid, default, _MAX_META_CACHE)
        return result

    # Batch query: up to 100 gene_ids per request (Turso pipeline limit)
    CHUNK = 80
    for i in range(0, len(missing_ids), CHUNK):
        chunk = missing_ids[i:i + CHUNK]
        placeholders = ",".join("?" * len(chunk))
        sql = (
            f"SELECT gene_id, description, biotype "
            f"FROM genes_meta WHERE gene_id IN ({placeholders})"
        )
        args = [{"type": "text", "value": gid} for gid in chunk]
        try:
            resp = _turso_request([
                {"type": "execute", "stmt": {"sql": sql, "args": args}},
                {"type": "close"},
            ])
            rows = (
                resp.get("results", [{}])[0]
                    .get("response", {})
                    .get("result", {})
                    .get("rows", [])
            )
            found = set()
            for row in rows:
                gid  = row[0]["value"]
                desc = row[1]["value"]
                bio  = row[2]["value"]
                meta = {"description": desc, "biotype": bio}
                result[gid] = meta
                found.add(gid)
                _lru_set(_meta_cache, gid, meta, _MAX_META_CACHE)

            # Fill defaults for ids not in Turso
            for gid in chunk:
                if gid not in found:
                    default = {"description": gid, "biotype": "transcript"}
                    result[gid] = default
                    _lru_set(_meta_cache, gid, default, _MAX_META_CACHE)

        except Exception as e:
            logger.error("Metadata batch fetch error: %s", e)
            for gid in chunk:
                if gid not in result:
                    result[gid] = {"description": gid, "biotype": "transcript"}

    return result


# ─── MODELS ───────────────────────────────────────────────────────────────────

class AnalysisRequest(BaseModel):
    mirna_fasta:     str
    target_fasta:    str   = ""
    search_mode:     str   = "manual"
    mfe_threshold:   float = -17.0
    max_mismatches:  int   = 4
    strict_cleavage: bool  = True
    rule_set:        str   = "animal"
    organism:        str   = DEFAULT_ORGANISM


class BatchRequest(BaseModel):
    mirnas:          list[str]
    search_mode:     str   = "automatic"
    mfe_threshold:   float = -17.0
    max_mismatches:  int   = 4
    strict_cleavage: bool  = True
    rule_set:        str   = "animal"
    organism:        str   = DEFAULT_ORGANISM


# ─── CORE DATABASE SEARCH (R2 + Turso Hybrid) ─────────────────────────────────

def run_database_search(mirna_id: str, mirna_seq: str,
                        mfe_threshold: float, max_mismatches: int,
                        strict_cleavage: bool,
                        rule_set: str = "animal",
                        organism: str = DEFAULT_ORGANISM) -> list:
    """
    3-layer hybrid search:
      1. R2 kmer lookup  → gene_id candidates
      2. Turso batch     → description + biotype metadata
      3. R2 seq fetch    → full sequences for thermodynamic analysis

    NOTE: For animal rule_set, strict_cleavage is automatically False.
    """
    org_info           = ORGANISM_DATABASES.get(organism, ORGANISM_DATABASES[DEFAULT_ORGANISM])
    org_display        = org_info["display_name"]
    effective_cleavage = strict_cleavage if rule_set == "plant" else False

    cache_key = _cache_key(mirna_id, mirna_seq, mfe_threshold, max_mismatches, rule_set, organism)
    cached    = _lru_get(_search_cache, cache_key)
    if cached is not None:
        logger.info("✨ Cache HIT: %s [%s]", mirna_id, org_display)
        return cached

    logger.info("🚀 Scan: %s | %s | R2+Turso hybrid", mirna_id, org_display)

    # ── Step 1: Seed RC calculation ───────────────────────────────────────────
    mirna_seed = mirna_seq[1:8] if len(mirna_seq) >= 8 else mirna_seq
    rc_map     = {"A": "U", "U": "A", "G": "C", "C": "G"}
    seed_rc    = "".join(rc_map.get(c, c) for c in reversed(mirna_seed))

    # ── Step 2: R2 kmer lookup → candidate gene_id list ──────────────────────
    gene_ids = fetch_kmer_list(seed_rc)
    logger.info("  kmer '%s' → %d candidates", seed_rc, len(gene_ids))

    if not gene_ids:
        logger.info("  No candidates found for seed %s", seed_rc)
        _lru_set(_search_cache, cache_key, [], _MAX_CACHE)
        return []

    # Cap candidates — configurable via env MAX_CANDIDATES (default 100)
    MAX_CANDIDATES = int(os.environ.get("MAX_CANDIDATES", "100"))
    if len(gene_ids) > MAX_CANDIDATES:
        import random
        gene_ids = random.sample(gene_ids, MAX_CANDIDATES)
        logger.info("  Sampled %d / original candidates", MAX_CANDIDATES)

    # ── Step 3: Parallel R2 sequence fetch + thermodynamic analysis ───────────
    raw_results = []   # list of {gene_id, prediction}
    scanned     = 0
    deep        = 0

    MAX_SEQ_WORKERS = 32
    MAX_TARGETS     = 15

    def _analyze_gene(gene_id: str) -> Optional[dict]:
        """Fetch sequence and run thermodynamic analysis for one gene."""
        target_seq = fetch_sequence(gene_id)
        if not target_seq or len(target_seq) <= 10 or "N" in target_seq:
            return None

        # Locate seed match window (±50 nt around first hit)
        match_idx = target_seq.find(seed_rc)
        if match_idx < 0:
            short_target = target_seq[:200]
        else:
            start_idx    = max(0, match_idx - 50)
            end_idx      = min(len(target_seq), match_idx + 50)
            short_target = target_seq[start_idx:end_idx]

        engine_result = find_targets(
            mirna_id, mirna_seq,
            gene_id, short_target,
            compute_pvalue=False,
            rule_set=rule_set,
            max_mismatches=max_mismatches,
            strict_cleavage=effective_cleavage,
            mfe_threshold=mfe_threshold,
        )

        if engine_result.get("Status") != "PASS":
            return None
        if engine_result.get("Similarity_Percent", 0) < 55.0:
            return None

        # Return intermediate result — metadata fetched later in batch
        return {"gene_id": gene_id, "prediction": engine_result}

    # Submit all genes to thread pool; collect until we have 15 matches
    with ThreadPoolExecutor(max_workers=MAX_SEQ_WORKERS) as pool:
        future_map = {pool.submit(_analyze_gene, gid): gid for gid in gene_ids}
        for future in as_completed(future_map):
            scanned += 1
            try:
                hit = future.result()
            except Exception as e:
                logger.debug("Analysis error for %s: %s", future_map[future], e)
                continue

            if hit is not None:
                deep += 1
                raw_results.append(hit)
                logger.info(
                    "  ✅ MATCH: %s | Confidence: %s | MFE: %s kcal/mol",
                    future_map[future],
                    hit["prediction"].get("Confidence_Score", "?"),
                    hit["prediction"].get("MFE_kcal_mol", "?"),
                )
                if len(raw_results) >= MAX_TARGETS:
                    for f in future_map:
                        f.cancel()
                    break

    # ── Step 4: Turso metadata batch fetch (only for matched genes) ─────────
    matched_ids = [h["gene_id"] for h in raw_results]
    meta_map    = fetch_metadata_batch(matched_ids)

    results = []
    for h in raw_results:
        gid  = h["gene_id"]
        meta = meta_map.get(gid, {"description": gid, "biotype": "transcript"})
        results.append({
            "prediction": h["prediction"],
            "biological_context": {
                "description":   meta["description"],
                "biotype":       meta["biotype"],
                "database":      "NCBI/Ensembl (R2 S3 + Turso)",
                "organism":      org_info["common_name"],
                "organism_full": org_info["display_name"],
                "genome_build":  org_info.get("genome_build", "N/A"),
                "taxon_id":      org_info.get("taxon_id", 0),
            },
        })

    logger.info(
        "Scan done — Candidates: %d, Fetched: %d, Matches: %d",
        len(gene_ids), scanned, len(results),
    )

    results.sort(key=lambda x: x["prediction"].get("Confidence_Score", 0), reverse=True)
    _lru_set(_search_cache, cache_key, results, _MAX_CACHE)
    return results


# ─── HELPER: PARSE FASTA ──────────────────────────────────────────────────────

def parse_fasta(fasta_str: str) -> tuple:
    lines  = fasta_str.strip().split("\n")
    seq_id = lines[0].lstrip(">").strip() if lines[0].startswith(">") else "Unknown"
    seq    = "".join(lines[1:]).replace(" ", "").upper()
    return seq_id, seq


# ─── ENDPOINT: /api/analyze ───────────────────────────────────────────────────

@app.post("/api/analyze")
async def run_analysis(request: AnalysisRequest):
    try:
        mirna_id, mirna_seq = parse_fasta(request.mirna_fasta)

        if request.search_mode == "automatic":
            results = await asyncio.to_thread(
                run_database_search,
                mirna_id, mirna_seq,
                request.mfe_threshold,
                request.max_mismatches,
                request.strict_cleavage,
                request.rule_set,
                request.organism,
            )
        else:
            target_id, target_seq = parse_fasta(request.target_fasta)
            effective_cleavage = request.strict_cleavage if request.rule_set == "plant" else False
            engine_result = find_targets(
                mirna_id, mirna_seq, target_id, target_seq,
                compute_pvalue=True,
                rule_set=request.rule_set,
                max_mismatches=request.max_mismatches,
                strict_cleavage=effective_cleavage,
                mfe_threshold=request.mfe_threshold,
            )
            results = [{
                "prediction": engine_result,
                "biological_context": {
                    "description":   "Manual Target Sequence",
                    "biotype":       "N/A",
                    "database":      "Manual Input",
                    "organism":      "Manual",
                    "organism_full": "Manual Input (user-provided sequence)",
                    "genome_build":  "N/A",
                    "taxon_id":      0,
                },
            }]

        return {"status": "success", "data": results}

    except Exception as e:
        logger.error("Analysis failed: %s", e)
        raise HTTPException(status_code=500, detail=str(e))


# ─── ENDPOINT: /api/batch-analyze ─────────────────────────────────────────────

@app.post("/api/batch-analyze")
async def batch_analysis(request: BatchRequest):
    if len(request.mirnas) > 20:
        raise HTTPException(status_code=400,
                            detail="Maximum 20 miRNAs per batch request.")
    try:
        tasks = [
            asyncio.to_thread(
                run_database_search,
                *parse_fasta(fasta),
                request.mfe_threshold,
                request.max_mismatches,
                request.strict_cleavage,
                request.rule_set,
                request.organism,
            )
            for fasta in request.mirnas
        ]

        all_results = await asyncio.gather(*tasks, return_exceptions=True)

        batch_out = []
        for fasta, res in zip(request.mirnas, all_results):
            mirna_id, _ = parse_fasta(fasta)
            if isinstance(res, Exception):
                batch_out.append({"mirna_id": mirna_id, "error": str(res), "data": []})
            else:
                batch_out.append({"mirna_id": mirna_id, "data": res, "match_count": len(res)})

        return {"status": "success", "batch_results": batch_out}

    except Exception as e:
        logger.error("Batch analysis failed: %s", e)
        raise HTTPException(status_code=500, detail=str(e))


# ─── ENDPOINT: /api/mirna-lookup/{name} ───────────────────────────────────────

def _fetch_mirbase_sync(mirna_name: str) -> dict:
    try:
        encoded = urllib.parse.quote(mirna_name)
        url = f"https://mirbase.org/api/v1/mirna/?name={encoded}"
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode())
            if data and isinstance(data, list) and len(data) > 0:
                entry = data[0]
                return {
                    "found":    True,
                    "source":   "miRBase",
                    "id":       entry.get("id", mirna_name),
                    "name":     entry.get("name", mirna_name),
                    "sequence": entry.get("sequence", ""),
                    "organism": entry.get("organism", ""),
                }
    except Exception:
        pass

    try:
        term   = urllib.parse.quote(f"{mirna_name}[Gene Name]")
        search = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            f"?db=gene&term={term}&retmode=json&retmax=1"
        )
        with urllib.request.urlopen(search, timeout=10) as resp:
            s_data = json.loads(resp.read().decode())
            ids    = s_data.get("esearchresult", {}).get("idlist", [])
            if ids:
                return {"found": True, "source": "NCBI",
                        "id": ids[0], "name": mirna_name, "sequence": ""}
    except Exception:
        pass

    return {"found": False, "source": None,
            "id": None, "name": mirna_name, "sequence": None}


@app.get("/api/mirna-lookup/{mirna_name}")
async def mirna_lookup(mirna_name: str):
    result = await asyncio.to_thread(_fetch_mirbase_sync, mirna_name)
    if not result["found"]:
        raise HTTPException(status_code=404,
                            detail=f"miRNA '{mirna_name}' not found in miRBase or NCBI.")
    return result


# ─── ENDPOINT: /api/cache-stats ───────────────────────────────────────────────

@app.get("/api/cache-stats")
async def cache_stats():
    return {
        "search_cache":  len(_search_cache),
        "kmer_cache":    len(_kmer_cache),
        "seq_cache":     len(_seq_cache),
        "meta_cache":    len(_meta_cache),
        "limits": {
            "search": _MAX_CACHE,
            "kmer":   _MAX_KMER_CACHE,
            "seq":    _MAX_SEQ_CACHE,
            "meta":   _MAX_META_CACHE,
        },
    }


# ─── ENDPOINT: /api/organisms ─────────────────────────────────────────────────

@app.get("/api/organisms")
async def list_organisms():
    return {
        "organisms": [
            {
                "id":           org_id,
                "display_name": info["display_name"],
                "common_name":  info["common_name"],
                "taxon_id":     info.get("taxon_id", 0),
                "genome_build": info.get("genome_build", "N/A"),
                "notes":        info.get("notes", ""),
                "available":    True,  # R2 is always online
                "db_size_mb":   0,     # No local DB in hybrid mode
            }
            for org_id, info in ORGANISM_DATABASES.items()
        ],
        "default": DEFAULT_ORGANISM,
    }


# ─── ENDPOINT: /api/health ────────────────────────────────────────────────────

@app.get("/api/health")
async def health():
    return {
        "status":           "ok",
        "engine":           "v4.0",
        "database_backend": "r2_s3+turso_hybrid",
        "turso_url":        TURSO_URL or None,
        "turso_active":     USE_TURSO,
        "r2_account":       (R2_ACCOUNT_ID[:6] + "...") if R2_ACCOUNT_ID else None,
        "r2_kmer_bucket":   R2_KMER_BUCKET,
        "r2_seq_bucket":    R2_SEQ_BUCKET,
        "r2_active":        USE_R2,
        "organisms":        {"human": True},
        "default_organism": DEFAULT_ORGANISM,
        "cache": {
            "search": len(_search_cache),
            "kmer":   len(_kmer_cache),
            "seq":    len(_seq_cache),
            "meta":   len(_meta_cache),
        },
    }
