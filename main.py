# -*- coding: utf-8 -*-
# backend/main.py -- NanoSynapse API v2.1
# New in v2:
#   [4] In-memory result cache (MD5 keyed)
#   [7] Batch analysis endpoint  /api/batch-analyze
#   [8] miRBase lookup endpoint  /api/mirna-lookup/{name}

import sqlite3
import os
import asyncio
import hashlib
import json
import logging
import urllib.request
import urllib.parse
from typing import Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from mirna_engine import find_targets

# ─── LOGGING ──────────────────────────────────────────────────────────────────
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# ─── APP & CORS ───────────────────────────────────────────────────────────────
app = FastAPI(title="NanoSynapseTech API", version="2.0")

_raw_origins = os.environ.get("ALLOWED_ORIGINS", "http://localhost:3000")
ALLOWED_ORIGINS = [o.strip() for o in _raw_origins.split(",")]

app.add_middleware(
    CORSMiddleware,
    allow_origins=ALLOWED_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ─── ORGANISM DATABASE REGISTRY ──────────────────────────────────────────────
# Each entry maps an organism_id → { db_path, display_name, common_name, taxon_id }
# To add a new organism: add a new entry here and rebuild its SQLite database.
ORGANISM_DATABASES: dict = {
    "human": {
        "db_path":     "database.db",
        "display_name": "Homo sapiens (Human)",
        "common_name":  "Human",
        "taxon_id":    9606,
        "genome_build": "GRCh38/hg38",
        "notes":       "Human transcriptome — NCBI RefSeq mRNA sequences",
    },
    # Future organisms — add DB, register here, rebuild FTS index:
    # "mouse": {
    #     "db_path":     "mouse.db",
    #     "display_name": "Mus musculus (Mouse)",
    #     "common_name":  "Mouse",
    #     "taxon_id":    10090,
    #     "genome_build": "GRCm39/mm39",
    #     "notes":       "Mouse transcriptome",
    # },
    # "arabidopsis": {
    #     "db_path":     "arabidopsis.db",
    #     "display_name": "Arabidopsis thaliana",
    #     "common_name":  "Arabidopsis",
    #     "taxon_id":    3702,
    #     "genome_build": "TAIR10",
    #     "notes":       "Arabidopsis thaliana transcriptome (for plant rule set)",
    # },
}

DEFAULT_ORGANISM = "human"

# ─── DATABASE ─────────────────────────────────────────────────────────────────
DB_PATH = ORGANISM_DATABASES[DEFAULT_ORGANISM]["db_path"]

logger.info("--- STARTING BIOLOGICAL ENGINE v2.0 ---")
for org_id, org_info in ORGANISM_DATABASES.items():
    if os.path.exists(org_info["db_path"]):
        logger.info(f"✅ Organism DB ready: {org_info['display_name']} → {org_info['db_path']}")
    else:
        logger.warning(f"⚠️  Organism DB missing: {org_info['display_name']} → {org_info['db_path']} (run csv_to_sqlite.py)")

# ─── IN-MEMORY CACHE ──────────────────────────────────────────────────────────
_search_cache = {}
_MAX_CACHE = 100   # Keep last 100 unique queries


def _cache_key(mirna_id: str, mirna_seq: str,
               mfe_threshold: float, max_mismatches: int, rule_set: str,
               rule_max_mismatch: bool = True, rule_cleavage_site: bool = True,
               rule_seed_mismatch: bool = True, rule_consecutive: bool = True,
               rule_seed_mfe: bool = False, seed_mfe_threshold: float = -20.0) -> str:
    raw = (
        f"{mirna_id}:{mirna_seq}:{mfe_threshold}:{max_mismatches}:{rule_set}"
        f":{rule_max_mismatch}:{rule_cleavage_site}:{rule_seed_mismatch}"
        f":{rule_consecutive}:{rule_seed_mfe}:{seed_mfe_threshold}"
    )
    return hashlib.md5(raw.encode()).hexdigest()


def _cache_get(key: str) -> Optional[list]:
    return _search_cache.get(key)


def _cache_set(key: str, value: list) -> None:
    if len(_search_cache) >= _MAX_CACHE:
        # Evict oldest entry
        oldest = next(iter(_search_cache))
        del _search_cache[oldest]
    _search_cache[key] = value

# ─── MODELS ───────────────────────────────────────────────────────────────────

class AnalysisRequest(BaseModel):
    mirna_fasta:        str
    target_fasta:       str   = ""
    search_mode:        str   = "manual"
    mfe_threshold:      float = -17.0
    max_mismatches:     int   = 4
    strict_cleavage:    bool  = True
    rule_set:           str   = "animal"
    organism:           str   = DEFAULT_ORGANISM
    # v3.0: Granüler kural kontroller (her biri bağımsız toggle)
    rule_max_mismatch:  bool  = True    # Kural i:   Global max mismatch
    rule_cleavage_site: bool  = True    # Kural ii:  Pos 10-11 koruması
    rule_seed_mismatch: bool  = True    # Kural iii: Pos 1-9 max 1 mismatch
    rule_consecutive:   bool  = True    # Kural iv:  Ardışık 2+ yasak
    rule_seed_mfe:      bool  = False   # Kural v:   Seed region MFE (default kapalı)
    seed_mfe_threshold: float = -20.0   # Kural v eşiği (kcal/mol)


class BatchRequest(BaseModel):
    mirnas:             list[str]
    search_mode:        str   = "automatic"
    mfe_threshold:      float = -17.0
    max_mismatches:     int   = 4
    strict_cleavage:    bool  = True
    rule_set:           str   = "animal"
    organism:           str   = DEFAULT_ORGANISM
    # v3.0: Granüler kural kontroller
    rule_max_mismatch:  bool  = True
    rule_cleavage_site: bool  = True
    rule_seed_mismatch: bool  = True
    rule_consecutive:   bool  = True
    rule_seed_mfe:      bool  = False
    seed_mfe_threshold: float = -20.0


# ─── CORE DATABASE SEARCH ─────────────────────────────────────────────────────

def run_database_search(mirna_id: str, mirna_seq: str,
                         mfe_threshold: float, max_mismatches: int,
                         strict_cleavage: bool,
                         rule_set: str = "animal",
                         organism: str = DEFAULT_ORGANISM,
                         # v3.0: Granüler kural parametreleri
                         rule_max_mismatch: bool = True,
                         rule_cleavage_site: bool = True,
                         rule_seed_mismatch: bool = True,
                         rule_consecutive: bool = True,
                         rule_seed_mfe: bool = False,
                         seed_mfe_threshold: float = -20.0) -> list:
    """
    Seed-filter → ViennaRNA thermodynamic scan → top-15 results.
    Results are cached by query hash.
    
    NOTE: For animal mode, strict_cleavage is automatically set to False
    because AGO2-mediated slicing is exceptional in animals; most targeting
    occurs via translational repression (no perfect cleavage site needed).
    """
    # Resolve organism DB path
    org_info = ORGANISM_DATABASES.get(organism, ORGANISM_DATABASES[DEFAULT_ORGANISM])
    db_path  = org_info["db_path"]
    org_display = org_info["display_name"]

    if not os.path.exists(db_path):
        logger.error(f"Organism DB not found: {db_path}")
        return []

    # v3.0: strict_cleavage artık kullanıcı tercihine bırakıldı (animal override kaldırıldı)
    effective_cleavage = strict_cleavage
    key = _cache_key(mirna_id, mirna_seq, mfe_threshold, max_mismatches, rule_set,
                     rule_max_mismatch, rule_cleavage_site, rule_seed_mismatch,
                     rule_consecutive, rule_seed_mfe, seed_mfe_threshold) + f":{organism}"
    cached = _cache_get(key)
    if cached is not None:
        logger.info(f"✨ Cache HIT: {mirna_id} [{org_display}]")
        return cached

    results = []
    logger.info(f"🚀 Scan started: {mirna_id} | Organism: {org_display}")

    # Seed region (pos 2-8) → reverse complement for SQL LIKE filter
    mirna_seed = mirna_seq[1:8] if len(mirna_seq) >= 8 else mirna_seq
    rc_map     = {"A": "U", "U": "A", "G": "C", "C": "G"}
    seed_rc    = "".join(rc_map.get(c, c) for c in reversed(mirna_seed))

    scanned = 0
    deep    = 0

    try:
        conn   = sqlite3.connect(db_path)
        cursor = conn.cursor()

        query = (
            "SELECT g.* FROM genes g "
            "JOIN genes_fts f ON g.rowid = f.rowid "
            "WHERE f.Sequence MATCH ? "
            "AND (g.Biotype = 'transcript' OR g.Biotype LIKE '%mrna%' "
            "     OR g.Biotype LIKE '%cdna%') "
            "AND NOT (g.Gene_ID LIKE 'XM_%' OR g.Gene_ID LIKE 'XR_%' OR g.Description LIKE 'PREDICTED%')"
        )
        cursor.execute(query, (f'"{seed_rc}"',))
        columns = [col[0] for col in cursor.description]

        while True:
            row = cursor.fetchone()
            if not row:
                break

            row_dict    = dict(zip(columns, row))
            target_id   = str(row_dict.get("Gene_ID",     "Unknown"))
            target_seq  = str(row_dict.get("Sequence",    "")).strip().upper()
            description = str(row_dict.get("Description", "Target Gene"))
            biotype     = str(row_dict.get("Biotype",     "N/A"))

            scanned += 1
            if len(target_seq) <= 10 or "N" in target_seq:
                continue

            deep += 1

            # Trim to 100 nt window around seed match
            match_idx   = target_seq.find(seed_rc)
            start_idx   = max(0, match_idx - 50)
            end_idx     = min(len(target_seq), match_idx + 50)
            short_target = target_seq[start_idx:end_idx]

            engine_result = find_targets(mirna_id, mirna_seq,
                                          target_id, short_target,
                                          compute_pvalue=False,
                                          rule_set=rule_set,
                                          max_mismatches=max_mismatches,
                                          strict_cleavage=effective_cleavage,
                                          mfe_threshold=mfe_threshold,
                                          rule_max_mismatch=rule_max_mismatch,
                                          rule_cleavage_site=rule_cleavage_site,
                                          rule_seed_mismatch=rule_seed_mismatch,
                                          rule_consecutive=rule_consecutive,
                                          # Kural v: önce PASS'leri bul, sonra seed MFE uygula
                                          rule_seed_mfe=False,
                                          seed_mfe_threshold=seed_mfe_threshold)

            if engine_result["Status"] == "PASS":
                sim = engine_result.get("Similarity_Percent", 0)
                if sim >= 55.0:
                    # Kural v: PASS + similarity OK ise seed MFE post-filter uygula
                    if rule_seed_mfe:
                        from mirna_engine import calculate_region_mfe
                        bind_pos = int(engine_result["Binding_Position"].split("-")[0])
                        s_mfe_2_7, s_mfe_8_13 = calculate_region_mfe(
                            mirna_seq, short_target, bind_pos
                        )
                        engine_result["Seed_MFE_2_7"]  = s_mfe_2_7
                        engine_result["Seed_MFE_8_13"] = s_mfe_8_13
                        if not (s_mfe_2_7 <= seed_mfe_threshold or s_mfe_8_13 <= seed_mfe_threshold):
                            logger.info(
                                f"   Kural v FAIL: {target_id} seed MFE 2-7={s_mfe_2_7}, 8-13={s_mfe_8_13}"
                            )
                            continue  # Bu geni sonuca ekleme

                    results.append({
                        "prediction": engine_result,
                        "biological_context": {
                            "description":   description,
                            "biotype":        biotype,
                            "database":       "NCBI/Ensembl",
                            "organism":       org_info["common_name"],
                            "organism_full":  org_info["display_name"],
                            "genome_build":   org_info.get("genome_build", "N/A"),
                            "taxon_id":       org_info.get("taxon_id", 0),
                        }
                    })
                    logger.info(
                        f"✅ MATCH: {target_id} "
                        f"| Confidence: {engine_result['Confidence_Score']}/100 "
                        f"| MFE: {engine_result['MFE_kcal_mol']} kcal/mol"
                    )

            if len(results) >= 15:
                break

        logger.info(
            f"Scan done — SQL: {scanned} genes, "
            f"Thermodynamic: {deep}, Matches: {len(results)}"
        )

    except Exception as e:
        logger.error(f"Database Search Error: {e}")
    finally:
        if 'conn' in locals():
            conn.close()

    # Sort results by confidence score (highest first)
    results.sort(key=lambda x: x["prediction"].get("Confidence_Score", 0),
                 reverse=True)

    _cache_set(key, results)
    return results


# ─── HELPER: PARSE FASTA STRING ───────────────────────────────────────────────

def parse_fasta(fasta_str: str) -> tuple[str, str]:
    lines   = fasta_str.strip().split('\n')
    seq_id  = lines[0].lstrip('>').strip() if lines[0].startswith('>') else "Unknown"
    seq     = "".join(lines[1:]).replace(" ", "").upper()
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
                request.rule_max_mismatch,
                request.rule_cleavage_site,
                request.rule_seed_mismatch,
                request.rule_consecutive,
                request.rule_seed_mfe,
                request.seed_mfe_threshold,
            )
        else:
            # Manual: compute p-value since it's just 1 target
            target_id, target_seq = parse_fasta(request.target_fasta)
            # v3.0: kullanıcı tercihi geçerli (animal override yok)
            effective_cleavage = request.strict_cleavage
            engine_result = find_targets(
                mirna_id, mirna_seq, target_id, target_seq,
                compute_pvalue=True,
                rule_set=request.rule_set,
                max_mismatches=request.max_mismatches,
                strict_cleavage=effective_cleavage,
                mfe_threshold=request.mfe_threshold,
                rule_max_mismatch=request.rule_max_mismatch,
                rule_cleavage_site=request.rule_cleavage_site,
                rule_seed_mismatch=request.rule_seed_mismatch,
                rule_consecutive=request.rule_consecutive,
                rule_seed_mfe=request.rule_seed_mfe,
                seed_mfe_threshold=request.seed_mfe_threshold,
            )
            results = [{
                "prediction": engine_result,
                "biological_context": {
                    "description":  "Manual Target Sequence",
                    "biotype":      "N/A",
                    "database":     "Manual Input",
                    "organism":     "Manual",
                    "organism_full": "Manual Input (user-provided sequence)",
                    "genome_build": "N/A",
                    "taxon_id":     0,
                }
            }]

        return {"status": "success", "data": results}

    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# ─── ENDPOINT: /api/batch-analyze ─────────────────────────────────────────────

@app.post("/api/batch-analyze")
async def batch_analysis(request: BatchRequest):
    """
    Analyze up to 20 miRNAs in parallel against the full database.
    Returns a list of results keyed by miRNA ID.
    """
    if len(request.mirnas) > 20:
        raise HTTPException(status_code=400,
                            detail="Maximum 20 miRNAs per batch request.")
    try:
        tasks = []
        for fasta in request.mirnas:
            mirna_id, mirna_seq = parse_fasta(fasta)
            tasks.append(asyncio.to_thread(
                run_database_search,
                mirna_id, mirna_seq,
                request.mfe_threshold,
                request.max_mismatches,
                request.strict_cleavage,
                request.rule_set,
                request.organism,
                request.rule_max_mismatch,
                request.rule_cleavage_site,
                request.rule_seed_mismatch,
                request.rule_consecutive,
                request.rule_seed_mfe,
                request.seed_mfe_threshold,
            ))

        all_results = await asyncio.gather(*tasks, return_exceptions=True)

        batch_out = []
        for fasta, res in zip(request.mirnas, all_results):
            mirna_id, _ = parse_fasta(fasta)
            if isinstance(res, Exception):
                batch_out.append({"mirna_id": mirna_id,
                                  "error": str(res), "data": []})
            else:
                batch_out.append({"mirna_id": mirna_id, "data": res,
                                  "match_count": len(res)})

        return {"status": "success", "batch_results": batch_out}

    except Exception as e:
        logger.error(f"Batch analysis failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# ─── ENDPOINT: /api/mirna-lookup/{name} ───────────────────────────────────────

def _fetch_mirbase_sync(mirna_name: str) -> dict:
    """
    Fetch miRNA sequence from miRBase REST API (synchronous, runs in thread).
    Falls back to NCBI eUtils if miRBase returns no result.
    """
    # 1) Try miRBase REST API
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

    # 2) Fallback: NCBI eUtils
    try:
        term    = urllib.parse.quote(f"{mirna_name}[Gene Name]")
        search  = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                   f"?db=gene&term={term}&retmode=json&retmax=1")
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
    """
    Auto-fetch miRNA sequence by name (e.g. hsa-miR-155-5p).
    Returns sequence ready to paste into the predictor.
    """
    result = await asyncio.to_thread(_fetch_mirbase_sync, mirna_name)
    if not result["found"]:
        raise HTTPException(status_code=404,
                            detail=f"miRNA '{mirna_name}' not found in miRBase or NCBI.")
    return result


# ─── ENDPOINT: /api/cache-stats ───────────────────────────────────────────────

@app.get("/api/cache-stats")
async def cache_stats():
    return {
        "cached_queries": len(_search_cache),
        "max_cache_size": _MAX_CACHE,
        "keys": list(_search_cache.keys()),
    }


# ─── ENDPOINT: /api/organisms ─────────────────────────────────────────────────

@app.get("/api/organisms")
async def list_organisms():
    """
    Returns all registered organism databases and their availability status.
    Use this endpoint to populate organism selectors in the frontend.
    """
    result = []
    for org_id, info in ORGANISM_DATABASES.items():
        db_exists = os.path.exists(info["db_path"])
        result.append({
            "id":           org_id,
            "display_name": info["display_name"],
            "common_name":  info["common_name"],
            "taxon_id":     info.get("taxon_id", 0),
            "genome_build": info.get("genome_build", "N/A"),
            "notes":        info.get("notes", ""),
            "available":    db_exists,
            "db_size_mb":   round(os.path.getsize(info["db_path"]) / 1e6, 1)
                            if db_exists else 0,
        })
    return {"organisms": result, "default": DEFAULT_ORGANISM}


# ─── HEALTH CHECK ─────────────────────────────────────────────────────────────

@app.get("/api/health")
async def health():
    organisms_status = {
        org_id: os.path.exists(info["db_path"])
        for org_id, info in ORGANISM_DATABASES.items()
    }
    return {
        "status":           "ok",
        "engine":           "v2.1",
        "organisms":        organisms_status,
        "default_organism": DEFAULT_ORGANISM,
        "db_exists":        os.path.exists(DB_PATH),
        "db_size_mb":       round(os.path.getsize(DB_PATH) / 1e6, 1)
                            if os.path.exists(DB_PATH) else 0,
    }