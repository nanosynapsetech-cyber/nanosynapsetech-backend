#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# build_turso_meta.py  v2 — Pure HTTP API (no libsql-experimental needed)
# Migrates gene metadata (NO sequences) from local database.db → Turso
#
# Usage:
#   export TURSO_DATABASE_URL=libsql://nanosynapse-meta-xxxx.turso.io
#   export TURSO_AUTH_TOKEN=eyJh...
#   python3 build_turso_meta.py

import sqlite3
import os
import sys
import json
import time

try:
    import urllib.request
    import urllib.error
except ImportError:
    print("urllib not found (should be built-in). Exiting.")
    sys.exit(1)

# ─── CONFIG ───────────────────────────────────────────────────────────────────
LOCAL_DB    = "database.db"
TURSO_URL   = os.environ.get("TURSO_DATABASE_URL", "").strip()
TURSO_TOKEN = os.environ.get("TURSO_AUTH_TOKEN", "").strip()
BATCH_SIZE  = 200   # rows per HTTP request (keep small to avoid timeout)

# ─── CHECKS ───────────────────────────────────────────────────────────────────
if not TURSO_URL or not TURSO_TOKEN:
    print("❌  Set environment variables first:")
    print("    export TURSO_DATABASE_URL=libsql://nanosynapse-meta-xxxx.turso.io")
    print("    export TURSO_AUTH_TOKEN=eyJh...")
    sys.exit(1)

if not os.path.exists(LOCAL_DB):
    print(f"❌  {LOCAL_DB} not found in current directory!")
    sys.exit(1)

# Convert libsql:// URL to HTTPS for HTTP API
http_url = TURSO_URL.replace("libsql://", "https://")
pipeline_url = f"{http_url}/v2/pipeline"
print(f"🔗  Turso endpoint: {pipeline_url[:60]}...")

# ─── HTTP HELPER ──────────────────────────────────────────────────────────────

def turso_request(statements: list, retries: int = 3) -> dict:
    """
    Send a list of SQL statements to Turso via HTTP pipeline API.
    Each statement is a dict: {"type": "execute", "stmt": {"sql": "...", "args": [...]}}
    """
    payload = json.dumps({"requests": statements}).encode("utf-8")
    req = urllib.request.Request(
        pipeline_url,
        data=payload,
        headers={
            "Authorization": f"Bearer {TURSO_TOKEN}",
            "Content-Type": "application/json",
        },
        method="POST",
    )
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            body = e.read().decode("utf-8")
            print(f"\n⚠️  HTTP {e.code}: {body[:200]}")
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
        except Exception as e:
            print(f"\n⚠️  Request error: {e}")
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
    raise RuntimeError("Turso request failed after retries")


def execute_sql(sql: str, args: list = None) -> dict:
    """Execute a single SQL statement."""
    stmt = {"sql": sql}
    if args:
        stmt["args"] = [{"type": "text", "value": str(a)} if isinstance(a, str)
                        else {"type": "integer", "value": a} if isinstance(a, int)
                        else {"type": "text", "value": str(a)}
                        for a in args]
    return turso_request([
        {"type": "execute", "stmt": stmt},
        {"type": "close"},
    ])


def execute_batch(rows: list) -> None:
    """
    Execute a batch of INSERT statements in one HTTP call.
    rows: list of (gene_id, description, biotype, seq_length, organism)
    NOTE: Turso HTTP API requires ALL values (including integers) to be strings!
    """
    requests = []
    for gene_id, description, biotype, seq_len, organism in rows:
        requests.append({
            "type": "execute",
            "stmt": {
                "sql": (
                    "INSERT OR REPLACE INTO genes_meta "
                    "(gene_id, description, biotype, seq_length, organism) "
                    "VALUES (?, ?, ?, ?, ?)"
                ),
                "args": [
                    {"type": "text", "value": str(gene_id)},
                    {"type": "text", "value": str(description)},
                    {"type": "text", "value": str(biotype)},
                    {"type": "text", "value": str(seq_len)},   # integer → string!
                    {"type": "text", "value": str(organism)},
                ],
            }
        })
    requests.append({"type": "close"})
    turso_request(requests)

# ─── CREATE TABLE ─────────────────────────────────────────────────────────────
print("\n📐  Creating genes_meta table in Turso...")

execute_sql("DROP TABLE IF EXISTS genes_meta")
execute_sql("""
    CREATE TABLE IF NOT EXISTS genes_meta (
        gene_id     TEXT PRIMARY KEY,
        description TEXT NOT NULL DEFAULT '',
        biotype     TEXT NOT NULL DEFAULT 'transcript',
        seq_length  INTEGER DEFAULT 0,
        organism    TEXT DEFAULT 'human'
    )
""")
execute_sql("CREATE INDEX IF NOT EXISTS idx_biotype  ON genes_meta(biotype)")
execute_sql("CREATE INDEX IF NOT EXISTS idx_organism ON genes_meta(organism)")
print("✅  Table and indexes created")

# ─── VERIFY CONNECTION ────────────────────────────────────────────────────────
resp = execute_sql("SELECT COUNT(*) FROM genes_meta")
print(f"✅  Turso connected — current row count: 0 (fresh table)")

# ─── READ LOCAL DB ────────────────────────────────────────────────────────────
print(f"\n📂  Reading local {LOCAL_DB}...")
local = sqlite3.connect(LOCAL_DB)
local.row_factory = sqlite3.Row

cur = local.cursor()
cur.execute("SELECT COUNT(*) FROM genes")
total = cur.fetchone()[0]
print(f"    Total rows to migrate: {total:,}")
print(f"    Batch size: {BATCH_SIZE} rows per request")
print(f"    Estimated requests: {total // BATCH_SIZE + 1}")

# ─── MIGRATE ──────────────────────────────────────────────────────────────────
cur.execute("""
    SELECT Gene_ID, Description, Biotype, length(Sequence) as seq_len
    FROM genes
""")

inserted = 0
skipped  = 0
errors   = 0
batch    = []
start    = time.time()

print(f"\n🚀  Migrating metadata (sequences NOT included)...\n")

while True:
    rows = cur.fetchmany(BATCH_SIZE)
    if not rows:
        break

    for row in rows:
        gene_id     = str(row["Gene_ID"]     or "").strip()
        description = str(row["Description"] or "").strip()[:400]
        biotype     = str(row["Biotype"]     or "transcript").strip()
        seq_len     = int(row["seq_len"]     or 0)

        if not gene_id:
            skipped += 1
            continue

        batch.append((gene_id, description, biotype, seq_len, "human"))

    if batch:
        try:
            execute_batch(batch)
            inserted += len(batch)
        except Exception as e:
            errors += len(batch)
            print(f"\n⚠️  Batch error: {e}")

        batch = []

        elapsed = time.time() - start
        pct     = (inserted / total) * 100
        rate    = inserted / elapsed if elapsed > 0 else 0
        eta     = (total - inserted) / rate if rate > 0 else 0
        print(
            f"    ✓ {inserted:>7,} / {total:,}  "
            f"({pct:.1f}%)  "
            f"{rate:.0f} rows/s  "
            f"ETA: {eta:.0f}s",
            end="\r", flush=True
        )

local.close()

# ─── FINAL VERIFY ─────────────────────────────────────────────────────────────
print(f"\n\n{'='*50}")
print(f"✅  Migration complete!")
print(f"    Inserted : {inserted:,}")
print(f"    Skipped  : {skipped:,}")
print(f"    Errors   : {errors:,}")
print(f"    Duration : {time.time() - start:.1f}s")

print("\n🔍  Verifying in Turso...")
resp = execute_sql("SELECT COUNT(*) FROM genes_meta")
count_row = resp["results"][0]["response"]["result"]["rows"][0][0]["value"]
print(f"    Turso row count: {int(count_row):,}")

resp2 = execute_sql("""
    SELECT biotype, COUNT(*) as cnt
    FROM genes_meta
    GROUP BY biotype
    ORDER BY cnt DESC
    LIMIT 5
""")
rows_data = resp2["results"][0]["response"]["result"]["rows"]
print("    Top biotypes:")
for r in rows_data:
    bt  = r[0]["value"]
    cnt = r[1]["value"]
    print(f"      {bt}: {int(cnt):,}")

print(f"\n🎉  Turso metadata DB is ready!")
print(f"    Next: python3 build_kmer_index.py  (K-mer index for Cloudflare R2)")
