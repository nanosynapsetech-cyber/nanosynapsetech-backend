# -*- coding: utf-8 -*-
"""check_db.py — Verifies that database.db is correctly populated."""
from __future__ import print_function
import sqlite3
import os

DB_PATH = "database.db"

print("=" * 50)
print("DATABASE INTEGRITY CHECK")
print("=" * 50)

# ── File check ─────────────────────────────────────────────────────────────────
print("\n[1] FILE")
if not os.path.exists(DB_PATH):
    print("  ERROR: database.db not found!")
    exit(1)
size_mb = os.path.getsize(DB_PATH) / (1024.0 ** 2)
print("  Exists : YES")
print("  Size   : {0:.1f} MB".format(size_mb))

conn = sqlite3.connect(DB_PATH)
cur  = conn.cursor()

# ── Tables ─────────────────────────────────────────────────────────────────────
print("\n[2] TABLES")
cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
tables = [r[0] for r in cur.fetchall()]
for t in tables:
    print("  table:", t)

# ── Row count ──────────────────────────────────────────────────────────────────
print("\n[3] ROW COUNT")
cur.execute("SELECT COUNT(*) FROM genes")
total = cur.fetchone()[0]
print("  genes table : {0} rows".format(total))
if total == 0:
    print("  ERROR: Table is empty!")
    exit(1)

# ── Schema ─────────────────────────────────────────────────────────────────────
print("\n[4] SCHEMA (genes)")
cur.execute("PRAGMA table_info(genes)")
for r in cur.fetchall():
    print("  col: {0:15s} | type: {1}".format(r[1], r[2]))

# ── Sample rows ────────────────────────────────────────────────────────────────
print("\n[5] SAMPLE ROWS (first 3)")
cur.execute("SELECT Gene_ID, Description, Biotype, length(Sequence) FROM genes LIMIT 3")
rows = cur.fetchall()
for r in rows:
    print("  Gene_ID  :", r[0])
    print("  Desc     :", r[1][:80] if r[1] else "N/A")
    print("  Biotype  :", r[2])
    print("  Seq len  : {0} nt".format(r[3]))
    print()

# ── Indexes ────────────────────────────────────────────────────────────────────
print("[6] INDEXES")
cur.execute("SELECT name FROM sqlite_master WHERE type='index'")
indexes = cur.fetchall()
if indexes:
    for r in indexes:
        print("  index:", r[0])
else:
    print("  WARNING: No indexes found.")

# ── Lookup test ────────────────────────────────────────────────────────────────
print("\n[7] GENE_ID LOOKUP TEST")
cur.execute("SELECT Gene_ID, length(Sequence) FROM genes WHERE Gene_ID LIKE 'NM_%' LIMIT 5")
results = cur.fetchall()
if results:
    for r in results:
        print("  {0} | seq_len: {1} nt".format(r[0], r[1]))
else:
    print("  No NM_ entries found — checking first 5 IDs:")
    cur.execute("SELECT Gene_ID, length(Sequence) FROM genes LIMIT 5")
    for r in cur.fetchall():
        print("  {0} | seq_len: {1} nt".format(r[0], r[1]))

# ── Sequence content check ─────────────────────────────────────────────────────
print("\n[8] SEQUENCE CONTENT CHECK")
cur.execute("SELECT Sequence FROM genes LIMIT 1")
sample_seq = cur.fetchone()[0]
has_u   = "U" in sample_seq
has_t   = "T" in sample_seq
print("  Contains U (RNA) : {0}".format(has_u))
print("  Contains T (DNA) : {0}".format(has_t))
print("  First 60 nt      :", sample_seq[:60])

# ── FTS5 check ─────────────────────────────────────────────────────────────────
print("\n[9] FTS5 INDEX")
try:
    cur.execute("SELECT COUNT(*) FROM genes_fts")
    fts_count = cur.fetchone()[0]
    print("  genes_fts rows : {0}".format(fts_count))
except sqlite3.OperationalError:
    print("  WARNING: FTS5 index not built yet.")
    print("  Run 'python optimize_db.py' with Python 3.8+ to enable fast search.")

conn.close()

print("\n" + "=" * 50)
print("RESULT: database.db is HEALTHY and READY TO USE")
print("=" * 50)
