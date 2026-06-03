#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# build_kmer_index.py
# Reads database.db, extracts all 7-mers from each sequence,
# builds an inverted index: { kmer_rc -> [gene_id, ...] }
# Saves results to ./kmer_output/ as JSON files (ready for R2 upload)
# Also saves sequences to ./seq_output/ (one .txt per gene_id)
#
# Usage:
#   python3 build_kmer_index.py
#
# Output:
#   kmer_output/AAAAAAA.json  -> ["NM_000014.6", "NM_000015.3", ...]
#   seq_output/NM_000014.6.txt -> ATGCGT...  (full sequence)
#
# PERFORMANCE NOTE:
#   The entire kmer index is held in RAM (approx 1-2 GB for 186K genes).
#   This is MUCH faster than reading/writing 16K+ JSON files per batch.
#   If you are RAM-constrained, reduce BATCH_SIZE or use a chunked approach.

import sqlite3
import os
import sys
import json
import time
from collections import defaultdict

# --- CONFIG ------------------------------------------------------------------
LOCAL_DB    = "database.db"
KMER_K      = 7          # seed region size (miRNA pos 2-8)
BATCH_SIZE  = 5000       # rows to read at once from SQLite
KMER_DIR    = "kmer_output"
SEQ_DIR     = "seq_output"
CHECKPOINT  = "kmer_checkpoint.json"  # resume support
SAVE_EVERY  = 10000      # write checkpoint every N rows (no disk flush needed)

# --- CHECKS ------------------------------------------------------------------
if not os.path.exists(LOCAL_DB):
    print(f"ERROR: {LOCAL_DB} not found!")
    sys.exit(1)

os.makedirs(KMER_DIR, exist_ok=True)
os.makedirs(SEQ_DIR,  exist_ok=True)
print(f"Output dirs: {KMER_DIR}/ and {SEQ_DIR}/")

# --- REVERSE COMPLEMENT ------------------------------------------------------
RC_MAP = str.maketrans("AUGCT", "UACGA")   # RNA reverse complement

def rev_comp_rna(seq: str) -> str:
    """Reverse complement of an RNA sequence (A<->U, G<->C)."""
    return seq.translate(RC_MAP)[::-1]

def extract_kmers(seq: str, k: int = KMER_K) -> set:
    """
    Extract all unique k-mers from sequence.
    Returns reverse complements (what we search FOR in targets).
    Skips k-mers containing N or non-AUGC characters.
    """
    kmers = set()
    seq = seq.upper().replace("T", "U")  # Normalize to RNA
    valid = set("AUGC")
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if all(c in valid for c in kmer):
            kmers.add(rev_comp_rna(kmer))
    return kmers

# --- LOAD CHECKPOINT ---------------------------------------------------------
start_rowid   = 0
start_offset  = 0

# kmer_index lives entirely in RAM - much faster than per-batch disk I/O
kmer_index: dict[str, set] = defaultdict(set)

if os.path.exists(CHECKPOINT):
    print("Checkpoint found - resuming...")
    with open(CHECKPOINT) as f:
        cp = json.load(f)
        start_rowid  = cp.get("last_rowid", 0)
        start_offset = cp.get("processed", 0)

    # Reload existing kmer JSON files into RAM so we don't lose prior work
    existing_files = os.listdir(KMER_DIR)
    print(f"    Loading {len(existing_files):,} existing kmer files into RAM...")
    load_t = time.time()
    for fname in existing_files:
        kmer = fname[:-5]  # strip .json
        fpath = os.path.join(KMER_DIR, fname)
        with open(fpath) as f:
            kmer_index[kmer] = set(json.load(f))
    print(f"    Loaded {len(kmer_index):,} kmer entries in {time.time()-load_t:.1f}s")
    print(f"    Resuming from rowid>{start_rowid} | processed={start_offset:,}")

# --- CONNECT -----------------------------------------------------------------
print(f"\nOpening {LOCAL_DB}...")
conn = sqlite3.connect(LOCAL_DB)
conn.row_factory = sqlite3.Row

cur = conn.cursor()
cur.execute("SELECT COUNT(*) FROM genes")
total = cur.fetchone()[0]
print(f"    Total genes: {total:,}")
print(f"    K-mer size:  {KMER_K}-mer")
print(f"    Remaining:   {total - start_offset:,} genes to process")
print(f"    Output:      {KMER_DIR}/ + {SEQ_DIR}/\n")

# --- BUILD INDEX -------------------------------------------------------------
# Use rowid-based cursor instead of OFFSET for fast resume (O(1) seek vs O(N))
cur.execute("""
    SELECT rowid, Gene_ID, Sequence
    FROM genes
    WHERE rowid > ?
    ORDER BY rowid
""", (start_rowid,))

processed   = start_offset
last_rowid  = start_rowid
seq_saved   = 0
kmer_total  = 0
start_time  = time.time()

print("Building K-mer inverted index (all in RAM - fast mode)...\n")

while True:
    rows = cur.fetchmany(BATCH_SIZE)
    if not rows:
        break

    for row in rows:
        last_rowid = row["rowid"]
        gene_id    = str(row["Gene_ID"] or "").strip()
        seq        = str(row["Sequence"] or "").strip().upper().replace("T", "U")

        if not gene_id or len(seq) < KMER_K:
            processed += 1
            continue

        # -- Save sequence file ------------------------------------------------
        seq_path = os.path.join(SEQ_DIR, f"{gene_id}.txt")
        if not os.path.exists(seq_path):
            with open(seq_path, "w") as f:
                f.write(seq)
            seq_saved += 1

        # -- Extract k-mers and build inverted index ---------------------------
        kmers = extract_kmers(seq, KMER_K)
        for kmer in kmers:
            kmer_index[kmer].add(gene_id)   # set deduplicates automatically
        kmer_total += len(kmers)
        processed  += 1

    # -- Save checkpoint (no disk flush of kmer_index - stays in RAM) ----------
    with open(CHECKPOINT, "w") as f:
        json.dump({"processed": processed, "last_rowid": last_rowid}, f)

    # -- Progress --------------------------------------------------------------
    elapsed = time.time() - start_time
    pct     = (processed / total) * 100
    rate    = (processed - start_offset) / elapsed if elapsed > 0 else 0
    eta     = (total - processed) / rate if rate > 0 else 0
    print(
        f"    {processed:>7,} / {total:,}  "
        f"({pct:.1f}%)  "
        f"{rate:.0f} rows/s  "
        f"kmer keys: {len(kmer_index):,}  "
        f"ETA: {int(eta//60)}m{int(eta%60)}s",
        end="\r", flush=True
    )

conn.close()
print()  # newline after \r progress

# --- FLUSH KMER INDEX TO DISK (one final pass) --------------------------------
print(f"\nFlushing {len(kmer_index):,} kmer files to {KMER_DIR}/...")
flush_t = time.time()
for i, (kmer, gene_ids) in enumerate(kmer_index.items()):
    kmer_path = os.path.join(KMER_DIR, f"{kmer}.json")
    with open(kmer_path, "w") as f:
        json.dump(sorted(gene_ids), f)
    if i % 10000 == 0:
        print(f"    Flushed {i:,} / {len(kmer_index):,}...", end="\r", flush=True)

print(f"\n    Done! Flushed in {time.time()-flush_t:.1f}s")

# --- FINAL SUMMARY ------------------------------------------------------------
kmer_files = len(os.listdir(KMER_DIR))
seq_files  = len(os.listdir(SEQ_DIR))
elapsed    = time.time() - start_time

print(f"\n{'='*55}")
print("K-mer index build complete!")
print(f"    Sequences saved   : {seq_files:,}  -> {SEQ_DIR}/")
print(f"    Unique K-mer files: {kmer_files:,}  -> {KMER_DIR}/")
print(f"    Total K-mer hits  : {kmer_total:,}")
print(f"    Duration          : {elapsed/60:.1f} minutes")

# Estimate R2 upload size
kmer_size = sum(os.path.getsize(os.path.join(KMER_DIR, f))
                for f in os.listdir(KMER_DIR)) / 1e6
seq_size  = sum(os.path.getsize(os.path.join(SEQ_DIR, f))
                for f in os.listdir(SEQ_DIR)) / 1e6
print(f"\n    R2 upload estimate:")
print(f"      kmer_output/: {kmer_size:.1f} MB")
print(f"      seq_output/ : {seq_size:.1f} MB")
print(f"      Total       : {kmer_size + seq_size:.1f} MB")

# Remove checkpoint on success
if os.path.exists(CHECKPOINT):
    os.remove(CHECKPOINT)

print("\nReady to upload to Cloudflare R2!")
print("    Next: python3 upload_to_r2.py")
