# -*- coding: utf-8 -*-
"""
format_converter.py
-------------------
Converts a FASTA file (rna.fna) directly into a SQLite database (database.db).
Compatible with Python 2.7+.

Steps performed:
  1. Stream-parse the FASTA file record by record (no full file load into RAM).
  2. Insert rows in batches of CHUNK_SIZE for fast bulk-insert performance.
  3. Build B-tree indexes on Gene_ID and Biotype.
  4. Build an FTS5 trigram virtual table for fast sequence substring search.
"""

from __future__ import print_function  # Python 2/3 compatibility

import sqlite3
import os
import sys
import time

# ── Configuration ──────────────────────────────────────────────────────────────
INPUT_FASTA = "rna.fna"
OUTPUT_DB   = "database.db"
CHUNK_SIZE  = 5000   # Rows per transaction — balances speed vs. memory usage
# ───────────────────────────────────────────────────────────────────────────────


def parse_fasta(filepath):
    """
    Generator: yields (gene_id, sequence, description, biotype) tuples
    from an NCBI-style FASTA file, one record at a time.
    Sequences are returned in RNA format (T replaced with U).
    """
    current_id   = None
    current_desc = ""
    current_seq  = []

    # Python 2.7 open() with io.open for encoding support
    import io
    with io.open(filepath, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            if line.startswith(u">"):
                # Flush the previous record
                if current_id is not None:
                    full_seq = u"".join(current_seq).upper().replace(u"T", u"U")
                    yield (current_id, full_seq, current_desc, u"transcript")

                # Parse new header
                # NCBI format: >NM_000014.6 Homo sapiens alpha-2-macroglobulin...
                parts        = line[1:].split(u" ", 1)
                current_id   = parts[0]
                current_desc = parts[1] if len(parts) > 1 else u"Unknown Target"
                current_seq  = []
            else:
                current_seq.append(line)

    # Flush the last record
    if current_id is not None:
        full_seq = u"".join(current_seq).upper().replace(u"T", u"U")
        yield (current_id, full_seq, current_desc, u"transcript")


def create_schema(cursor):
    """Create the genes table."""
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS genes (
            Gene_ID     TEXT NOT NULL,
            Sequence    TEXT NOT NULL,
            Description TEXT,
            Biotype     TEXT
        )
    """)


def apply_performance_pragmas(cursor):
    """Apply SQLite PRAGMAs optimised for bulk-insert speed."""
    cursor.execute("PRAGMA journal_mode = OFF")
    cursor.execute("PRAGMA synchronous  = OFF")
    cursor.execute("PRAGMA cache_size   = 200000")
    cursor.execute("PRAGMA temp_store   = MEMORY")
    cursor.execute("PRAGMA locking_mode = EXCLUSIVE")


def build_indexes(cursor):
    """Create B-tree indexes for fast lookups."""
    print("  Building index on Gene_ID ...")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_gene_id ON genes(Gene_ID)")
    print("  Building index on Biotype ...")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_biotype ON genes(Biotype)")


def build_fts5_index(cursor):
    """
    Attempt to create an FTS5 virtual table with a trigram tokenizer.
    FTS5 with trigram support requires SQLite >= 3.35.0.
    Python 2.7's bundled SQLite is too old and will silently skip this step.
    Run optimize_db.py with a modern Python (3.8+) to build the FTS5 index.
    """
    print("  Checking FTS5 + trigram support ...")
    try:
        cursor.execute("""
            CREATE VIRTUAL TABLE IF NOT EXISTS genes_fts USING fts5(
                Sequence,
                content       = 'genes',
                content_rowid = 'rowid',
                tokenize      = 'trigram'
            )
        """)
        print("  Populating FTS5 index from genes table ...")
        cursor.execute("""
            INSERT INTO genes_fts(rowid, Sequence)
            SELECT rowid, Sequence
            FROM   genes
            WHERE  Biotype = 'transcript'
               OR  Biotype LIKE '%mrna%'
               OR  Biotype LIKE '%cdna%'
        """)
        print("  FTS5 trigram index built successfully.")
    except sqlite3.OperationalError as exc:
        print("  [SKIP] FTS5 trigram not supported by this SQLite build: {0}".format(exc))
        print("         Run 'python optimize_db.py' with Python 3.8+ to build the FTS5 index.")


def show_stats(cursor):
    """Print a quick summary of the newly created database."""
    cursor.execute("SELECT COUNT(*) FROM genes")
    total = cursor.fetchone()[0]
    print("\n  Total records in genes table : {0}".format(total))

    cursor.execute("SELECT DISTINCT Biotype FROM genes LIMIT 10")
    biotypes = [row[0] for row in cursor.fetchall()]
    print("  Unique biotype values        : {0}".format(", ".join(biotypes)))

    # FTS5 table may not exist on older SQLite — check first
    try:
        cursor.execute("SELECT COUNT(*) FROM genes_fts")
        fts_rows = cursor.fetchone()[0]
        print("  FTS5 indexed rows            : {0}".format(fts_rows))
    except sqlite3.OperationalError:
        print("  FTS5 table                   : not built (run optimize_db.py)")


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    if not os.path.exists(INPUT_FASTA):
        print("[ERROR] Input file '{0}' not found in the backend folder.".format(INPUT_FASTA))
        print("        Make sure rna.fna is placed inside the /backend directory.")
        sys.exit(1)

    # Remove a stale database so we always start clean
    if os.path.exists(OUTPUT_DB):
        os.remove(OUTPUT_DB)
        print("[INFO]  Existing '{0}' removed — starting fresh.".format(OUTPUT_DB))

    fasta_size_gb = os.path.getsize(INPUT_FASTA) / (1024.0 ** 3)
    print("[INFO]  Starting conversion: {0} -> {1}".format(INPUT_FASTA, OUTPUT_DB))
    print("        FASTA file size : {0:.2f} GB".format(fasta_size_gb))
    overall_start = time.time()

    conn   = sqlite3.connect(OUTPUT_DB)
    cursor = conn.cursor()

    apply_performance_pragmas(cursor)
    create_schema(cursor)

    # ── Stream FASTA and insert in chunks ──────────────────────────────────────
    print("\n[STEP 1/3] Parsing FASTA and inserting rows into SQLite ...")
    insert_sql = "INSERT INTO genes (Gene_ID, Sequence, Description, Biotype) VALUES (?, ?, ?, ?)"

    batch       = []
    total_rows  = 0
    last_report = time.time()

    for record in parse_fasta(INPUT_FASTA):
        batch.append(record)
        if len(batch) >= CHUNK_SIZE:
            cursor.executemany(insert_sql, batch)
            conn.commit()
            total_rows += len(batch)
            batch = []

            now = time.time()
            if now - last_report >= 5:
                elapsed = now - overall_start
                print("        {0:>10} records inserted  |  {1:.0f}s elapsed".format(
                    total_rows, elapsed))
                last_report = now

    # Flush remaining records
    if batch:
        cursor.executemany(insert_sql, batch)
        conn.commit()
        total_rows += len(batch)

    print("        Done — {0} records inserted in {1:.1f}s".format(
        total_rows, time.time() - overall_start))

    # ── Build indexes ──────────────────────────────────────────────────────────
    print("\n[STEP 2/3] Building B-tree indexes ...")
    build_indexes(cursor)
    conn.commit()

    # ── Build FTS5 index ───────────────────────────────────────────────────────
    print("\n[STEP 3/3] Building FTS5 trigram index for sequence search ...")
    build_fts5_index(cursor)
    conn.commit()

    # ── Summary ────────────────────────────────────────────────────────────────
    print("\n[STATS]")
    show_stats(cursor)

    conn.close()
    total_time = time.time() - overall_start
    db_size_mb = os.path.getsize(OUTPUT_DB) / (1024.0 ** 2)
    print("\n[DONE]  Conversion complete in {0:.1f}s".format(total_time))
    print("        Output : {0}  ({1:.0f} MB)".format(OUTPUT_DB, db_size_mb))


if __name__ == "__main__":
    main()