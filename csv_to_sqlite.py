# -*- coding: utf-8 -*-
"""
csv_to_sqlite.py
----------------
Adds FTS5 trigram index to an existing database.db.
If database.csv is present, it will rebuild database.db from scratch first.
Compatible with Python 2.7+.
"""
from __future__ import print_function
import sqlite3
import os

DB_PATH  = "database.db"
CSV_PATH = "database.csv"


def rebuild_from_csv(conn, cursor):
    """
    (Optional path) Rebuild genes table from database.csv using pandas.
    Only runs if database.csv exists.
    """
    try:
        import pandas as pd
    except ImportError:
        print("[ERROR] pandas is not installed. Run: pip install pandas")
        return False

    print("Reading {0} and writing into SQLite in chunks...".format(CSV_PATH))
    chunksize = 50000
    for i, chunk in enumerate(pd.read_csv(CSV_PATH, chunksize=chunksize)):
        chunk.to_sql("genes", conn,
                     if_exists="replace" if i == 0 else "append",
                     index=False)
        print("  {0} rows processed...".format((i + 1) * chunksize))

    cursor.execute("CREATE INDEX IF NOT EXISTS idx_gene_id ON genes(Gene_ID)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_biotype ON genes(Biotype)")
    conn.commit()
    print("Indexes created.")

    print("\nUnique Biotype values:")
    cursor.execute("SELECT DISTINCT Biotype FROM genes LIMIT 20")
    for row in cursor.fetchall():
        print("  ", row[0])

    return True


def build_fts5(cursor, conn):
    """
    Add FTS5 trigram virtual table to the existing genes table.
    Requires SQLite >= 3.35.0 (Python 3.8+ recommended).
    """
    print("\nCreating FTS5 virtual table with trigram tokenizer...")
    try:
        cursor.execute("""
            CREATE VIRTUAL TABLE IF NOT EXISTS genes_fts USING fts5(
                Sequence,
                content       = 'genes',
                content_rowid = 'rowid',
                tokenize      = 'trigram'
            )
        """)
        print("Populating FTS5 index from the genes table...")
        cursor.execute("""
            INSERT INTO genes_fts(rowid, Sequence)
            SELECT rowid, Sequence FROM genes
            WHERE  Biotype = 'transcript'
               OR  Biotype LIKE '%mrna%'
               OR  Biotype LIKE '%cdna%'
        """)
        conn.commit()
        cursor.execute("SELECT COUNT(*) FROM genes_fts")
        count = cursor.fetchone()[0]
        print("FTS5 index ready — {0} rows indexed.".format(count))
        return True
    except sqlite3.OperationalError as e:
        print("[SKIP] FTS5 not supported by this SQLite build: {0}".format(e))
        print("       Upgrade to Python 3.8+ for FTS5 trigram support.")
        return False


def convert():
    # ── Case 1: database.csv exists → rebuild DB from CSV ─────────────────────
    if os.path.exists(CSV_PATH):
        print("[INFO] Found {0} — rebuilding database from CSV.".format(CSV_PATH))
        if os.path.exists(DB_PATH):
            os.remove(DB_PATH)
            print("[INFO] Old {0} deleted.".format(DB_PATH))
        conn   = sqlite3.connect(DB_PATH)
        cursor = conn.cursor()
        cursor.execute("PRAGMA journal_mode = OFF")
        cursor.execute("PRAGMA synchronous  = OFF")
        cursor.execute("PRAGMA cache_size   = 100000")
        cursor.execute("PRAGMA temp_store   = MEMORY")
        rebuild_from_csv(conn, cursor)

    # ── Case 2: database.db already exists → just add FTS5 ───────────────────
    elif os.path.exists(DB_PATH):
        print("[INFO] database.csv not found.")
        print("[INFO] Using existing database.db ({0:.0f} MB) — adding FTS5 index only.".format(
            os.path.getsize(DB_PATH) / (1024.0 ** 2)))
        conn   = sqlite3.connect(DB_PATH)
        cursor = conn.cursor()
        cursor.execute("PRAGMA journal_mode = WAL")
        cursor.execute("PRAGMA synchronous  = OFF")
        cursor.execute("PRAGMA cache_size   = 200000")
        cursor.execute("PRAGMA temp_store   = MEMORY")

    # ── Case 3: nothing to work with ──────────────────────────────────────────
    else:
        print("[ERROR] Neither {0} nor {1} found.".format(CSV_PATH, DB_PATH))
        print("        Run format_converter.py first to create database.db.")
        return

    # Verify genes table exists
    cursor.execute("SELECT COUNT(*) FROM genes")
    row_count = cursor.fetchone()[0]
    print("[INFO] genes table has {0} rows.".format(row_count))

    # Build FTS5 on top of existing data
    build_fts5(cursor, conn)

    conn.close()
    print("\nDone! database.db is ready.")


if __name__ == "__main__":
    convert()
