# -*- coding: utf-8 -*-
"""
build_fts5.py
-------------
Adds FTS5 trigram virtual table to an existing database.db.
Requires Python 3.8+ with SQLite 3.35.0+ (FTS5 trigram support).
Run this once after format_converter.py has populated the database.
"""
import sqlite3
import os
import time

DB_PATH = "database.db"


def main():
    if not os.path.exists(DB_PATH):
        print("[ERROR] database.db not found. Run format_converter.py first.")
        return

    size_mb = os.path.getsize(DB_PATH) / (1024 ** 2)
    print("=" * 55)
    print("  FTS5 TRIGRAM INDEX BUILDER")
    print("=" * 55)
    print(f"  SQLite version : {sqlite3.sqlite_version}")
    print(f"  database.db    : {size_mb:.0f} MB")

    conn   = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # Performance PRAGMAs for index building
    cursor.execute("PRAGMA journal_mode = WAL")
    cursor.execute("PRAGMA synchronous  = OFF")
    cursor.execute("PRAGMA cache_size   = 200000")
    cursor.execute("PRAGMA temp_store   = MEMORY")

    # Check row count
    cursor.execute("SELECT COUNT(*) FROM genes")
    total = cursor.fetchone()[0]
    print(f"  genes rows     : {total:,}")

    # Drop old FTS5 table if exists (clean rebuild)
    print("\n[STEP 1/2] Dropping old FTS5 table if it exists...")
    cursor.execute("DROP TABLE IF EXISTS genes_fts")
    conn.commit()

    # Create FTS5 virtual table with trigram tokenizer
    print("[STEP 2/2] Creating FTS5 trigram virtual table...")
    cursor.execute("""
        CREATE VIRTUAL TABLE genes_fts USING fts5(
            Sequence,
            content       = 'genes',
            content_rowid = 'rowid',
            tokenize      = 'trigram'
        )
    """)
    print("           Populating index — this may take several minutes...")

    start = time.time()
    cursor.execute("""
        INSERT INTO genes_fts(rowid, Sequence)
        SELECT rowid, Sequence
        FROM   genes
        WHERE  Biotype = 'transcript'
           OR  Biotype LIKE '%mrna%'
           OR  Biotype LIKE '%cdna%'
    """)
    conn.commit()
    elapsed = time.time() - start

    cursor.execute("SELECT COUNT(*) FROM genes_fts")
    fts_count = cursor.fetchone()[0]

    conn.close()

    new_size_mb = os.path.getsize(DB_PATH) / (1024 ** 2)
    print("\n" + "=" * 55)
    print("  DONE!")
    print(f"  FTS5 indexed rows : {fts_count:,}")
    print(f"  Time taken        : {elapsed:.1f}s")
    print(f"  Final DB size     : {new_size_mb:.0f} MB")
    print("=" * 55)
    print("\n  Sequence search (seed-match) is now TURBO FAST.")
    print("  The prediction engine is fully operational.")


if __name__ == "__main__":
    main()
