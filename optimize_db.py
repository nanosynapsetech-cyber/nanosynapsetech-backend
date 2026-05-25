# -*- coding: utf-8 -*-
import sqlite3
import os
import sys

DB_PATH = "database.db"

def optimize():
    if not os.path.exists(DB_PATH):
        print(f"Error: {DB_PATH} not found!")
        sys.exit(1)

    print("Connecting to database...")
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # Try to proceed without changing journal_mode if database is locked
    # cursor.execute("PRAGMA journal_mode = WAL")
    # cursor.execute("PRAGMA synchronous = NORMAL")
    cursor.execute("PRAGMA cache_size = -100000") # ~100MB cache
    
    print("Checking for existing FTS table...")
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='genes_fts'")
    if cursor.fetchone():
        print("FTS table 'genes_fts' already exists. Dropping it to rebuild...")
        cursor.execute("DROP TABLE genes_fts")

    print("Creating FTS5 virtual table with trigram tokenizer...")
    try:
        cursor.execute("""
            CREATE VIRTUAL TABLE genes_fts USING fts5(
                Sequence, 
                content='genes', 
                content_rowid='rowid', 
                tokenize='trigram'
            )
        """)
    except sqlite3.OperationalError as e:
        print(f"Error creating FTS5 table: {e}")
        print("Make sure your Python SQLite library supports FTS5 and trigram tokenizer.")
        sys.exit(1)

    print("Populating FTS index from the genes table...")
    print("This may take several minutes depending on database size. Please wait...")
    
    # We only index genes that match our usual criteria to save space and time
    # (Biotype = 'transcript' OR Biotype LIKE '%mrna%' OR Biotype LIKE '%cdna%')
    cursor.execute("""
        INSERT INTO genes_fts(rowid, Sequence) 
        SELECT rowid, Sequence FROM genes 
        WHERE Biotype = 'transcript' 
           OR Biotype LIKE '%mrna%' 
           OR Biotype LIKE '%cdna%'
    """)
    
    conn.commit()
    print("FTS index built successfully!")
    
    print("Database size after indexing: {} MB".format(round(os.path.getsize(DB_PATH) / (1024 * 1024), 2)))
    conn.close()

if __name__ == "__main__":
    optimize()
