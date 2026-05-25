# -*- coding: utf-8 -*-
import sqlite3
import pandas as pd
import os

DB_PATH = "database.db"
CSV_PATH = "database.csv"

def convert():
    if not os.path.exists(CSV_PATH):
        print(f"{CSV_PATH} not found!")
        return

    # Delete existing (likely empty/corrupt) database
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)
        print(f"Old {DB_PATH} deleted.")

    print("Creating SQLite database...")
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    
    # Performance: disable WAL and set large cache for bulk inserts
    cursor.execute("PRAGMA journal_mode = OFF")
    cursor.execute("PRAGMA synchronous = OFF")
    cursor.execute("PRAGMA cache_size = 100000")
    cursor.execute("PRAGMA temp_store = MEMORY")
    
    chunksize = 50000
    print(f"Reading {CSV_PATH} and writing in chunks of {chunksize}...")
    
    for i, chunk in enumerate(pd.read_csv(CSV_PATH, chunksize=chunksize)):
        chunk.to_sql("genes", conn, if_exists="replace" if i == 0 else "append", index=False)
        print(f"  {(i+1) * chunksize} rows processed...")

    print("Creating indexes (this makes searches much faster)...")
    # Index on Gene_ID for lookups
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_gene_id ON genes(Gene_ID)")
    # Index on Biotype for filtering (transcript, protein_coding etc)
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_biotype ON genes(Biotype)")
    
    # Show the actual Biotype values in the DB so we know what to filter on
    print("\nUnique Biotype values in database:")
    cursor.execute("SELECT DISTINCT Biotype FROM genes LIMIT 20")
    for row in cursor.fetchall():
        print(" ", row[0])

    print("\nCreating FTS5 virtual table with trigram tokenizer...")
    cursor.execute("""
        CREATE VIRTUAL TABLE IF NOT EXISTS genes_fts USING fts5(
            Sequence, 
            content='genes', 
            content_rowid='rowid', 
            tokenize='trigram'
        )
    """)
    print("Populating FTS index from the genes table...")
    cursor.execute("""
        INSERT INTO genes_fts(rowid, Sequence) 
        SELECT rowid, Sequence FROM genes 
        WHERE Biotype = 'transcript' 
           OR Biotype LIKE '%mrna%' 
           OR Biotype LIKE '%cdna%'
    """)

    conn.commit()
    conn.close()
    print("\nConversion complete! database.db and FTS5 index are ready.")

if __name__ == "__main__":
    convert()
