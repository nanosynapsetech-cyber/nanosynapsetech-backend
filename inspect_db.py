import sqlite3, os, sys

# Try sqlite3 version info
print("Python sqlite3 version:", sqlite3.sqlite_version)

db = "database.db"
print("DB size:", round(os.path.getsize(db)/1e6, 1), "MB")

# Connect without loading FTS5 catalog
conn = sqlite3.connect(db)
cur = conn.cursor()

# Try reading just the genes table directly  
try:
    cur.execute("PRAGMA table_info(genes)")
    cols_info = cur.fetchall()
    print("genes columns:", [c[1] for c in cols_info])
except Exception as e:
    print("PRAGMA error:", e)

try:
    cur.execute("SELECT Gene_ID, Description, Biotype FROM genes LIMIT 5")
    rows = cur.fetchall()
    for r in rows:
        print("Gene_ID:", str(r[0])[:60])
        print("Desc:", str(r[1])[:80])
        print("Biotype:", r[2])
        print("---")
except Exception as e:
    print("SELECT error:", e)

conn.close()
