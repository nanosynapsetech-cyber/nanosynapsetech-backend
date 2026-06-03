import sqlite3, os

db = "database.db"
if not os.path.exists(db):
    print("DB not found!")
    exit(1)

conn = sqlite3.connect(db)
cur = conn.cursor()

cur.execute("SELECT name, type FROM sqlite_master WHERE type IN ('table','index') ORDER BY type, name")
items = cur.fetchall()
print("=== Tables & Indexes ===")
for row in items:
    print(f"  [{row[1]}] {row[0]}")

print("\n=== Row Counts ===")
for row in items:
    if row[1] == "table" and not row[0].startswith("sqlite_"):
        try:
            cur.execute(f"SELECT COUNT(*) FROM [{row[0]}]")
            count = cur.fetchone()[0]
            print(f"  {row[0]}: {count:,} rows")
        except Exception as e:
            print(f"  {row[0]}: ERROR - {e}")

print("\n=== Gene Table Schema ===")
cur.execute("PRAGMA table_info(genes)")
for col in cur.fetchall():
    print(f"  {col[1]} ({col[2]})")

print("\n=== Sample Row ===")
cur.execute("SELECT Gene_ID, Description, Biotype, length(Sequence) as seq_len FROM genes LIMIT 3")
for row in cur.fetchall():
    print(f"  {row}")

print("\n=== DB File Size ===")
size_mb = os.path.getsize(db) / 1e6
print(f"  database.db: {size_mb:.1f} MB")

conn.close()
print("\nScan complete.")
