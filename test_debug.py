import time
import sqlite3
from mirna_engine import find_targets

mirna_id = "hsa-miR-155-5p"
mirna_seq = "UUAAUGCUAAUCGUGAUAGGGGU"
mirna_seed = mirna_seq[1:8]
rc_map     = {"A": "U", "U": "A", "G": "C", "C": "G"}
seed_rc    = "".join(rc_map.get(c, c) for c in reversed(mirna_seed))

conn = sqlite3.connect("database.db")
cursor = conn.cursor()

query = (
    "SELECT g.* FROM genes g "
    "JOIN genes_fts f ON g.rowid = f.rowid "
    "WHERE f.Sequence MATCH ? "
    "AND (g.Biotype = 'transcript' OR g.Biotype LIKE '%mrna%' "
    "     OR g.Biotype LIKE '%cdna%')"
)
print("Running query...")
start = time.time()
cursor.execute(query, (f'"{seed_rc}"',))
print(f"Query executed in {time.time()-start:.2f}s")

row = cursor.fetchone()
print(f"First fetch in {time.time()-start:.2f}s")
