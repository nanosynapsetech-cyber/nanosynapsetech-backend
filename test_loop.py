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
)
cursor.execute(query, (f'"{seed_rc}"',))
columns = [col[0] for col in cursor.description]

print("Looping rows...")
start = time.time()
processed = 0
results = []
while True:
    row = cursor.fetchone()
    if not row: break
    row_dict    = dict(zip(columns, row))
    target_id   = str(row_dict.get("Gene_ID", "Unknown"))
    target_seq  = str(row_dict.get("Sequence", "")).strip().upper()
    
    if len(target_seq) <= 10 or "N" in target_seq: continue
    
    match_idx = target_seq.find(seed_rc)
    start_idx = max(0, match_idx - 50)
    end_idx   = min(len(target_seq), match_idx + 50)
    short_target = target_seq[start_idx:end_idx]
    
    t_start = time.time()
    res = find_targets(mirna_id, mirna_seq, target_id, short_target)
    if res["Status"] == "PASS" and res.get("Similarity_Percent", 0) >= 65.0:
        results.append(res)
    processed += 1
    
    if processed % 100 == 0:
        print(f"Processed {processed} rows in {time.time()-start:.2f}s, found {len(results)} matches")
        
    if len(results) >= 15:
        break
        
print(f"Total time: {time.time()-start:.2f}s, Processed: {processed}")
