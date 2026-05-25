from main import run_database_search
import time

mirna_id = "hsa-miR-155-5p"
mirna_seq = "UUAAUGCUAAUCGUGAUAGGGGU"
mfe_threshold = -15.0
max_mismatches = 4
strict_cleavage = True

start = time.time()
results = run_database_search(mirna_id, mirna_seq, mfe_threshold, max_mismatches, strict_cleavage)
end = time.time()

print(f"Time: {end-start:.2f}s")
print(f"Matches found: {len(results)}")
