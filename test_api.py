import requests
import time
url = 'http://localhost:8000/api/analyze'
payload = {
    'mirna_fasta': '>hsa-miR-155-5p\nUUAAUGCUAAUCGUGAUAGGGGU',
    'search_mode': 'automatic',
    'mfe_threshold': -15.0,
    'max_mismatches': 4,
    'strict_cleavage': True
}
start = time.time()
try:
    res = requests.post(url, json=payload)
    data = res.json()
    print(f"Time: {time.time()-start:.2f}s")
    print(f"Status: {res.status_code}")
    if "data" in data:
        print(f"Match Count: {len(data['data'])}")
    else:
        print(f"Error: {data}")
except Exception as e:
    print(e)
