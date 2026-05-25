import csv, sys

csv_path = "database.csv"
print("Reading CSV headers...")

with open(csv_path, "r", encoding="utf-8", errors="replace") as f:
    reader = csv.DictReader(f)
    print("Columns:", reader.fieldnames)
    for i, row in enumerate(reader):
        if i >= 5:
            break
        print(f"\n--- Row {i+1} ---")
        for k in reader.fieldnames:
            val = str(row.get(k, ""))[:100]
            if val.strip():
                print(f"  {k}: {val}")
