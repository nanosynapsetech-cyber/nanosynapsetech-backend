import sqlite3
import os

db = r'database.db'
print(f'DB size: {os.path.getsize(db)} bytes')
conn = sqlite3.connect(db)
cur = conn.cursor()
cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
print('Tables:', cur.fetchall())
try:
    cur.execute("SELECT * FROM genes LIMIT 2")
    cols = [d[0] for d in cur.description]
    print('Columns:', cols)
    rows = cur.fetchall()
    for row in rows:
        row_dict = dict(zip(cols, row))
        print('Biotype:', row_dict.get('Biotype'))
        print('Gene_ID:', row_dict.get('Gene_ID'))
        print('Seq len:', len(str(row_dict.get('Sequence', ''))))
except Exception as e:
    print('Error querying genes table:', e)
conn.close()
