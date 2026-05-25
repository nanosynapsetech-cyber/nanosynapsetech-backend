import sqlite3
conn = sqlite3.connect(':memory:')
conn.execute('CREATE TABLE genes (id INTEGER PRIMARY KEY, seq TEXT)')
conn.execute("INSERT INTO genes (seq) VALUES ('AGCUUAACCGG'), ('CCGGAGCUUAA'), ('AAAAA')")
conn.execute("CREATE VIRTUAL TABLE genes_fts USING fts5(seq, content='genes', content_rowid='id', tokenize='trigram')")
conn.execute("INSERT INTO genes_fts(genes_fts) VALUES('rebuild')")
res = conn.execute("SELECT genes.id, genes.seq FROM genes JOIN genes_fts ON genes.id = genes_fts.rowid WHERE genes_fts.seq MATCH '\"AGCUUAA\"'").fetchall()
print(res)
