import sqlite3
c = sqlite3.connect('database.db').cursor()
c.execute("SELECT COUNT(1) FROM genes WHERE Biotype LIKE '%mrna%' OR Biotype = 'transcript' OR Biotype LIKE '%cdna%'")
print(c.fetchone()[0])
