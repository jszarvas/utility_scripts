#!/usr/bin/env python

import sys, os, time
import sqlite3
from ete3 import Tree
import subprocess


if len(sys.argv) < 2:
    sys.exit('Put metadata onto given tree.\n Usage: <program> <db/tsv path> <newick file> <output dir>')

db_source = sys.argv[1]
nw_filename = sys.argv[2]
outdir = sys.argv[3]

if not os.path.exists(db_source):
    sys.exit("Data source doesn't exists at {}".format(db_source))
if not os.path.exists(nw_filename):
    sys.exit("Newick path not correct.")

db_path = ""
conn = None
if db_source.rsplit(".",1)[-1] not in ["db", "sqlite"]:
    db_path = "{}.tmp.sqlite".format(db_source.rsplit(".", 1)[0])
    if not os.path.exists(db_path):
        inputstr = ".mode tabs\n.import {} metadata".format(db_source)
        proc = subprocess.Popen(["sqlite3", "{}".format(db_path)], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        proc.stdin.write(inputstr.encode())
        (stdoutdata, stderrdata) = proc.communicate()
        if proc.returncode:
            sys.exit("Sqlite import failed.")
        conn = sqlite3.connect(db_path)

else:
    conn = sqlite3.connect(db_source)

conn.commit()
cur = conn.cursor()

try:
    tree = Tree(nw_filename)
except NewickError as e:
    sys.ex("Couldn't open {0}".format(nw_filename))

for node in tree.traverse():
    if node.name and node.name != "template":
        if node.name.startswith('*'):
            cur.execute("select * from metadata where sra_id=?", (node.name[1:],))
            rec = cur.fetchone()
            if rec:
                node.name = '*{}'.format(" ".join([x.replace(":", ".") for x in rec if x is not None]))
        else:
            cur.execute("select * from metadata where sra_id=?", (node.name,))
            rec = cur.fetchone()
            if rec:
                node.name = " ".join([x.replace(":", ".") for x in rec if x is not None])

conn.close()
if db_path:
    os.unlink(db_path)

try:
    filename = os.path.join(outdir, "{}.m.nw".format(os.path.split(nw_filename)[-1].split(".")[0]))
    outfile = open(filename, "w")
    outfile.write(tree.write(format=0))
    outfile.close()
except IOError as e:
    sys.exit("Output path not correct.")

sys.exit()
