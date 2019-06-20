#!/usr/bin/env python3

import sys
import os
from Bio import SeqIO

if (len(sys.argv) == 1):
    print("program <name/name>.fsa <entrez len range>")
    sys.exit(1)
else:
    filename = os.path.join("/home/s151038/project/phyloKVIT/kvit_db", sys.argv[1], "{}.fsa".format(sys.argv[1]))
    # 6390[SLEN]:9790[SLEN]
    froml, tol = [int(x[:-6]) for x in sys.argv[2].split(":")]

seq_iter = SeqIO.parse(filename, "fasta")

failed = 0
okd = 0
for rec in seq_iter:
    if len(rec) < froml or len(rec) > tol:
        print(rec.id, len(rec), file=sys.stderr)
        failed += 1
    else:
        okd += 1
print("Final tally:", failed, okd, failed+okd)
