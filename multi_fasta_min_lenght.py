#!/usr/bin/env python3

import sys
import os
from Bio import SeqIO

if (len(sys.argv) == 1):
    print("program <multi fasta> <output fasta> <min length>")
    sys.exit(1)
else:
    filename = sys.argv[1]
    min_length = int(sys.argv[3])

seq_iter = SeqIO.parse(filename, "fasta")

longer_records = []
for rec in seq_iter:
    if len(rec) >= min_length:
        longer_records.append(rec)

out_fastafilename = sys.argv[2]
with open(out_fastafilename, "w") as of:
    SeqIO.write(longer_records, of, "fasta")
