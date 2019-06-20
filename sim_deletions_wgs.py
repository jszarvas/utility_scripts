#!/data/tools/anaconda3/bin/python3

import os
import sys
from Bio import SeqIO

if len(sys.argv) < 3:
    print("Usage: program <reference fasta> <deletion start - length tsv> <output file>")
    sys.exit()

deli_pos = []
with open(sys.argv[2], "r") as fp:
    for line in fp:
        if line.startswith("#"):
            continue
        tmp = line.strip().split()
        start = int(tmp[0])
        deli_pos.append([start, start + int(tmp[1])])

abs_pos = 0
chunks = []
start = 0
records = list(SeqIO.parse(sys.argv[1], "fasta"))
for rec in records:
    for dl in deli_pos:
        # start of the deletion is the end of the previous chunk of genome
        chunks.append(rec.seq[start:dl[0] - abs_pos])
        # the end + 1 of the deletion is going to be the start of the next chunk
        start = dl[1] - abs_pos
    # last chunk after the last deletion
    chunks.append(rec.seq[start:len(rec)])
    abs_pos += len(rec)
    # replace current fasta entry with the perforated one
    rec.seq = "".join(chunks)
    rec.id += "_dels"
    chunks = []

with open(sys.argv[3], "w") as op:
    SeqIO.write(records, op, "fasta")
