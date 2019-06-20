#!/usr/bin/env python3

import sys
import os
from Bio import SeqIO

if (len(sys.argv) == 1):
    print("program <multi fasta> <output dir> <seq id1> <opt seq id2>")
    sys.exit(1)
else:
    filename = sys.argv[1]

seq_iter = SeqIO.parse(filename, "fasta")

for rec in seq_iter:
    if rec.id in sys.argv[3:]:
        ofn = os.path.join(sys.argv[2], "{0}_{1}.fa".format(os.path.basename(filename).split(".")[0], rec.id.replace(".", "_").replace(":","_")))
        rec.id = os.path.basename(filename).split(".")[0]
        with open(ofn, "w") as of:
            SeqIO.write(rec, of, "fasta") 

