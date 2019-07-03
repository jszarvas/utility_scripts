#!/usr/bin/env python3

import sys
import os
from Bio import SeqIO

if (len(sys.argv) == 1):
    print("program <multi fasta> <output dir> <id split character>")
    sys.exit(1)
else:
    filename = sys.argv[1]

seq_iter = SeqIO.parse(filename, "fasta")
splitchar = sys.argv[3]

for rec in seq_iter:
    tmp = rec.id.split(splitchar)
    rec.id = tmp[0]
    if len(tmp) > 1:
        rec.description = splitchar.join(tmp[1:])
    ofn = os.path.join(sys.argv[2], "{0}.fa".format(rec.id))
    with open(ofn, "w") as of:
        SeqIO.write(rec, of, "fasta")
