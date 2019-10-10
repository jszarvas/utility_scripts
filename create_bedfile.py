#!/usr/bin/env python3

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq

if len(sys.argv) < 2:
    sys.exit("program <reference> <opt output path> \n It takes the contig accession as chromosome ID")


records = list(SeqIO.parse(sys.argv[1], "fasta"))

contig = []
lens = []
for rec in records:
    print(rec.id)
    #print(vcf_i.seqnames)
    contig.append(rec.id)
    lens.append(len(rec.seq))

if contig:
    outfilename = "{}.bed".format(sys.argv[1].rsplit(".", 1)[0])
    if len(sys.argv) == 3:
        outfilename = sys.argv[2]
    with open(outfilename, "w") as pop:
        for i, part_len in enumerate(lens):
            print(contig[i], "0", part_len, sep="\t", file=pop)
