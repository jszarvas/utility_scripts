#!/usr/bin/env python3

import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

if (len(sys.argv) == 1):
    print("program <output file> <description> <mfa 1> <mfa 2> <mfa 3 etc>")
    sys.exit(1)
else:
    filenames = sys.argv[3:]

rec_ords = {}
lens = []
no_taxa = 0
# mfa file for each locus, containing multiple organisms
for filename in filenames:
    seqs = list(SeqIO.parse(filename, "fasta"))
    # append other loci to the 1st one, each to their own organism
    if not rec_ords:
        lens.append(len(seqs[0].seq))
        for rec in seqs:
            rec_ords[rec.id] = rec.seq
        no_taxa = len(rec_ords)
    else:
        if len(seqs) != no_taxa:
            print("{} missing taxa".format(filename))
        else:
            lens.append(len(seqs[0].seq))
            for rec in seqs:
                rec_ords[rec.id] += rec.seq

# put it out together
records = []
for recordid, concat_seq in rec_ords.items():
    concatted = SeqRecord(concat_seq)
    concatted.id = recordid
    concatted.description = sys.argv[2]
    records.append(concatted)
with open(sys.argv[1], "w") as of:
    SeqIO.write(records, of, "fasta")

# partitions file to go with it
cummulated_len = 0
with open("{}.partitions".format(sys.argv[1].rsplit(".",1)[0]), "w") as pop:
    for i, part_len in enumerate(lens):
        end_len = cummulated_len + part_len
        print("DNA, part{0} = {1}-{2}".format(i+1, cummulated_len + 1, end_len), file=pop)
        cummulated_len = end_len
