#!/usr/bin/env python3

import os
import sys
from Bio import SeqIO

if len(sys.argv) < 3:
    print("Usage: program <reference fasta> <insert start - fasta id tsv> <gene fasta> <output file>")
    sys.exit()

# insert positions
gene_posi = {}
gene_insidx = []
with open(sys.argv[2], "r") as fp:
    for line in fp:
        if line.startswith("#"):
            continue
        tmp = line.strip().split()
        start = int(tmp[0])
        gene_insidx.append(start)
        gene_posi[start] = tmp[1]

# genes to be added
gins_records = SeqIO.to_dict(SeqIO.parse(sys.argv[3], "fasta"))

cummprevlen = 0
records = list(SeqIO.parse(sys.argv[1], "fasta"))
for rec in records:
    chunks = None
    start = 0
    centry_id = rec.id
    for ins_p in gene_insidx:
        # add the genes meant for this scaffold
        if ins_p > cummprevlen:
            if ins_p <= (cummprevlen + len(rec)):
                n_shift = 0
                # avoid putting a gene in the scaffold links
                if rec.seq[ins_p - 1- cummprevlen] == "N":
                    while rec.seq[ins_p + n_shift- cummprevlen] == "N":
                        n_shift += 1 
                # start of the insertion is the end of the previous chunk of genome
                if chunks is not None:
                    chunks = chunks + rec.seq[start:ins_p + n_shift - cummprevlen]
                else:
                    chunks = rec.seq[start:ins_p + n_shift - cummprevlen]
                # adding the gene as chunk accessed via its entry_id
                chunks = chunks + gins_records[gene_posi.get(ins_p)]
                # the start of the next chunk
                start = ins_p + n_shift - cummprevlen
            else:
                break
    # last chunk after the last deletion
    chunks = chunks + rec.seq[start:len(rec)]
    cummprevlen += len(rec)
    # replace current fasta entry with the perforated one
    rec.seq = chunks.seq
    rec.id = f"{centry_id}_{cummprevlen}"

print(records)

with open(sys.argv[4], "w") as op:
    SeqIO.write(records, op, "fasta")
