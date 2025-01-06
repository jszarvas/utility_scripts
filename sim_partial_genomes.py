#!/usr/bin/env python3

import os
import sys
from random import randint
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

MIN_LEN = 10000
MAX_LEN = 30000

def get_random_len(min_len = MIN_LEN, max_len = MAX_LEN):
    return randint(min_len, max_len)

if len(sys.argv) < 4:
    print("Usage: program <genome fasta> <relative length [0,1]> <output file prefix> <opt: random chunks>")
    sys.exit()

desired_coverage = float(sys.argv[2])

chunk_it_more = False
if len(sys.argv) == 5:
    chunk_it_more = True

records = list(SeqIO.parse(sys.argv[1], "fasta"))
total_orig_len = len(records)
new_seqrecords = []
for rec in records:
    start = 0
    new_total_len = 0
    # simplistic way: each chr will be reduced with the same ratio
    orig_len = len(rec)
    target_len = int(orig_len * desired_coverage)
    if not chunk_it_more:
        # for metagenomic assembly inputs
        new_seqrecords.append(SeqRecord(rec.seq[start:target_len+1],
                                            id=f"{rec.id}_cov_{desired_coverage:.2f}",
                                            description=f"shortened contig"))
    else:
        # chop it into smaller chunks, recommended for complete chromosomal inputs
        missing_len = orig_len - target_len
        target_contig_no = int(target_len / ((MIN_LEN+MAX_LEN) / 2.5))
        skipping_len_avrg = int(missing_len / (target_contig_no - 1))
        for contig_no in range(target_contig_no):
            contig_len = get_random_len(min_len = min(MIN_LEN, target_len - new_total_len), max_len = min(MAX_LEN, target_len - new_total_len))
            new_seqrecords.append(SeqRecord(rec.seq[start:start + contig_len],
                                            id=f"{rec.id}_{len(rec.seq[start:start + contig_len])}",
                                            description=f"simulated contig"))
            start += contig_len
            new_total_len += contig_len
            skip_len = get_random_len(min_len = skipping_len_avrg, max_len = int(skipping_len_avrg * 1.5))
            start += skip_len
            if start >= orig_len or new_total_len >= target_len:
                break

with open(sys.argv[3], "w") as op:
    SeqIO.write(new_seqrecords, op, "fasta")
