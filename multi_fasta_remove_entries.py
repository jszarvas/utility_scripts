#!/usr/bin/env python3

import sys
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(
    description='Remove entries from multi-fasta file by word filters')
parser.add_argument(
    '-i',
    dest="filename",
    default=None,
    help='Input multi-fasta')
parser.add_argument(
    '-o',
    dest="ofile",
    default=None,
    help='Output multi-fasta')
parser.add_argument(
    '-w',
    dest="filter",
    default=None,
    help='Comma separated words to filter on')
parser.add_argument(
    '--partition',
    dest="part",
    action="store_true",
    help='Create partitions file')
parser.add_argument(
    '--concat',
    dest="concat",
    action="store_true",
    help='Concatenate the entries')
args = parser.parse_args()

filter_words = args.filter.split(",")
entries_kept = []
lens = []
sample_id = ""
seq_iter = SeqIO.parse(args.filename, "fasta")
for rec in seq_iter:
    keep = True
    if rec.id in filter_words:
        keep = False
    for word in filter_words:
        if rec.description.find(word) != -1:
            keep = False
    if keep:
        entries_kept.append(rec)
        lens.append(len(rec.seq))
    if not sample_id:
        sample_id = rec.id

with open(args.ofile, "w") as of:
    SeqIO.write(entries_kept, of, "fasta")

if args.concat:
    concat = []
    for entry in entries_kept:
        concat.append(str(entry.seq))
    full_seq = "".join(concat)
    with open("{}.flat".format(args.ofile.rsplit(".",1)[0]), "w") as fop:
        print(">{}".format(sample_id), file=fop)
        for i in range(0, len(full_seq), 60):
            print(full_seq[i:i+60], file=fop)

if args.part:
    cummulated_len = 0
    with open("{}.partitions".format(args.ofile.rsplit(".",1)[0]), "w") as pop:
        for i, part_len in enumerate(lens):
            end_len = cummulated_len + part_len
            print("DNA, part{0} = {1}-{2}".format(i+1, cummulated_len + 1, end_len), file=pop)
            cummulated_len = end_len
