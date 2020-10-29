#!/usr/bin/env python3

import sys
import os
import argparse
from copy import deepcopy
from Bio import SeqIO

if (len(sys.argv) == 1):
    print("program <multi fasta> <opt N-threshold for print>")
    sys.exit(1)
else:
    filename = sys.argv[1]

parser = argparse.ArgumentParser(
    description='Caclulates pairwise kmer identity between input sequences in multifasta, or list of fastas, lower triangle is template, upper is query id')
parser.add_argument(
    '-m',
    dest="mfsa",
    default=None,
    help='Input multifasta')
parser.add_argument(
    '-t',
    dest="filter",
    default=100.0,
    type=float,
    help='Optional upper N-threshold for printing stats')
parser.add_argument(
    '-a',
    dest="ambig",
    action="store_true",
    help='Count ambigous IUPAC bases too')
args = parser.parse_args()

# init nuc dicts
nucs = {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0}
ambigous = {"M": 0, "R": 0, "W": 0 , "S": 0, "Y": 0, "K": 0, "B": 0 , "V": 0, "D": 0}
iupac = ["A", "T", "C", "G", "N", "M", "R", "W" , "S", "Y", "K", "B" , "V", "D"]
b = iupac.index("M")

if args.ambig:
    nucs.update(ambigous)
    b = len(iupac)
init_nucs = deepcopy(nucs)

# print header
print("# description", sep="\t", end="\t")
print("\t".join(iupac[:b]), end="\t")
print("Total")

# iterate through Seq records
seq_iter = SeqIO.parse(args.mfsa, "fasta")
for rec in seq_iter:
    for char in rec.seq:
        try:
            nucs[char.upper()] += 1
        except KeyError:
            pass

    total = sum(nucs.values())
    nucs_perc = [str(round(nucs.get(key)/float(total) * 100 , 3)) for key in iupac[:b]]
    if float(nucs_perc[-1]) < args.filter:
        print(rec.description, "\t".join(nucs_perc), total, sep="\t")
    nucs.update(init_nucs)
