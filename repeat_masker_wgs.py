#!/data/tools/anaconda3/bin/python3

import os
import sys
import argparse
from Bio import SeqIO
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def exiting(msg):
    print(msg, file=sys.stderr)
    sys.exit(1)


parser = argparse.ArgumentParser(description='Find repeating regions in a genome')
parser.add_argument("-g", "--genome", dest="gen_f", help="Input genome fasta")
parser.add_argument("-k", "--ksize", dest="ksize", default=21, type=int, help="K-mer size")
parser.add_argument("-f", "--feature", dest="gff", help="Optional feature file")
args = parser.parse_args()

if args.gen_f is None or not os.path.exists(args.gen_f):
    exiting("Input file not found")

k = args.ksize

kmers = {}
abs_pos = 0
for rec in SeqIO.parse(args.gen_f, "fasta"):
    for i in range(len(rec)-k+1):
        abs_pos += 1
        kmer = rec.seq[abs_pos:abs_pos+k]
        if kmer not in kmers:
            kmers[kmer] = [abs_pos]
        else:
            kmers[kmer].append(abs_pos)
    abs_pos += k - 1

print(len(rec), abs_pos)

genome = np.zeros(shape=(1,abs_pos), dtype=np.int8)
for kmer in kmers.keys():
    if len(kmers.get(kmer)) > 1:
        for pos in kmers.get(kmer):
            np.put(genome, range(pos-1,pos+k-1), 1)

plt.figure(figsize = (30,3))
plt.imshow(genome, cmap='Greys',  interpolation='none', aspect='auto')
plt.savefig(os.path.join(os.getcwd(), '{}.png'.format(os.path.basename(args.gen_f))))

if args.gff is not None:
    _, positions = genome.nonzero()
    features = []
    with open(args.gff, "r") as fp:
        for line in fp:
            # ctg123  .  exon  1300  1500  .  +  .  ID=exon00001
            if not line.startswith("#"):
                cols = line.split("\t")
                start = int(cols[3])
                end = int(cols[4])
                if cols[2] in ['exon', 'gene', 'region']:
                    # chr line in the beginning
                    continue
                if start > positions[-1]:
                    # past the repeats
                    break

                i = 0
                while i < np.shape(positions)[0] and positions[i] < end:
                    if positions[i] >= start:
                        if line not in features:
                            features.append(line)
                            i = np.shape(positions)[0]
                    i += 1

    print("".join(features))
