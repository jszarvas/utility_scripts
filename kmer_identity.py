#!/usr/bin/env python3

import sys
import os
import argparse
import hashlib
from Bio import SeqIO

def hash_seqid(seqid):
    hash = hashlib.sha1(seqid.encode("UTF-8")).hexdigest()[:16]
    fasta_ids[hash] = seqid
    return hash

def split_to_kmers(k, rec, hash):
    length = len(rec)
    for i in range(0,length - k):
        kmer_str = str(rec.seq[i:i+k])
        if kmer_str not in kmers:
            kmers[kmer_str] = [hash]
        else:
            if hash not in kmers[kmer_str]:
                kmers[kmer_str].append(hash)

parser = argparse.ArgumentParser(
    description='Caclulates pairwise kmer identity between input sequences in multifasta, or list of fastas, lower triangle is template, upper is query id')
parser.add_argument(
    '-m',
    dest="mfsa",
    default=None,
    help='Input multifasta')
parser.add_argument(
    '-l',
    dest="fsa_list",
    default=None,
    help='Input list of fastas')
parser.add_argument(
    '-o',
    dest="ofile",
    default=None,
    help='Output similance matrix')
parser.add_argument(
    '-k',
    dest="kmer",
    default=16,
    type=int,
    help='Kmer size')
args = parser.parse_args()


fasta_ids = {}
kmers = {}
if args.mfsa is not None and os.path.exists(args.mfsa):
    seq_iter = SeqIO.parse(args.mfsa, "fasta")
    for rec in seq_iter:
        hashid = hash_seqid(rec.id)
        split_to_kmers(args.kmer, rec, hashid)

elif args.fsa_list is not None and os.path.exists(args.fsa_list):
    with open(args.fsa_list, "r") as fp:
        for line in fp:
            parser = SeqIO.parse(line.strip(), "fasta")
            for rec in parser:
                hashid = hash_seqid(rec.id)
                split_to_kmers(args.kmer, rec, hashid)

#print stats
print(len(kmers))

matrix = []
seqlist = fasta_ids.keys()
# template
for i_seq in seqlist:
    matrix.append([])
    for j_seq in seqlist:
        if i_seq != j_seq:
            i_mers = 0.0
            common = 0
            for val in kmers.values():
                if i_seq in val:
                    i_mers += 1
                    if j_seq in val:
                        common += 1
            simil = common / i_mers
            if simil > 0.20:
                print(fasta_ids[i_seq], fasta_ids[j_seq], i_mers, common, simil)
            matrix[-1].append('{:.6f}'.format(simil))
        else:
            matrix[-1].append('{:.6f}'.format(1.0))

# print as tsv
with open(args.ofile, "w") as matfile:
    sids = [fasta_ids[x] for x in seqlist]
    print("\t".join(sids), file=matfile)
    for r, row in enumerate(matrix):
        print("{}".format(sids[r]), end = "\t", file=matfile)
        print('\t'.join(row), file=matfile)
