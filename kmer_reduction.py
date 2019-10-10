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
    dest="odir",
    default=None,
    help='Output directory')
parser.add_argument(
    '-k',
    dest="kmer",
    default=16,
    type=int,
    help='Kmer size')
parser.add_argument(
    '-t',
    dest="thr",
    default=1.00,
    type=float,
    help='Homology reduction threshold, discard above')
parser.add_argument(
    '-s',
    dest="split",
    default=None,
    help='Split entry id on this character')
args = parser.parse_args()

del_ids = set()
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

# list of hashes
seqlist = fasta_ids.keys()
# template
for i_seq in seqlist:
    if i_seq not in del_ids:
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
                if simil > args.thr:
                    print(fasta_ids[i_seq], fasta_ids[j_seq], i_mers, common, simil)
                    del_ids.add(j_seq)


# print as separate entries
splitchar = ""
if args.split is not None:
    splitchar = args.split

# get the ids to keep
output_ids = set(seqlist) - del_ids

if args.mfsa is not None and os.path.exists(args.mfsa):
    seq_iter = SeqIO.parse(args.mfsa, "fasta")
    for rec in seq_iter:
        if hash_seqid(rec.id) in output_ids:
            tmp = rec.id.split(splitchar)
            rec.id = tmp[0]
            if len(tmp) > 1:
                rec.description = splitchar.join(tmp[1:])
            ofn = os.path.join(args.odir, "{0}.fa".format(rec.id))
            with open(ofn, "w") as of:
                SeqIO.write(rec, of, "fasta")

elif args.fsa_list is not None and os.path.exists(args.fsa_list):
    with open(args.fsa_list, "r") as fp:
        for line in fp:
            parser = SeqIO.parse(line.strip(), "fasta")
            for rec in parser:
                if hash_seqid(rec.id) in output_ids:
                    tmp = rec.id.split(splitchar)
                    rec.id = tmp[0]
                    if len(tmp) > 1:
                        rec.description = splitchar.join(tmp[1:])
                    ofn = os.path.join(args.odir, "{0}.fa".format(rec.id))
                    with open(ofn, "w") as of:
                        SeqIO.write(rec, of, "fasta")
