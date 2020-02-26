#!/usr/bin/env python3

import sys
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from cyvcf2 import VCF

parser = argparse.ArgumentParser(
    description='Calculate LOH for sets of VCF files')
parser.add_argument(
    '-l',
    dest="vcf_list",
    default=None,
    help='Input list of VCFs')
parser.add_argument(
    '-o',
    dest="ofile",
    default=None,
    help='Output LOH stat table, windows, chrom, saples')
parser.add_argument(
    '-r',
    dest="ref_bed",
    help='Bed file with features for chromosome info')
parser.add_argument(
    '-w',
    dest="window",
    default=10000,
    type=int,
    help='Window size, default 10kbp')
args = parser.parse_args()

# chrom1 : []
window = {}
# sample1 : [], sample2 : []
stats_per_window = {}
sample_id = None
samples = []
vcf_i = None
hetcount = 0

if args.vcf_list is None or not os.path.exists(args.vcf_list):
    sys.exit("Input file not found.")

with open(args.vcf_list, "r") as fp:
    for line in fp:
        variant = None
        vcf_i = VCF(line.strip())
        if len(vcf_i.samples) == 1:
            sample_id = vcf_i.samples[0]
            samples.append(sample_id)
            print(sample_id)
        stats_per_window[sample_id] = []
        for contig, clen in zip(vcf_i.seqnames, vcf_i.seqlens):
            # window index
            i = 0
            # make windows in the first round
            if not window or window.get(contig) is None:
                window[contig] = []
                for j in range(args.window,clen, args.window):
                    window[contig].append(j)
                window[contig].append(clen)

            try:
                if variant is None:
                    variant = next(vcf_i)
                while (variant.CHROM == contig):
                    # PASS
                    if variant.FILTER is None:
                        # count the heterozygous 0/1, 1/2
                        if variant.genotypes[0][0] == 0 or variant.genotypes[0][-2] == 2:
                            # is the variant outside our current window?
                            # change block
                            if variant.start >= window[contig][i]:
                                stats_per_window[sample_id].append(hetcount)
                                i += 1
                                hetcount = 0
                            hetcount += 1
                    # read next variant
                    variant = next(vcf_i)
                # change contig
                stats_per_window[sample_id].append(hetcount)
                hetcount = 0
            except StopIteration:
                stats_per_window[sample_id].append(hetcount)
                hetcount = 0
                # fill rest of blocks with zero
                for j in range(i, len(window[contig])):
                    stats_per_window[sample_id].append(hetcount)


with open(args.ofile, "w") as op:
    print("contig", "block", "\t".join(samples), sep="\t", file=op)
    i = 0
    for contig in vcf_i.seqnames:
        for block in window[contig]:
            print(contig, str(block), sep="\t", end="\t", file=op)
            print("\t".join([str(stats_per_window[sid][i]) for sid in samples]), file=op)
            i += 1
