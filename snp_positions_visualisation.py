#!/usr/bin/env python3.6

import sys, time
import os
import argparse
from matplotlib import pyplot as plt

def timing(message):
    t1 = time.time()
    print("{0} Time used: {1} seconds".format(message, int(t1-t0)), file=sys.stdout)
    return

def exiting(message):
    print(message, file=sys.stderr)
    print("FAIL", file=sys.stderr)
    sys.exit(1)


# Start time to keep track of progress
t0 = time.time()

# Parse command line options
parser = argparse.ArgumentParser(
    description='Visualises list of VCF positions')
parser.add_argument("-vl", dest="vcf_list", help="list of VCFs")
parser.add_argument("-o", dest="ofile", help="write to pdf")
parser.add_argument("-r", dest="relative", action="store_true", help="Relative numbers")
parser.add_argument("-t", dest="txt", help="write to txt")
args = parser.parse_args()


# File with vcfs
if args.vcf_list is None:
    exiting('No input file(list) was provided')
elif args.ofile is None and args.txt is None:
    exiting('Output target is needed')

# SNP counts
positions = []
counts = []
file_count = 0
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
#Chromosome	29884	.	A	n	.	PASS	.
with open(args.vcf_list) as fp:
    for line in fp:
        filename = line.strip()
        if not filename:
            next
        try:
            vcf_p = open(filename, "r")
        except IOError as e:
            print("{} not found".format(filename), file=sys.stderr)
        else:
            file_count += 1
            for vcf_line in vcf_p:
                if not vcf_line.startswith("#"):
                    values = vcf_line.split("\t")
                    pos = int(values[1])
                    try:
                        array_i = positions.index(pos)
                    except ValueError:
                        positions.append(pos)
                        counts.append(1)
                    else:
                        counts[array_i] += 1
            vcf_p.close()
        #print(len(positions), len(counts))

if args.relative and file_count:
    for i in range(len(counts)):
        counts[i] = counts[i] / file_count

if args.ofile is not None:
    plt.title("VCF position frequencies")
    plt.plot(positions, counts, 'b+', markersize=2)
    plt.xlabel("Position")
    plt.ylabel("Count")
    plt.savefig(args.ofile, format="pdf")

if args.txt is not None:
    s_pos = sorted(positions)
    with open(args.txt, "w") as op:
        for i in s_pos:
            print(positions[positions.index(i)], counts[positions.index(i)], sep="\t", file=op)

timing("# Figure creation is finished.")
print("DONE", file=sys.stderr)
sys.exit(0)
