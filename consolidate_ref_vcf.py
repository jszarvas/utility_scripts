#!/usr/bin/env python3

import sys
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from cyvcf2 import VCF

parser = argparse.ArgumentParser(
    description='Consensus sequence from VCF file and reference fasta')
parser.add_argument(
    '-r',
    dest="ref_file",
    default=None,
    help='Input fasta reference')
parser.add_argument(
    '-i',
    dest="vcf_file",
    default=None,
    help='Input VCF.gz')
parser.add_argument(
    '-o',
    dest="ofile",
    default=None,
    help='Output consensus (multi) fasta')
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

vcf_i = None
seq_iter = None
if args.vcf_file is not None and os.path.exists(args.vcf_file):
    vcf_i = VCF(args.vcf_file)
else:
    sys.exit("VCF needed")

if args.ref_file is not None and os.path.exists(args.ref_file):
    records = list(SeqIO.parse(args.ref_file, "fasta"))
else:
    sys.exit("Ref needed")

# references = vcf_i.seqnames
concat = []
contig = []
variant = None
prev_end = 0
#contig_len = 0
lens = []
for rec in records:
    if rec.id in vcf_i.seqnames:
        print(rec.id)
        #print(vcf_i.seqnames)
        print(len(rec.seq))
        try:
            if variant is None:
                variant = next(vcf_i)
            # right contig
            print(variant.CHROM, rec.id)
            while (variant.CHROM == rec.id):
                #print(prev_end, variant.start, variant.end, variant.REF, variant.ALT)
                # overarching deletions
                if variant.start <= prev_end or '*' in variant.ALT:
                    variant = next(vcf_i)
                    continue
                contig.append(str(rec.seq[prev_end:variant.start]))
                #contig_len += variant.start - prev_end + len(variant.REF)
                # .end is 1 based, .start is 0
                prev_end = variant.end
                # passed VariantFiltration, and not deletion or struct var
                if variant.FILTER is None:
                    if not variant.is_deletion and not variant.is_sv:
                        # homozygous
                        # insertion A -> ACT gets trimmed
                        if sum(variant.genotypes[0][:2]) == 2:
                            contig.append(variant.ALT[0][:len(variant.REF)])
                        elif variant.is_indel:
                            contig.append(variant.ALT[0][0])
                            contig.append("n" * (len(variant.REF) - 1))
                            #print(contig[-2], contig[-1])
                            #print(len("".join(contig)), contig_len,  prev_end, variant.start, variant.end, variant.REF, variant.ALT, variant.genotypes, variant.is_deletion, variant.is_sv, variant.is_indel,  variant.FILTER)
                        else:
                            contig.append("n" * len(variant.REF))
                    elif variant.is_deletion:
                        # TTG -> T
                        # first (or preceding) base
                        contig.append(variant.ALT[0][0])
                        #print(contig[-1])
                        # deletion gets filled with Ns
                        contig.append("n" * (len(variant.REF) - 1))
                else:
                    # low conf. variants are n-s
                    contig.append("n" * len(variant.REF))

                # read next variant
                variant = next(vcf_i)
            contig.append(str(rec.seq[prev_end:]))
        except StopIteration:
            #print("Run out")
            contig.append(str(rec.seq[prev_end:]))
        rec.seq = Seq("".join(contig))
        if args.concat:
            concat.append("".join(contig))
        print("Cons len:", len(rec.seq))
        lens.append(len(rec.seq))
        contig = []
        #contig_len = 0
        prev_end = 0
        debug = True


with open(args.ofile, "w") as op:
    SeqIO.write(records, op, "fasta")

if args.concat:
    full_seq = "".join(concat)
    with open("{}.flat".format(args.ofile.rsplit(".",1)[0]), "w") as fop:
        print(">{}".format(os.path.basename(args.ofile.split(".")[0])), file=fop)
        for i in range(0, len(full_seq), 60):
            print(full_seq[i:i+60], file=fop)

if args.part:
    cummulated_len = 1
    with open("{}.partitions".format(args.ofile.rsplit(".",1)[0]), "w") as pop:
        for i, part_len in enumerate(lens):
            print("DNA, part{0} = {1}-{2}".format(i+1, cummulated_len, part_len), file=pop)
            cummulated_len += part_len
