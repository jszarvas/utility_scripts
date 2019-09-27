#!/usr/bin/env python3

import os
import sys
import argparse
from cyvcf2 import VCF
import numpy as np
from matplotlib import pyplot as plt

PLH = 500.0
N_BINS = 50

class VCF_stats:
    """docstring for VCF_stats."""

    def __init__(self, params, odir, filepath):
        self.path = filepath
        self.samplename = os.path.basename(filepath).split(".")[0]
        self.prefix = os.path.join(odir, self.samplename)
        self.filter = False
        self.stats = {}
        self.contigs = []
        self.tests = params.split(",")
        for p in self.tests:
            self.stats[p] = []

    def parse(self, filter):
        if filter:
            self.filter = True
        vcf_iter = VCF(self.path)
        for v in vcf_iter:
            self.contigs.append(v.CHROM)
            if not filter or v.FILTER is None:
                for p in self.tests:
                    if v.INFO.get(p) is not None:
                        self.stats[p].append(float(v.INFO.get(p)))
                    else:
                        self.stats[p].append(PLH)
        return

    def visualise(self):
        if len(self.tests) > 1:
            self.bookletize()
        else:
            self.create_histo(self.tests[0])
        return

    def create_histo(self, param):
        pdf_path = "{}_{}.pdf".format(self.prefix, param)
        if self.filter:
            pdf_path = "{}_{}_PASS.pdf".format(self.prefix, param)
        # make it into numpy
        npmat = np.array(self.stats[param])
        plt.rcParams["figure.figsize"] = [15, 4 * len(self.tests)]
        plt.hist(np.extract(npmat != PLH, npmat), bins=N_BINS, density=True)
        plt.title("{} {}".format(self.samplename, param))
        plt.savefig(pdf_path, format="pdf")
        return

    def bookletize(self):
        pdf_path = "{}_{}.pdf".format(self.prefix, "BL")
        if self.filter:
            pdf_path = "{}_{}_PASS.pdf".format(self.prefix, "BL")
        fig, axes = plt.subplots(len(self.tests), figsize=(15, 4 * len(self.tests)))
        fig.suptitle(self.samplename)
        for i, p in enumerate(self.tests):
            npmat = np.array(self.stats[p])
            axes[i].set_title(p)
            axes[i].hist(np.extract(npmat != PLH, npmat), bins=N_BINS, density=True)
        fig.savefig(pdf_path, format="pdf")
        return


parser = argparse.ArgumentParser(
    description='Visual output from alignment and variant calling staistics in vcf files')
parser.add_argument(
    '-i',
    dest="vcf_file",
    default=None,
    help='Input VCF.gz')
parser.add_argument(
    '-l',
    dest="vcf_list",
    default=None,
    help='Input list of VCFs')
parser.add_argument(
    '-o',
    dest="odir",
    default=None,
    help='Output pdfs to dir')
parser.add_argument(
    '-p',
    dest="params",
    default="DP,QD,FS,SOR,MQ,MQRankSum",
    help='List of INFO fields to capture')
parser.add_argument(
    '--pass',
    dest="filter",
    action="store_true",
    help='FILTER for PASS')
args = parser.parse_args()

if args.vcf_file is not None and os.path.exists(args.vcf_file):
    vcf = VCF_stats(args.params, args.odir, args.vcf_file)
    vcf.parse(args.filter)
    if len(vcf.tests) == 1:
        vcf.visualise()
    else:
        vcf.bookletize()
    sys.exit(0)

elif args.vcf_list is not None and os.path.exists(args.vcf_list):
    sys.exit("No input")
    with open(args.vcf_list, "r") as fp:
        for line in fp:
            vcf = VCF_stats(args.params, args.odir, line.strip())
            vcf.parse(args.filter)
            if len(vcf.tests) == 1:
                vcf.visualise()
            else:
                vcf.bookletize()

else:
    sys.exit("No input")
