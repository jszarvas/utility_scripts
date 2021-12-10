#!/usr/bin/env python3

import os
import sys
import shutil
import argparse
import subprocess
import shlex
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import seaborn as sns

MAXDEPTH = 100

parser = argparse.ArgumentParser(
    description='Feature coverage using bedtools')
parser.add_argument(
    '-i',
    dest="hist_file",
    default=None,
    help='Input hist.txt')
parser.add_argument(
    '-b',
    dest="bam_file",
    default=None,
    help='Input bam')
parser.add_argument(
    '-o',
    dest="odir",
    help='Output pdf to dir')
parser.add_argument(
    '-e',
    dest="bedtools_path",
    help='Bedtools executable path')
parser.add_argument(
    '-l',
    dest="minlength",
    default=2000,
    type=int,
    help='Min contig size to consider, bp')
args = parser.parse_args()

if args.bam_file is not None and args.ref_bed is not None:
    if args.bedtools_path is None:
        if shutil.which("bedtools") is not None:
            args.bedtools_path = shutil.which("bedtools")
        else:
            sys.exit("Path to bedtools binary is needed")
    else:
        sys.exit("Path to bedtools binary is needed")

    #bedtools genomecov -ibam ${F} >  ${S}.bam.hist.txt

    args.hist_file = "{}.hist.txt".format(args.bam_file)
    cmd = "{0} genomecov -ibam {1}".format(args.bedtools_path, args.bam_file)
    with open(args.hist_file, "w") as ofp:
        ex = subprocess.run(shlex.split(cmd), stdout=ofp)
        if ex.returncode != 0:
            args.hist_file = None

if args.hist_file is not None:
    samplename = os.path.basename(args.hist_file).split(".")[0]
    pdf_path = "{}_coverage.pdf".format(os.path.join(args.odir, samplename))

    hist_data = pd.read_table(args.hist_file, header=None, names=["feature", "depth", "bases", "feature_size", "fraction"])
    hist_data_feat = hist_data[hist_data.feature != "all"]
    # 1 - cumsum to show the fraction that is covered with at least that depth
    hist_data_feat['frac_cumsum'] = hist_data_feat.groupby(['feature'])['fraction'].apply(lambda x: 1 - x.cumsum())
    # remove cumsum values on the 1 bases long high depths
    hist_data_feat.loc[hist_data_feat.frac_cumsum < 0,'frac_cumsum'] = 0

    sns.set(style="whitegrid", rc={"xtick.bottom" : True, "ytick.left" : True})
    plt.suptitle(samplename)
    ax = sns.lineplot(x="depth", y="frac_cumsum", hue="feature", data=hist_data_feat[hist_data_feat.feature_size > args.minlength])
    plt.xscale('log')
    formatter = plt.ScalarFormatter()
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(plt.ScalarFormatter())
    plt.xlabel("Depth of coverage")
    plt.ylabel("Fraction covered by at least this depth")
    plt.legend([],[], frameon=False)
    plt.savefig(pdf_path, format="pdf")

    print("# PDF saved to {}".format(pdf_path))
