#!/usr/bin/env python3

import sys
from cyvcf2 import VCF
import numpy as np
import pandas as pd
from glob import glob

def parse_multisamplevcf_into_df(vcf_filename):
    vcf_i = VCF(vcf_filename)
    # vcf_i.samples
    vcf_d = {'CONTIG': [], 'POS': [], 'REF': [], 'ALT': [], 'sample':[], 'gt_type': [], 'gt_ref_d': [], 'gt_alt_d': []}
    variant = next(vcf_i)
    try:
        while True:
            if variant.FILTER is None:
                if not variant.is_deletion and not variant.is_sv:
                    for j, sample in enumerate(vcf_i.samples):
                        if variant.gt_ref_depths[j] != 0 or variant.gt_alt_depths[j] != 0:
                            vcf_d['sample'].append(vcf_i.samples[j])
                            vcf_d['CONTIG'].append(variant.CHROM)
                            vcf_d['POS'].append(variant.start)
                            vcf_d['REF'].append(variant.REF)
                            vcf_d['ALT'].append(variant.ALT[variant.genotypes[j][0] - 1])
                            vcf_d['gt_type'].append(variant.genotypes[j][0])
                            vcf_d['gt_ref_d'].append(variant.gt_ref_depths[j])
                            vcf_d['gt_alt_d'].append(variant.gt_alt_depths[j])
            variant = next(vcf_i)
    except StopIteration:
        pass
    df = pd.DataFrame.from_dict(vcf_d)
    return df

def parse_annotation_into_df(vcf_filename):
    vcf_i = VCF(vcf_filename)
    # vcf_i.samples
    vcf_d = {'CONTIG': [], 'POS': [], 'ALT': [], 'ANN': []}
    variant = next(vcf_i)
    try:
        while True:
            if variant.FILTER is None:
                if not variant.is_deletion and not variant.is_sv:
                    if variant.INFO.get('ANN') is not None:
                        for alt_anns in variant.INFO.get('ANN').split(","):
                            ann_cols = alt_anns.split("|")
                            vcf_d['CONTIG'].append(variant.CHROM)
                            vcf_d['POS'].append(variant.start)
                            vcf_d['ALT'].append(ann_cols[0])
                            vcf_d['ANN'].append("|".join([ann_cols[1],ann_cols[3],ann_cols[10]]))
            variant = next(vcf_i)
    except StopIteration:
        pass
    df = pd.DataFrame.from_dict(vcf_d)
    return df

if (len(sys.argv) < 3):
    print("program <multi-sample vcf> <output tsv> <opt: annonation TRUE>")
    sys.exit(1)

joint_calling_df = parse_multisamplevcf_into_df(sys.argv[1])
joint_calling_df.to_csv(path_or_buf=sys.argv[2], sep='\t', index=False)
if len(sys.argv) == 4 and sys.argv[3].lower() in ["t", "true"]:
    ann_df = parse_annotation_into_df(sys.argv[1])
    ann_df.to_csv(path_or_buf="{}.ann".format(sys.argv[2]), sep='\t', index=False)
