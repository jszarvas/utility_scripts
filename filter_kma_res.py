#!/usr/bin/env python3

import sys
import os
import argparse
import numpy as np
import pandas as pd
from Bio import Entrez
from time import sleep
from pandas.core.common import SettingWithCopyWarning
import warnings

warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

parser = argparse.ArgumentParser(
    description='Filtering KMA (metagenomic) read alignment results for chromosomes')
parser.add_argument(
    '-r','--res',
    dest="input_res",
    required=True,
    help='Path to .res output from KMA')
parser.add_argument(
    '-o',
    dest="output_folder",
    required=True,
    help='Path to output folder for filtered results')
parser.add_argument(
    '-s',
    dest="opt_suffix",
    help="Optional suffix"
)
parser.add_argument(
    '--min-cov',
    dest="min_cov",
    type=float,
    default=10.0,
    help="Minimum genome coverage"
)
parser.add_argument(
    '--min-depth',
    dest="min_depth",
    type=float,
    default=0.20,
    help="Minimum mean depth"
)
parser.add_argument(
    '--model',
    dest="use_model",
    action="store_true",
    help="Use non-linear model (depth ~ coverage) for thresholding"
)
parser.add_argument(
    '--min-frac',
    dest="min_frac",
    type=float,
    default=0.40,
    help="Minimum genome fraction sequenced, decimals"
)
parser.add_argument(
    '--coeff_a',
    dest="coeff_a",
    type=float,
    default=0.9987808,
    help="Coefficient a, in the a * exp(b * x) equation"
)
parser.add_argument(
    '--coeff_b',
    dest="coeff_b",
    type=float,
    default=-0.8670195,
    help="Coefficient b, in the a * exp(b * x) equation"
)
args = parser.parse_args()

# read res file
df_res_raw = pd.read_csv(args.input_res, delimiter='\t')

# filter on min coverage and min depth
df_min_filter = df_res_raw[(df_res_raw["Template_Coverage"] > args.min_cov) & (df_res_raw["Depth"] > args.min_depth)]

#df_min_filter = df_res_raw[(df_res_raw["Template_Coverage"] > 10.0) & (df_res_raw["Depth"] > 0.20)]

df_pass = None
if args.use_model:
    #df_min_filter.loc[:,"Exp_Cov"] = (-0.9987808 * np.exp(-0.8670195 * df_min_filter["Depth"] * (1 / 0.60)) + 1) * 100 * 0.60
    df_min_filter.loc[:,"Exp_Cov"] = (-1 * args.coeff_a * np.exp(args.coeff_b * df_min_filter["Depth"] * (1 / args.min_frac)) + 1) * 100 * args.min_frac

    df_pass = df_min_filter[df_min_filter["Template_Coverage"] > df_min_filter["Exp_Cov"]]

else:
    df_pass = df_min_filter

# acquire full taxonomy lineage from Entrez
taxonomy = []
Entrez.email = os.environ.get("NCBI_EMAIL")
Entrez.api_key = os.environ.get("NCBI_API_KEY")
if not df_pass.empty:
    df_pass.loc[:,"Accession"] = df_pass["#Template"].str.split(" ").str.get(0)
    batch = 8
    for i in range(0,df_pass.shape[0], batch):
        query_str = "[Accession] OR ".join(df_pass["#Template"].str.split(" ").str.get(0)[i:i+batch])
        esearch_reply = Entrez.esearch(db="Nuccore", term=query_str)
        esearch_record = Entrez.read(esearch_reply)
        for recid in esearch_record['IdList']:
            efetch_reply = Entrez.efetch(db="Nuccore", id=recid, retmode="xml")
            efetch_record = Entrez.read(efetch_reply)
            #print(efetch_record[0]['GBSeq_organism'])
            esearch_tax_reply = Entrez.esearch(db="Taxonomy", term=efetch_record[0]['GBSeq_organism'])
            esearch_tax_record = Entrez.read(esearch_tax_reply)
            efetch_tax_reply = Entrez.efetch(db="Taxonomy", id=esearch_tax_record['IdList'][0], retmode="xml")
            efetch_tax_record = Entrez.read(efetch_tax_reply)
            #print(efetch_tax_record[0]['Lineage'])
            taxonomy.append({'Accession': efetch_record[0]['GBSeq_accession-version'], 'TaxId': efetch_tax_record[0]['TaxId'],
                efetch_tax_record[0]['Rank']: efetch_tax_record[0]['ScientificName'],
                'Lineage': efetch_tax_record[0]['Lineage'],
                'TaxNote': "",
                'ParentTaxNoRank': ""
                })
            for tax_level in efetch_tax_record[0]['LineageEx']:
                if tax_level['Rank'] in ['superkingdom', 'phylum', 'family', 'genus', 'species']:
                    if taxonomy[-1].get(tax_level['Rank']) is None:
                        taxonomy[-1][tax_level['Rank']] = tax_level['ScientificName']
            parent_tax = efetch_tax_record[0]['Lineage'].split("; ")[-1]
            if taxonomy[-1]['genus'] and parent_tax != taxonomy[-1]['genus']:
                if parent_tax == taxonomy[-1]['species']:
                    taxonomy[-1]['TaxNote'] = taxonomy[-1].get(efetch_tax_record[0]['Rank'])
                elif parent_tax.split(" ")[0] == "unclassified":
                    taxonomy[-1]['TaxNote'] = parent_tax
                else:
                    taxonomy[-1]['ParentTaxNoRank'] = parent_tax
        sleep(1)

    taxonomy_data = {'Accession': [], 'TaxId': [],
    'superkingdom': [], 'phylum': [], 'family': [], 'genus': [], 'species': [],
    'Lineage': [],
    'TaxNote': [],
    'ParentTaxNoRank': []
    }

    # data = {'col_1': [3, 2, 1, 0], 'col_2': ['a', 'b', 'c', 'd']}
    # pd.DataFrame.from_dict(data)
    for rec in taxonomy:
        for k in taxonomy_data.keys():
            taxonomy_data[k].append(rec.get(k))
    df_pass = df_pass.merge(pd.DataFrame.from_dict(taxonomy_data))

if args.opt_suffix is None:
    args.opt_suffix = os.path.basename(args.input_res).replace("res", "filter.tsv")
output_filename = os.path.join(args.output_folder, args.opt_suffix)
with open(output_filename, "w") as op:
    arg_list = ["##", sys.argv[0]]
    arg_list += ["{0}={1}".format(k,v) for k,v in sorted(vars(args).items())]
    print(" ".join(arg_list), file=op)
    if args.use_model:
        print(f"## (-1 * {args.coeff_a} * np.exp({args.coeff_b} * Depth * (1 / {args.min_frac})) + 1) * 100 * {args.min_frac}", file=op)
df_pass.to_csv(path_or_buf=output_filename, sep='\t', mode='a', index=False)
