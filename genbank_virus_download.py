#!/usr/bin/env python3.6

import os
import sys
import argparse
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json
import sqlite3

BATCH_SIZE = 100

def exiting(message):
    print(message, file=sys.stderr)
    print("FAIL", file=sys.stderr)
    sys.exit(1)

def fit_folder(folder, bdir):
    """Avoid having more than 10000 files in one folder."""
    full_folder = os.path.join(bdir, folder)
    i = 1
    while len(os.listdir(full_folder)) > 10000:
        i += 1
        full_folder = os.path.join(bdir, "{0}-{1}".format(folder, i))
        if not os.path.exists(full_folder):
            os.mkdir(full_folder)
    return full_folder

parser = argparse.ArgumentParser(
    description='Makes Entrez query to Genbank, downloads queried gb and extracted nucleotide fasta')
parser.add_argument(
    '-q',
    dest="cf",
    help='Path to the query file')
parser.add_argument(
    '-o',
    dest="out",
    help='Data output directory absolute path')
parser.add_argument(
    '-d',
    dest="database",
    help="Sample database"
)
args = parser.parse_args()

odir = os.path.realpath(args.out)
base = os.path.dirname(odir)
folder = os.path.basename(odir)
odir = fit_folder(folder, base)


#start up the db during first run_id
conn = sqlite3.connect(args.database)
conn.execute("PRAGMA foreign_keys = 1")
cur = conn.cursor()

cur.execute('''CREATE TABLE IF NOT EXISTS samples
    (accession TEXT PRIMARY KEY,
    path TEXT,
    feature_path TEXT,
    dl_date TEXT DEFAULT CURRENT_DATE,
    included INTEGER DEFAULT NULL);''')
conn.commit()

queries = []
if args.cf is not None:
    if not os.path.exists(args.cf):
        exiting("Queries file doesnt exist.")
    else:
        with open(args.cf, "r") as fp:
            for line in fp:
                #ORF2   txid11983[Organism] AND biomol_genomic[PROP] AND ("450"[SLEN] : "8000"[SLEN])
                tmp = line.strip().split("\t")
                if not tmp[0].isspace():
                    tmp[0] = tmp[0].split("|")
                    queries.append(tmp)
else:
    exiting("No query supplied")

print(queries)


# use Entrez for searching genbank
if os.environ.get('NCBI_API_KEY') is not None:
    Entrez.api_key = os.environ.get('NCBI_API_KEY')
else:
    Entrez.api_key = "2a617493db26b60c5603b19744b4e60b3d08"

if os.environ.get('NCBI_EMAIL') is not None:
    Entrez.email = os.environ.get('NCBI_EMAIL')
else:
    Entrez.email = "jusz@dtu.dk"

for query in queries:
    query_str = query[1]
    # 'txid11983[Organism] AND biomol_genomic[PROP] AND ("450"[SLEN] : "8000"[SLEN])'
    handle = Entrez.esearch(db="nucleotide", term=query_str, idtype="acc", retmax=50000)
    record = Entrez.read(handle)
    handle.close()

    print(query_str, int(record["Count"]))

    acc = []
    for r in record['IdList']:
        cur.execute('''SELECT * FROM samples WHERE accession=?''', (r,))
        if cur.fetchone() is None:
            acc.append(r)

    n_acc = len(acc)
    batches = [acc[x:x + BATCH_SIZE] for x in range(0, n_acc, BATCH_SIZE)]

    sample_insert = []
    for batch in batches:
        acc_str = ",".join(batch)
        fetch_handle = Entrez.efetch(db="nuccore", id=acc_str, rettype="gb", retmode="text")
        records = SeqIO.parse(fetch_handle, "gb")
        for r in records:
            extracted_fasta_path = None
            for f in r.features:
                # save metadata in json
                if f.type == "source":
                    metadata = {'id': r.id,
                          'description': r.description }
                    metadata.update(f.qualifiers)
                    with open(os.path.join(odir, "{}.json".format(r.id)), "w") as op:
                        json.dump(metadata, op)
                # correct ORF or gene name
                if query[0][0] != "-":
                    if f.type == "gene" or f.type == "CDS":
                        if (f.qualifiers.get('gene') is not None and f.qualifiers.get('gene')[0] in query[0]) or (f.qualifiers.get('product') is not None and f.qualifiers.get('product')[0] in query[0]) or (f.qualifiers.get('note') is not None and f.qualifiers.get('note')[0] in query[0]):
                            feature_record = SeqRecord(f.extract(r.seq), id=r.id, description = " ".join(query[0]))
                            extracted_fasta_path = os.path.join(odir, "{}_{}.fsa".format(query[0][0], r.id))
                            SeqIO.write(feature_record, extracted_fasta_path, "fasta")
                else:
                    break
            fasta_path = os.path.join(odir, "{}.fsa".format(r.id))
            SeqIO.write(r, fasta_path, "fasta")
            sample_insert.append((r.id, fasta_path, extracted_fasta_path))
        fetch_handle.close()

        # insert to db, then re-init
        cur.executemany('''INSERT OR REPLACE INTO samples (accession, path, feature_path) VALUES (?,?,?)''', sample_insert)
        conn.commit()
        sample_insert = []
conn.close()
