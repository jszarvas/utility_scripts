#!/usr/bin/env python3

import sys, time
import os
import gc
import numpy as np
import argparse
import shutil

if len(sys.argv) < 3:
    sys.exit("Usage: program <lower triangular strict phylip> <output full relaxed phylip>\nINTS ONLY, NO FLOATS")

inputmat = sys.argv[1]
outputmat = sys.argv[2]
matrix = []
taxa = []
firstline = None
# read it in
try:
    f = open(inputmat, "r")
except IOError as e:
    exiting("Phylip matrix not found.")
else:
    firstline = f.readline()
    for l in f:
        if l and not l.isspace():
            taxaname, distances = l.split(None,1)
            taxa.append(taxaname)
            matrix.append(distances.split("\t"))
f.close()

# transform
a = np.array(matrix, dtype=np.uint16)
b = np.tril(a, k = -1)
del a
gc.collect()
full_mat = (b.T + b).tolist()
del b
gc.collect()
# write out
with open(outputmat, "w") as matfile:
    matfile.write(firstline)
    for i, row in enumerate(full_mat):
        print(f"{taxa[i]}", end="\t", file=matfile)
        for e in row[:-1]:
            print(f"{e:.0f}", end = "\t", file=matfile)
        print(f"{row[-1]:.0f}", file=matfile)

print("Done", file=sys.stderr)
