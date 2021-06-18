#!/usr/bin/env python3.6

import sys
import os
import time
import argparse
from ete3 import Tree

# Start time to keep track of progress
t0 = time.time()

# Parse command line options
parser = argparse.ArgumentParser(
    description='Trim tree')
parser.add_argument(
   '-o',
   dest="ofix",
   required=True,
   help='Output prefix')
parser.add_argument(
    '-s',
    dest="filelist",
    required=True,
    help='List of samples')
parser.add_argument(
    '-t',
    dest="treefile",
    required=True,
    help='Input phylogenic tree')
args = parser.parse_args()

# read sample accessions from first col of a tsv/lst
user_samples = []
with open(args.filelist, "r") as fp:
    for line in fp:
        user_samples.append(line.strip().split("\t")[0])

# load the tree to prune
tree = Tree(args.treefile)

# find the nodes to keep, both reduced and asteriks version
kept_nodes = []
for node in tree.traverse():
    if node.name:
        split_name = node.name .split(" ")
        if len(split_name) > 1:
            kept_name = []
            for acc in split_name:
                if acc in user_samples:
                    kept_name.append(acc)
            if kept_name:
                node.name = " ".join(kept_name)
                kept_nodes.append(node.name)
        else:
            node.name = node.name.split("*")[-1]
            if node.name in user_samples:
                kept_nodes.append(node.name)

# prune tree
tree.prune(kept_nodes)

# save tree
tree.write(format=0, outfile=f"{args.ofix}.nwk")

print(f"# Tree pruned. Time used: {int(time.time()-t0)} seconds", file=sys.stdout)

sys.exit(0)
