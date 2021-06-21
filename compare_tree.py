#!/usr/bin/env python3.6

import sys
import argparse
from ete3 import Tree

def sort_taxa_names(tree):
    for node in tree.traverse():
        if node.name:
            split_name = node.name.split(" ")
            if len(split_name) > 1:
                node.name = " ".join(sorted(split_name))
    return tree

# Parse command line options
parser = argparse.ArgumentParser(
    description='Compare trees RF distance')
parser.add_argument(
   '-t1',
   dest="tree_one",
   required=True,
   help='Tree 1')
parser.add_argument(
    '-t2',
    dest="tree_two",
    required=True,
    help='Tree 2')
args = parser.parse_args()

# load the trees to compare
tree1 = Tree(args.tree_one)
tree2 = Tree(args.tree_two)

# sort words in taxa names
tree1 = sort_taxa_names(tree1)
tree2 = sort_taxa_names(tree2)

# compare trees
rf, max_rf, common_leaves, parts_t1, parts_t2 = tree1.robinson_foulds(tree2)
print(f"{rf}\t{max_rf}\t{rf/max_rf:.3f}")

sys.exit(0)
