#!/usr/bin/env python3

import re
import sys

if len(sys.argv) < 3:
    sys.exit("Usage: program <bootstrapped newick> <nhx>")

filecont = []
with open(sys.argv[1], "r") as fp:
    for line in fp:
        filecont.append(line.strip())
bs_nw = "".join(filecont)

nhx_label = "[&&NHX:B={}]"
bs_nhx = []
b = bs_nw
match = re.search(r"\)(\d+\.?\d*\/\d+\.?\d*)(:\d+\.?\d*)", bs_nw)
while match:
    bs_nhx.append(b[:match.start()+1])
    split = re.search(r"(:\d+\.?\d*)[,\(\)]", bs_nhx[-1])
    while split:
        tmp = [bs_nhx[-1][:split.end()-1], nhx_label.format(""), bs_nhx[-1][split.end()-1:]]
        bs_nhx[-1] = "".join(tmp)
        split = re.search(r"(:\d+\.?\d*)[,\(\)]", bs_nhx[-1])
    bs_nhx.append(match.group(2))
    bs_nhx.append(nhx_label.format(match.group(1)))
    b = b[match.end():]
    match = re.search(r"\)(\d+\.?\d*\/\d+\.?\d*)(:\d+\.?\d*)", b)

bs_nhx.append(b)
split = re.search(r"(:\d+\.?\d*)[,\(\)]", bs_nhx[-1])
while split:
    tmp = [bs_nhx[-1][:split.end()-1], nhx_label.format(""), bs_nhx[-1][split.end()-1:]]
    bs_nhx[-1] = "".join(tmp)
    split = re.search(r"(:\d+\.?\d*)[,\(\)]", bs_nhx[-1])

bs_nhx[-1] = bs_nhx[-1][:-1]
bs_nhx.append(nhx_label.format(""))
bs_nhx.append(";")

with open(sys.argv[2], "w") as op:
    print("".join(bs_nhx), file=op)
