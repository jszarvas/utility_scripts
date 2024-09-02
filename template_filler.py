#!/usr/bin/env python3

import sys, os
import numpy as np
import pandas as pd

if len(sys.argv) < 3:
    print("Usage: program <base script> <params tsv> <output folder> <script prefix>")
    print("column names need to correspond to the wildcards in the template")
    sys.exit()
else:
    print(" ".join(sys.argv), file=sys.stderr)

script_main_text = None
with open(sys.argv[1], "r") as fp:
    script_main_text = fp.read()

params_df = pd.read_csv(sys.argv[2], sep="\t")
filled_df = params_df.apply(script_main_text.format_map, axis=1)

for i in range(filled_df.shape[0]):
    output_script_path = os.path.join(sys.argv[3], "{}_{}.sh".format(sys.argv[4], i + 1))
    with open(output_script_path, "w") as of:
        print(filled_df.iloc[i,], file=of)

