#!/usr/bin/env python3

import sys
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

if len(sys.argv) < 2:
    sys.exit("program <mash_dist_mat> <bacth size>")

df = None
if os.path.exists(sys.argv[1]):
    df = pd.read_csv(sys.argv[1], sep="\t", header=0)
    print("# Matrix loaded")
else:
   sys.exit("No inputfile")

#df = df.set_index('#query')
df_i = df.set_index('#query')

print("# Creating the png")

batch_size = 4000
if len(sys.argv) > 2:
    batch_size = int(sys.argv[2])

for i in range(0,df_i.shape[0],batch_size):
    f, ax = plt.subplots()
    sns.heatmap(df_i.iloc[i:i+(batch_size+int(batch_size/2)+1),i:i+batch_size+int(batch_size/2)+1], xticklabels=False, yticklabels=False, center=0.5)
    plt.savefig("{}_{}.png".format(sys.argv[1],i), transparent=True)

print("Done")
