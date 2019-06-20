#!/usr/bin/env python3

import sys
import os
import shlex
import subprocess
from joblib import Parallel,delayed

MAX_JOBS = 20

def call_cmd(cmd):
    p = subprocess.run(shlex.split(cmd))
    return p.returncode

if len(sys.argv) < 2:
    print("Usage: program <list of commands to run in parallel> <no. parallel proc>")
    sys.exit()

if os.path.exists(sys.argv[1]):
    try:
        sod = open(sys.argv[1], "r")
    except IOError:
        sys.exit("Input script not found")

cmd_list = []
for line in sod:
    cmd_list.append(line.strip())
sod.close()

n_jobs = 20
if len(sys.argv) > 2:
    try:
        n_jobs=int(sys.argv[2])
    except TypeError:
        pass

if cmd_list:
    jobs = Parallel(n_jobs=min(MAX_JOBS, n_jobs))(delayed(call_cmd)(c) for c in cmd_list)
print(len(cmd_list), sum(jobs))
