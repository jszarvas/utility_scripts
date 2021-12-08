#!/usr/bin/env python

import sys
import os

if len(sys.argv) < 2:
    print("Usage: program <A> <B> <a:A-B|b:intersect|c:difference|d:union>")
    sys.exit()

Alst = []
Blst = []

with open(sys.argv[1], "r") as fp:
    for line in fp:
        Alst.append(line.strip())

with open(sys.argv[2], "r") as fp:
    for line in fp:
        Blst.append(line.strip())

if sys.argv[3] == "a":
    for item in Alst:
        if item not in Blst:
            print(item)
elif sys.argv[3] == "b":
    for item in Alst:
        if item in Blst:
            print(item)
elif sys.argv[3] == "c":
    diff = []
    for item in Alst:
        if item not in Blst:
            diff.append(item)
    for item in Blst:
        if item not in Alst:
            diff.append(item)
    print("\n".join(diff))
elif sys.argv[3] == "d":
    for item in Alst:
        print(item)
    for item in Blst:
        if item not in Alst:
            print(item)
else:
    print("Not valid choice")
