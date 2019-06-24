#!/usr/bin/env python

# Copyright (c) 2014, Ole Lund, Technical University of Denmark
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#
# Import libraries
#
import sys, time
import os
import gc
import numpy as np
import array
from optparse import OptionParser
from operator import itemgetter
import re
import subprocess

#
# Parse command line options
#
parser = OptionParser()
parser.add_option("-i", "--inputfile", dest="inputfilename", help="read from INFILE", metavar="INFILE")
parser.add_option("-o", "--outputfile", dest="outputfilename", help="write to OUTFILE", metavar="OUTFILE")
parser.add_option("-s", "--separate", dest="separate", action="store_true", help="Sim reads for fasta entries independently")
parser.add_option("-z", "--zip", dest="zip", action="store_true", help="Gzip")
(options, args) = parser.parse_args()
#
# Open files
#
#
# File with reads
#
if options.inputfilename != None:
  inputfile = open(options.inputfilename,"r")
else:
  inputfile = sys.stdin
#
# File for general output
#
if options.outputfilename != None and not options.separate:
  outputfile = open(options.outputfilename,"w")
elif options.separate:
  outputfile = None
else:
  outputfile = sys.stdout

#
# Read Input fasta file
#
inputseq = []
inputseqsegments = []
inputname = []
inputdesc = []
Ninputs=0
i=0
if options.inputfilename != None:
  t1 = time.time()
  for line in inputfile:
    fields=line.split()
    if len(line)>1:
      if fields[0][0] == ">":
        if (i>0):
          inputseq[-1] = ''.join(inputseqsegments)
          del inputseqsegments
          inputseqsegments = []
          i=0
        inputseq.append("")
        inputname.append(fields[0][1:])
        inputdesc.append(re.sub(r"^[^\s]+\s","",line.strip()))
      else:
        inputseqsegments.append(fields[0])
        i+=1
  inputseq[-1] = ''.join(inputseqsegments)

del inputseqsegments

#
# parameters for the reads
#
nseq = len(inputseq)
readlength = 150
# lower coverage for the perfect reads
stepsize = 5
#stepsize = 30

coverage=30
#coverage=5

filename = options.outputfilename
for i in xrange(0, nseq):
  if options.separate:
    filename = os.path.join(options.outputfilename, "{}.fastq".format(re.sub(r"\W+","_",inputname[i])))
    outputfile = open(filename, "w")
  start = 0
  end = start + readlength
  while (end < len(inputseq[i])):
    outputfile.write("%s\n" % ("@"))
    outputfile.write("%s\n" % (inputseq[i][start:end]))
    outputfile.write("%s\n" % ("+"))
    for j in xrange(0,readlength):
      outputfile.write("%s" % ("6"))
    outputfile.write("\n")
    start += stepsize
    end += stepsize
  if (len(inputseq[i])>readlength):
    for k in(0,coverage):
      #
      # make sure to get coverage also in the ends (the overlapping windows
      # tapper off in the ends)
      #
      end = len(inputseq[i])
      outputfile.write("%s\n" % ("@"))
      outputfile.write("%s\n" % (inputseq[i][0:readlength]))
      outputfile.write("%s\n" % ("+"))
      for j in xrange(0,readlength):
        outputfile.write("%s" % ("6"))
      outputfile.write("\n")
      outputfile.write("%s\n" % ("@"))
      outputfile.write("%s\n" % (inputseq[i][end-readlength:end]))
      outputfile.write("%s\n" % ("+"))
      for j in xrange(0,readlength):
        outputfile.write("%s" % ("6"))
      outputfile.write("\n")
  if options.separate:
    outputfile.close()
    if options.zip:
      p = subprocess.call(['gzip', filename])
