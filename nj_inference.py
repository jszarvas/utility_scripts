#!/usr/bin/env python3

import sys, time
import os, gc
import numpy as np
import array
import argparse
from operator import itemgetter
import gzip
import math
import subprocess
from multiprocessing import cpu_count
from joblib import Parallel, delayed


# quick hack
no_jobs = 20
MEM_AVAIL = 80 # Gb
DTREE = "neighbor"

############# ITERATORS #############
def SeqsFromFile(filename, exit_on_err=False):
   '''Extract sequences from a file

   Name:
      SeqsFromFile
   Author(s):
      Martin CF Thomsen
   Date:
      18 Jul 2013
   Description:
      Iterator which extract sequence data from the input file
   Args:
      filename: string which contain a path to the input file
   Supported Formats:
      fasta, fastq

   USAGE:
   >>> import os, sys
   >>> # Create fasta test file
   >>> file_content = ('>head1 desc1\nthis_is_seq_1\n>head2 desc2\n'
                       'this_is_seq_2\n>head3 desc3\nthis_is_seq_3\n')
   >>> with open('test.fsa', 'w') as f: f.write(file_content)
   >>> # Parse and print the fasta file
   >>> for seq, name, desc in SeqsFromFile('test.fsa'):
   ...    print ">%s %s\n%s"%(name, desc, seq)
   ...
   >head1 desc1
   this_is_seq_1
   >head2 desc2
   this_is_seq_2
   >head3 desc3
   this_is_seq_3
   '''

   # EXTRACT DATA
   with open_(filename,"r") as f:
      queryseqsegments = []
      seq, name, desc = '', '', ''
      line = ''
      nextline = f.__next__
      addsegment = queryseqsegments.append
      for line in f:
         if len(line.strip()) == 0: continue
         #print("%s\n"%line, file=sys.stderr)
         fields=line.strip().split()
         if line[0] == ">":
            # FASTA HEADER FOUND
            if queryseqsegments != []:
               # YIELD SEQUENCE AND RESET
               seq = ''.join(queryseqsegments)
               yield (seq, name, desc)
               seq, name, desc = '', '', ''
               del queryseqsegments[:]
            name = fields[0][1:]
            desc = ' '.join(fields[1:])

         elif line[0] == "@":
            # FASTQ HEADER FOUND
            name = fields[0][1:]
            desc = ' '.join(fields[1:])
            try:
               # EXTRACT FASTQ SEQUENCE
               line = nextline()
               seq  = line.strip().split()[0]
               # SKIP SECOND HEADER LINE AND QUALITY SCORES
               line = nextline()
               line = nextline() # Qualities
            except:
               break
            else:
               # YIELD SEQUENCE AND RESET
               yield (seq, name, desc)
               seq, name, desc = '', '', ''

         elif len(fields[0])>0:
            # EXTRACT FASTA SEQUENCE
            addsegment(fields[0])

      # CHECK FOR LAST FASTA SEQUENCE
      if queryseqsegments != []:
         # YIELD SEQUENCE
         seq = ''.join(queryseqsegments)
         yield (seq, name, desc)


############# FUNCTIONS #############
def open_(filename, mode=None, compresslevel=9):
   """Switch for both open() and gzip.open().

   Determines if the file is normal or gzipped by looking at the file
   extension.

   The filename argument is required; mode defaults to 'rb' for gzip and 'r'
   for normal and compresslevel defaults to 9 for gzip.

   """
   if filename[-3:] == '.gz':
      if mode is None: mode = 'rb'
      return gzip.open(filename, mode, compresslevel)
   else:
      if mode is None: mode = 'r'
      return open(filename, mode)

def timing(message):
    if not args.quiet:
        t1 = time.time()
        print("{0} Time used: {1} seconds".format(message, int(t1-t0)), file=sys.stdout)
    return

def exiting(message):
    print(message, file=sys.stderr)
    print("FAIL", file=sys.stderr)
    sys.exit(1)

def read_encode_univ(fp, tot_len):
    if os.path.exists(fp):
        entries = list(zip(*[[seq, name, desc] for seq, name, desc in SeqsFromFile(fp)]))
        strain = "".join(entries[0])
        if tot_len is None:
            tot_len = len(strain)

        encodedinput = np.zeros((tot_len), dtype=np.int8)
        for i in range(tot_len):
            try:
                encodedinput[i] = nuc2num[strain[i].upper()]
            except KeyError:
                pass
        return encodedinput
    else:
        return None

def dist_calc_pw(s1, s2, i, S, no_compared_strains):
    if slens < i + S:
        S = slens - i
    dist_t = np.zeros(shape=(no_compared_strains,S), dtype=np.int32)
    for j in range(S):
        l = i+j
        for k in range(no_compared_strains):
            #print(k,l)
            dist_t[k,j] = np.not_equal(s1[l,], s2[k,]).sum(0) - np.not_equal(s1[l,]!= 0, s2[k,]!= 0).sum(0)
    return dist_t

def dist_calc_all(s1, s2, i, S, no_compared_strains):
    if slens < i + S:
        S = slens - i
    dist_t = np.zeros(shape=(no_compared_strains,S), dtype=np.int32)
    for j in range(S):
        l = i+j
        for k in range(no_compared_strains):
            #print(k,l)
            dist_t[k,j] = np.not_equal(s1[l,], s2[k,]).sum(0)
    return dist_t

def change2subdir(mode, method, suffix, bdir):
    if mode in ["all", "pw"]:
        subdir = os.path.join(bdir, "{0}_{1}_{2}{3}".format(mode, method, ctime, suffix))
    else:
        subdir = os.path.join(bdir, "{0}_{1}{2}".format(method, ctime, suffix))

    try:
        os.mkdir(subdir)
        os.chdir(subdir)
    except:
        exiting("Couldn't make {0}".format(subdir))
    return subdir

def decode_nw_file(seqid2name, nw_file):
    newick = []
    with open(nw_file, "r") as nw_f:
        for line in nw_f:
            newick.append(line)

    newick_str = "".join(newick)
    for key in seqid2name.keys():
        newick_str = newick_str.replace(key, seqid2name.get(key))

    with open(nw_file, "w") as nw_f:
        nw_f.write(newick_str)

    return

def decode_dist_matrix(seqid2name, dist_mat_file, full_dist_mat_file):
    mat_cont = []
    with open(dist_mat_file, "r") as mat_f:
        for line in mat_f:
            first_word = line.split(" ")[0]
            try:
                mat_cont.append(line.replace(first_word, seqid2name[first_word]))
            except KeyError:
                mat_cont.append(line)

    with open(full_dist_mat_file, "w") as op:
        print("".join(mat_cont), file=op)

    return

# Start time to keep track of progress
t0 = time.time()
ctime = int(t0)

# Parse command line options
parser = argparse.ArgumentParser(
    description='Calculates genetic distance between a set of input sequences')
parser.add_argument("-seq", "--seqlist", dest="seq_file", help="Read non-homologous sequences from")
parser.add_argument("-pre", "--prefix", dest="prefix", help="Output path prefix")
parser.add_argument("-a", "--allcalled", dest="allcalled", action="store_true", help="Removal of all uncertain positions")
parser.add_argument("-s", "--probable", dest="probabilistic", action="store_true", help="Probabilistic removal of uncertain positions")
parser.add_argument('-k', '--keep', dest="keep", action="store_true", help='Keep temporary subdirectory')
parser.add_argument("-q", "--quiet", dest="quiet", action="store_true", help="Quiet")
args = parser.parse_args()

suffix = ""


# File with consensus
if args.seq_file is None:
    exiting('No input filelist was provided')
wdir = os.path.dirname(args.seq_file)
if not wdir:
    wdir = os.getcwd()

# Output prefix
prefix = os.path.realpath(args.seq_file)
if args.prefix is not None:
    prefix = args.prefix


# Define nucleotides as numbers
nuc2num = {
   # A  adenosine       C  cytidine          G  guanine
   # T  thymidine       N  A/G/C/T (any)     U  uridine
   # K  G/T (keto)      S  G/C (strong)      Y  T/C (pyrimidine)
   # M  A/C (amino)     W  A/T (weak)        R  G/A (purine)
   # B  G/T/C           D  G/A/T             H  A/C/T
   # V  G/C/A           -  gap of indeterminate length
   #   A C G T
   # A A M R W
   # C M C S Y
   # G R S G K
   # T W Y K T
   'A' : 1,
   'T' : 2,
   'C' : 3,
   'G' : 4,
   'M' : 5,
   'R' : 6,
   'W' : 7,
   'S' : 8,
   'Y' : 9,
   'K' : 10,
   'B' : 11,
   'V' : 12,
   'D' : 13
}


# homology reduced old isolates
inseqs = []
try:
    f = open(args.seq_file, "r")
except IOError as e:
    exiting("Sequence file list not found.")
for l in f:
    l = l.strip()
    inseqs.append(l)
f.close()

timing("# Read inputfiles from the file lists.")

#TODO DONE
slens = len(inseqs) # Number of strains in each list
inputseqmat = None
tot_len = None

# load and encode isolates
# might not need to be parallel (overhead)
if slens > 30:
    arrays = Parallel(n_jobs=no_jobs)(delayed(read_encode_univ)(isolatefile, None) for isolatefile in inseqs)
    # dump as a memmap and concat
    tot_len = np.shape(arrays[0])[0]
    for i, arr in enumerate(arrays):
        if i == 0:
            inputseqmat = np.zeros((slens, tot_len), dtype=np.int8)
            inputseqmat[0,:] = arr[:]
        else:
            inputseqmat[i,:] = arr[:]
    del arrays[:]
else:
    for i, isolatefile in enumerate(inseqs):
        if i == 0:
            tmp_np = read_encode_univ(isolatefile, tot_len)
            tot_len = np.shape(tmp_np)[0]
            inputseqmat = np.zeros((slens, tot_len), dtype=np.int8)
            inputseqmat[0,:] = tmp_np[:]
        else:
            inputseqmat[i,:] = read_encode_univ(isolatefile, tot_len)[:]

timing("# Loaded sequence matrix into memory.")

## Remove uncertain regions
# determine masks, probabilistic, remove: False
m = None
if args.allcalled:
    m = (inputseqmat != 0).all(axis=0)
elif args.probabilistic:
    known_frac = np.round_((inputseqmat != 0).sum(0) / vol_len, 5)
    std_dev = np.std(known_frac, dtype=np.float64)
    mean_a = np.mean(known_frac, dtype=np.float64)
    threshold = mean_a - std_dev
    m = (known_frac >= threshold)
else:
    pass

if m is not None:
    inputseqmat = inputseqmat.T[m].T

# conserved positions for both type of distance calculation, if there is only one volume (ie limited set size)
cons_m = None
if slens > 1:
    cons_m = np.logical_not(np.all(inputseqmat == inputseqmat[0,:], axis = 0))
    inputseqmat = inputseqmat.T[cons_m].T
del cons_m

# update lengths and isolate counts
slens = len(inseqs)
tot_len = np.shape(inputseqmat)[1]

timing("# Removed non-informative positions from matrix.")
# print(inseqs)


if not args.quiet:
    print("# Total length: %s" % (tot_len), file=sys.stdout)
    print("# Number of strains: %s" % (slens), file=sys.stdout)

# calculate genetic distance between isolates
no_jobs = min(slens, no_jobs)
batch_size = int(math.ceil(slens/float(no_jobs)))
dist_calc_func = dist_calc_pw
if args.allcalled:
    dist_calc_func = dist_calc_all

if slens > 1:
    dist_arr = []
    dist_arr = Parallel(n_jobs=no_jobs)(delayed(dist_calc_func)(inputseqmat, inputseqmat, i, batch_size, slens) for i in range(0,slens,batch_size))

    # put it together, do a np.where on it
    dist_o_o = dist_arr[0]
    for np_arr in dist_arr[1:]:
        dist_o_o = np.concatenate((dist_o_o,np_arr), axis=1)

    matrix = dist_o_o.tolist()

else: # only one old seq, dist 0
    matrix = [[0]]

# print(matrix)
timing("# Calculated distances between sequences.")

# print dist matrix in phylip
outputmat = "{}.mat".format(prefix)
seqid2name = {}
seq_id = ""
seqnames = [os.path.basename(x).split(".")[0] for x in inseqs]
with open(outputmat, "w") as matfile:
    print("  {0}".format(len(matrix)), file=matfile)
    for r, row in enumerate(matrix):
        seq_id = "I{:06}".format(r+1)
        seqid2name[seq_id] = seqnames[r]
        print("{:<10}".format(seq_id), end = "", file=matfile)
        for e in row[:-1]:
            print('{0:.0f}'.format(e), end = "\t", file=matfile)
        print('{0:.0f}'.format(row[-1]), file=matfile)

timing("# Constructed distance matrix.")

timing("# Distance calculation is finished.")

mode = 'pw'
if args.allcalled:
    mode = 'all'
elif args.probabilistic:
    mode = 'prob'

method = "dist"
subwdir = change2subdir(mode, method, suffix, wdir)
matpath = os.path.relpath(outputmat)
treefilename = "outtree"

inputstr = "{0}\nL\nY\n".format(matpath)
proc = subprocess.Popen(DTREE, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
proc.stdin.write(inputstr.encode())
# wait for it. it gets stuck in simple wait() if fails
(stdoutdata, stderrdata) = proc.communicate()
if proc.returncode:
    exiting("Neighbor program failed.")

timing("# Distance based tree constructed.")

decode_nw_file(seqid2name, treefilename)
decode_dist_matrix(seqid2name, outputmat, outputmat)

try:
    shutil.copy(treefilename, "{}.nw".format(prefix))
except:
    pass

os.chdir(wdir)
if not args.keep:
    try:
        shutil.rmtree(subwdir)
    except:
        pass

print("Done.", file=sys.stderr)
sys.exit(0)
