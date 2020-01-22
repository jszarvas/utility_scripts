#!/usr/bin/env python3.6

import sys, time
import os
import math
import numpy as np
import argparse
from operator import itemgetter
import gzip
from joblib import Parallel, delayed

NO_JOBS = 20

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
   # VALIDATE INPUT
   if not isinstance(filename, str):
      msg = 'Filename has to be a string.'
      if exit_on_err: sys.exit('Error: '+msg)
      else: raise IOError(msg)
   if not os.path.exists(filename):
      msg = 'File "%s" does not exist.'%filename
      if exit_on_err: sys.exit('Error: '+msg)
      else: raise IOError(msg)

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
    t1 = time.time()
    print("{0} Time used: {1} seconds".format(message, int(t1-t0)), file=sys.stdout)
    return

def exiting(message):
    print(message, file=sys.stderr)
    print("FAIL", file=sys.stderr)
    sys.exit(1)


def create_vcf(ref, ind, mask):
    """Create VCF file from redundant sequences."""
    #ref: path, vcfname: seqname from path

    date = time.strftime("%Y%m%d")
    version = "VCFv4.2"
    source = sys.argv[0].split("/")[-1]

    reference = refseqs[ref[1]]
    reffilename = refseqlist[ref[1]]
    inputseq = queryseqs[ind]
    inputseqname = querylist[ind].split(".")[0]
    nchrom = len(reference)
    clens = [len(chrom) for chrom in reference]
    nseq = len(inputseq)

    if not args.allcalled:
        mask = [np.logical_and(np.asarray(mask[0][j][ref[1]]), np.asarray(mask[1][j][ind])) for j in range(nchrom)]
    #print(mask)
    vcffile = []
    vcffile.append("##fileformat={0}".format(version))
    vcffile.append("##source={0}".format(source))
    vcffile.append("##reference={0}".format(os.path.join(args.odir, reffilename)))
    for i in range(nchrom):
        vcffile.append("##contig=<ID={0},length={1}>".format(i+1, clens[i]))
    vcffile.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

    for j in range(nchrom):
        for k in range(clens[j]):
            if reference[j][k] != inputseq[j][k] and mask[j][k]:
                if args.context:
                    for c in range(-5,0,1):
                        vcffile.append("{0}\t{1}\t.\t{2}\t{3}\t.\tPASS\t.".format("Chromosome", k+1+c, reference[j][k+c], inputseq[j][k+c]))
                vcffile.append("{0}\t{1}\t.\t{2}\t{3}\t.\tPASS\t.".format("Chromosome", k+1, reference[j][k], inputseq[j][k]))
                if args.context:
                    for c in range(1,6):
                        vcffile.append("{0}\t{1}\t.\t{2}\t{3}\t.\tPASS\t.".format("Chromosome", k+1+c, reference[j][k+c], inputseq[j][k+c]))
    if args.context:
        inputseqname += ".context"
    with open(os.path.join(args.odir, inputseqname + ".vcf"), "w") as of:
        print("\n".join(vcffile), file=of)
    return

def encode_seq_pw(strain):
   encodedinput = [np.zeros((clen), dtype=np.int8) for clen in clens]
   non_nuc_mask = [np.ones((clen), dtype=np.bool) for clen in clens]
   for j in range(len(clens)):
      for i in range(clens[j]):
         try:
            encodedinput[j][i] = nuc2num[strain[j][i]]
         except KeyError:
            non_nuc_mask[j][i] = False
   print(np.sum(non_nuc_mask[0]))
   return encodedinput, non_nuc_mask

def encode_seq_all(strain):
   encodedinput = [np.zeros((clen), dtype=np.int8) for clen in clens]
   non_nuc_mask = [np.ones((clen), dtype=np.bool) for clen in clens]
   for j in range(len(clens)):
      for i in range(clens[j]):
         try:
            encodedinput[j][i] = nuc2num[strain[j][i]]
         except KeyError:
            non_nuc_mask[j][i] = False
   return encodedinput, non_nuc_mask

# Start time to keep track of progress
t0 = time.time()
etta = 0.001

# Parse command line options
parser = argparse.ArgumentParser(
    description='Creates VCF files from reference and query file pairs')
parser.add_argument("-rl", dest="old", help="list of reference sequences", metavar="INFILE_HR")
parser.add_argument("-ql", dest="new", help="list of query sequences", metavar="INFILE_NEW")
parser.add_argument("-r", dest="ref", help="reference")
parser.add_argument("-q", dest="query", help="query")
parser.add_argument("-o", "--outputdir", dest="odir", help="write to DIR", metavar="DIR")
parser.add_argument("-c", "--context", dest="context", action="store_true", help="Print 5 pos around snp")
parser.add_argument("-a", "--allcalled", dest="allcalled", action="store_true", help="Only use positions called in all strains")
args = parser.parse_args()


# File with reads
listbased = False
if args.ref is None or args.query is None:
    if args.old is None or args.new is None:
        exiting('No input file(list) was provided')
    else:
        listbased = True

# path to results dir ../results_db/template/
if args.odir is None:
    exiting('Output directory is needed')

# New strains
querylist = []
queryseqs = []

if listbased:
   with open(args.new) as f:
      for l in f:
         l = l.strip()
         if l == "":
           next #avoid empty line
         fp = os.path.join(args.odir, l)
         if os.path.exists(fp) and os.path.getsize(fp) > 0:
            entries = list(zip(*[[seq, name, desc] for seq, name, desc in SeqsFromFile(fp)]))
            querylist.append(l)
            queryseqs.append(list(entries[0])) #tuple with number of chromosomes
else:
   entries = list(zip(*[[seq, name, desc] for seq, name, desc in SeqsFromFile(args.query)]))
   querylist.append(os.path.basename(args.query))
   queryseqs.append(list(entries[0]))

refseqs = []
refseqlist = []

if listbased:
   with open(args.old) as f:
      for l in f:
         # First strain should be the the template
         l = l.strip()
         if l == "":
            next #avoid empty line
         fp = os.path.join(args.odir, l)
         if os.path.exists(fp):
            refseqlist.append(l)
            entries = list(zip(*[[seq, name, desc] for seq, name, desc in SeqsFromFile(fp)]))
            refseqs.append(list(entries[0]))
else:
   entries = list(zip(*[[seq, name, desc] for seq, name, desc in SeqsFromFile(args.ref)]))
   refseqlist.append(os.path.basename(args.ref))
   refseqs.append(list(entries[0]))

timing("# Read inputfiles from the file lists.")

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
   'K' : 10
}

clens = [len(chrom) for chrom in refseqs[0]] # Length of chromosomes
tot_len = sum(clens) #total length of sequences, for the pw calc
nchrom = len(clens)
slens = [len(refseqs), len(queryseqs)] # Number of strains in each list
mask_lens = clens
if not False:
    print("# Length of chromosomes %s" %(clens), file=sys.stdout)
    print("# Number of strains: %s" % (slens), file=sys.stdout)

# Adjust number of processes for memory consumption
no_jobs = min(slens[1], NO_JOBS)
refseqsmat = []
queryseqsmat = []
non_nuc_masks = []
mask_nw = []
mask_hr = []

arrays = []
if args.allcalled:
    for inputseqmat, inputseq, mask in ([refseqsmat, refseqs, mask_hr], [queryseqsmat, queryseqs, mask_nw]):

        arrays = Parallel(n_jobs=no_jobs)(delayed(encode_seq_all)(isolate) for isolate in inputseq)

        inputseqmat.extend([np.asarray([item[0][i] for item in arrays]) for i in range(nchrom)])
        mask.extend([np.asarray([item[1][i] for item in arrays]).all(axis=0) for i in range(nchrom)]) # take AND for mask matrix by evaling if all True
    non_nuc_masks = [np.logical_and(np.asarray(mask_hr[i]), np.asarray(mask_nw[i])) for i in range(nchrom)] # take AND for the hr and nw chrom arrays

else:
    for inputseqmat, inputseq in ([refseqsmat, refseqs], [queryseqsmat, queryseqs]):

        arrays = Parallel(n_jobs=no_jobs)(delayed(encode_seq_pw)(isolate) for isolate in inputseq)

        inputseqmat.extend([np.asarray([item[0][i] for item in arrays]) for i in range(nchrom)])
        non_nuc_masks.append([np.asarray([item[1][i] for item in arrays]) for i in range(nchrom)])

timing("# Encoded sequences into numpy array.")

# create the vcf files
for i in range(slens[1]):
    # loop over references if less than query
    ref = (0, i % slens[0])
    create_vcf(ref, i, non_nuc_masks)


timing("# VCF creation is finished.")
print("DONE", file=sys.stderr)
sys.exit(0)
