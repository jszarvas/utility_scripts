#!/usr/bin/env python3

import sys, time
import os
import math
import numpy as np
import argparse
import gzip
import subprocess
import multiprocessing
from joblib import Parallel, delayed
import string
from copy import deepcopy

def combine_dict_array(dict_array):
    '''Combine/merge an array of dicts into a single dict by summing the values in the dict'''
    combined_dict = deepcopy(dict_array[0])
    for cont in dict_array[1:]:
        for bc_hash in cont.keys():
            if bc_hash in combined_dict:
                combined_dict[bc_hash] += cont.get(bc_hash)
            else:
                combined_dict[bc_hash] = cont.get(bc_hash)
    return combined_dict

def split_whitelist_into_8mers(filename):
    '''Split the 24bps long whitelisted barcodes into three 8-mers'''
    eightmers = [[],[],[]]
    bound = [0,8,16]
    with open(filename, "r") as fp:
        for line in fp:
            tmp = line.strip()
            for i in range(3):
                eightmers[i].append(tmp[bound[i]:bound[i]+8])
    return eightmers

def bc_quality_check(barcode, mer_boundaries = [0,8,16]):
    '''Return True if there are no more than one Ns in a k-mer'''
    i = 0
    no_mers = len(mer_boundaries)
    k = mer_boundaries[1] - mer_boundaries[0]
    while i < no_mers:
        if barcode[mer_boundaries[i]:mer_boundaries[i]+k].find("N", barcode[mer_boundaries[i]:mer_boundaries[i]+k].find("N", 0) + 1) > -1:
            return False
        else:
            i += 1
    return True

def read_R3_bc_offset(filename, no_bcs = 500000, skip = 0, offset = 0):
    '''Parse gzip-ed fastq of barcodes from an offset and by skipping skip amount of reads'''
    bcs_strings = []
    offset_lines = offset * 4
    line_count = 0
    max_lines = no_bcs * 4 * (skip + 1) + offset_lines
    with gzip.open(filename, "rt") as fp:
        for line in fp:
            line_count += 1
            if line_count < offset_lines:
                continue
            elif line_count > max_lines:
                break
            else:
                if (line_count % (4 * (skip + 1))) == 2:
                    line = line.strip()
                    if bc_quality_check(line):
                        bcs_strings.append(line)
    return bcs_strings

def reduce_eightmer(eightmers):
    '''Return a list of unique k-mers'''
    unique_mers = []
    for i in range(3):
        unique_mers.append(list(set(eightmers[i])))
    return unique_mers

def encode_eightmer(umer_list, k_size = 8):
    '''Encode k-mers in an numpy array'''
    arrays = []
    for i in range(3):
        arrays.append(np.zeros((len(umer_list[i]), k_size), dtype=np.int8))
        for j in range(len(umer_list[i])):
            for k in range(8):
                # triplet_i : 8mer_j : letter_k
                arrays[i][j,k] = letter_conversion.get(umer_list[i][j][k])
    return arrays

def encode_raw_bc(bc_string, len_umers, k_size = 8):
    '''Encode a raw barcode string into a numpy array'''
    arrays = []
    for i in range(3):
        single_arr = np.zeros((1,k_size), dtype=np.int8)
        for k in range(k_size):
            single_arr[0,k] = letter_conversion.get(bc_string[i*k_size+k])
        arrays.append(np.tile(single_arr, (len_umers,1)))
    return arrays

def decode_original_bc(wcs, ca, cb, cc):
    '''Return the corresponding whitelisted barcode from the combined umers'''
    # from array indeces it can go straight to umers as the order is preserved
    a = wcs[0].index(ca)
    b = wcs[1].index(cb)
    c = wcs[2].index(cc)
    wl_index = a + b + c
    return wcs[0][wl_index] + wcs[1][wl_index] + wcs[2][wl_index]

def clustering_wcs_bc(wc_arr, bc_strings, dist_thr = 1):
    '''Cluster the raw barcodes to the whitelisted ones'''
    wcs_distro = {}
    for bcs in bc_strings:
        string_arr = encode_raw_bc(bcs, 96)
        dist_t = np.zeros(shape=(3,96), dtype=np.uint8)
        for i in range(3):
            for k in range(96):
                dist_t[i,k] = np.not_equal(wc_arr[i][k,], string_arr[i][k,]).sum(0)
        closest_wc = np.where(dist_t < dist_thr + 1)
        if closest_wc[0].shape[0] == 3 and closest_wc[0].sum() == 3:
            wcs_hash = ";".join([str(x) for x in closest_wc[1]])
            if wcs_hash in wcs_distro:
                wcs_distro[wcs_hash] += 1
            else:
                wcs_distro[wcs_hash] = 1
    return wcs_distro

def get_abundant_wcs(wcs_clusters, min_count, whitelist, uniqmers):
    '''Return whitelisted barcodes above min. abundance'''
    abundant_wcs = []
    for index_hash, count in wcs_clusters.items():
        if count > min_count:
            ind = [int(x) for x in index_hash.split(";")]
            wbc = decode_original_bc(whitelist, uniqmers[0][ind[0]], uniqmers[1][ind[1]], uniqmers[2][ind[2]])
            # print(wbc, count)
            abundant_wcs.append(wbc)
    return abundant_wcs

def weave_fastqs(fastqs):
    '''Read fastq into streams'''
    no_files = len(fastqs)
    processes = [subprocess.Popen("gzip --stdout -d {}".format(fn), shell=True, text=True, stdout=subprocess.PIPE) for fn in fastqs]
    streams = [r.stdout for r in processes]

    try:
        while True:
            identifiers = [next(s)[:-1] for s in streams]
            names = [x.split()[0] for x in identifiers]
            seqs = [next(s)[:-1] for s in streams]
            blanks = [next(s)[:-1]  for s in streams]
            quals = [next(s)[:-1]  for s in streams]
            assert all(name==names[0] for name in names)
            reads = [identifiers[i]+'\n'+seqs[i]+'\n+\n'+quals[i]+'\n' for i in range(no_files)]
            yield seqs, reads
    except StopIteration:
        pass

    for s in streams:
        s.close()

def identify_wcs_bc(wc_arr, orphan_bc_list, dist_thr = 1):
    '''Cluster the raw barcodes that did not match exactly to the whitelisted barcodes'''
    identified_bcs = {}
    no_mers = len(wc_arr)
    for [bc_i, bc] in orphan_bc_list:
        # assuming equal amount of mers in each part
        string_arr = encode_raw_bc(bc, wc_arr[0].shape[0])
        dist_t = np.zeros(shape=(no_mers,wc_arr[0].shape[0]), dtype=np.uint8)
        for i in range(no_mers):
            for k in range(wc_arr[0].shape[0]):
                dist_t[i,k] = np.not_equal(wc_arr[i][k,], string_arr[i][k,]).sum(0)
        closest_wc = np.where(dist_t < dist_thr + 1)
        if closest_wc[0].shape[0] == 3 and closest_wc[0].sum() == 3:
            wcs_hash = ";".join([str(x) for x in closest_wc[1]])
            if wcs_hash in identified_bcs:
                identified_bcs[wcs_hash].append(bc_i)
            else:
                identified_bcs[wcs_hash] = [bc_i]
    return identified_bcs

def output_fastqgzip(name):
    for i in range(no_ends):
        with gzip.open(os.path.join(args.odir, f"{name}_R{i+1}.fastq.gz"), "wt") as op:
            op.write("".join(split_libraries[i][name_to_bc.get(name)]))
    return

parser = argparse.ArgumentParser(
    description='De-multiplex split-and-pool barcoded paired-end fastqs')
parser.add_argument(
    '-i',
    dest="fastq_files",
    required=True,
    help='Comma separated read fastq path(s)')
parser.add_argument(
    '-b',
    dest="bc_file",
    required=True,
    help='Barcode fastq path')
parser.add_argument(
    '-o',
    dest="odir",
    required=True,
    help='Output folder for the demultiplexed fastqs')
parser.add_argument(
    '-w',
    dest="whitelist",
    required=True,
    help='List of whitelisted barcodes')
parser.add_argument(
    '--bc',
    dest="sample_barcodes",
    help="Path to file of the found barcodes")
parser.add_argument(
    '-n',
    dest="raw_bc_num",
    default=2000000,
    type=int,
    help='Number of initially sampled raw barcode sequences')
parser.add_argument(
    '-l',
    dest="offset",
    default=10,
    type=int,
    help='Size of the offsets between batches in million reads')
parser.add_argument(
    '-t',
    dest="parallel",
    default=4,
    type=int,
    help='Number of parallel processes')
parser.add_argument(
    '-c',
    dest="min_size",
    default=10,
    type=int,
    help='Minimum read count in initial sampling per barcode to be accepted')
parser.add_argument(
    '-m',
    dest="min_read_count",
    default=1000,
    type=int,
    help='Minimum read count per barcode to write to disk')
args = parser.parse_args()

t0 = time.time()

# load whitelisted barcodes
letter_conversion = {'A': 0, 'T': 1, 'C': 2, 'G': 3, 'N': -1}
wcs = split_whitelist_into_8mers(args.whitelist)
umers = reduce_eightmer(wcs)
umer_arrays = encode_eightmer(umers)

# get list of barcodes present in sequenced sample
present_wcs = []
if args.sample_barcodes is None or not os.path.exists(args.sample_barcodes):
    # predict the barcodes present in the sequencing run
    if os.path.exists(args.bc_file):
        # read in raw bcs
        r3_bcs = []
        bcs_batch = int(args.raw_bc_num / args.parallel)
        read_spacing = math.floor(args.offset * 10 ** 6 / (bcs_batch + 1))
        for i in range(args.parallel):
            offset_batch = args.offset * 10 ** 6 * (i + 1)
            r3_bcs.append(read_R3_bc_offset(args.bc_file, no_bcs = bcs_batch, skip = read_spacing, offset = offset_batch))
        wcs_distro_list = Parallel(n_jobs=args.parallel)(delayed(clustering_wcs_bc)(umer_arrays, bc_list) for bc_list in r3_bcs)
        combined_distro = combine_dict_array(wcs_distro_list)
        present_wcs = sorted(get_abundant_wcs(combined_distro, args.min_size, wcs, umers))
        if args.sample_barcodes is not None:
            with open(args.sample_barcodes, "w") as of:
                for item in present_wcs:
                    print(item, file=of)
                print("Whitelisted barcodes present written to file", file=sys.stderr)
        del r3_bcs
        del combined_distro
    else:
        sys.exit("Please provide the raw barcodes.")
else:
    # use the previously determined list of barcodes to demultiplex
    if os.path.exists(args.sample_barcodes):
        with open(args.sample_barcodes, "r") as fp:
            for line in fp:
                present_wcs.append(line.strip())
    else:
        sys.exit("Sample barcodes file not found.")

if not present_wcs:
    sys.exit("Barcodes not found in sample")
else:
    print("Collected whitelisted barcodes present", file=sys.stderr)

# demultiplex whole fastqs using the barcodes found
input_files = args.fastq_files.split(",")
no_ends = len(input_files)
split_libraries = []
for i in range(no_ends):
    split_libraries.append({barcode: [] for barcode in present_wcs})
input_files = input_files + [args.bc_file]
orphan_bcs = []
orphan_reads = [[] for x in range(no_ends)]
cnt = 0
for seqs, reads in weave_fastqs(input_files):
    if bc_quality_check(seqs[-1]):
        bc = seqs[-1]
        if bc in present_wcs:
            # append reads except last
            for i in range(no_ends):
                split_libraries[i][bc].append(reads[i])
        else:
            orphan_bcs.append([cnt, bc])
            for i in range(no_ends):
                orphan_reads[i].append(reads[i])
            cnt += 1

print("Exact-matching reads gathered", file=sys.stderr)

# calculate the nearest wcs to the non-exact-matching barcodes
if orphan_bcs:
    batch_size = math.ceil(len(orphan_bcs) / args.parallel)
    iwc_array = Parallel(n_jobs=args.parallel)(delayed(identify_wcs_bc)(umer_arrays, orphan_bcs[i:i+batch_size]) for i in range(0,len(orphan_bcs), batch_size))
    identified_reads = combine_dict_array(iwc_array)
    for wc_hash, indeces in identified_reads.items():
        ind = [int(x) for x in wc_hash.split(";")]
        wbc = decode_original_bc(wcs, umers[0][ind[0]], umers[1][ind[1]], umers[2][ind[2]])
        if wbc in split_libraries[0]:
            # we ignore the wcs not found already
            for i in range(no_ends):
                for j in indeces:
                    split_libraries[i][wbc].append(orphan_reads[i][j])

del orphan_reads
del orphan_bcs

print("Mis-matching reads gathered", file=sys.stderr)

# give each barcode a unique name
az_list = list(string.ascii_lowercase)
code4 = [a+b+c+d for a in az_list for b in az_list for c in az_list for d in az_list]
code4 = code4[:len(present_wcs)]
name_to_bc = {i:j for i,j in zip(code4,present_wcs)}

read_counts_of = open(os.path.join(args.odir, 'read_counts.txt'),'w')
fastqs_to_output = []
for name,bc in name_to_bc.items():
    total_read_count=sum([len(split_libraries[i][bc]) for i in range(no_ends)])
    print(f"{name}\t{total_read_count}", file = read_counts_of)

    if total_read_count > args.min_read_count:
        fastqs_to_output.append(name)
read_counts_of.close()

print("Read counts written", file=sys.stderr)

if fastqs_to_output:
    if __name__ == '__main__':

        p = multiprocessing.Pool(8)
        p.imap_unordered(output_fastqgzip, fastqs_to_output)
        p.close()
        p.join()

print("Time to finish:", time.time()-t0)
