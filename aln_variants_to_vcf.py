#!/usr/bin/env python3

import sys
import os
import argparse
import time

class SAMPLE_ALN:
    """Sample alignment to VCF"""

    def __init__(self, filepath):
        self.path = filepath
        self.samplename = os.path.basename(filepath).split(".")[0]

    def read_aln(self):
        template = []
        aln = []
        subject = []
        with open(self.path, "r") as aln_p:
            for line in aln_p:
                line = line.strip("\n")
                if not line:
                    continue

                tmp = line.split("\t")
                if tmp[0].startswith("#"):
                    continue

                if tmp[0].startswith("template"):
                    template.append(tmp[1])
                elif tmp[0].isspace():
                    aln.append(tmp[1])
                elif tmp[0].startswith("query"):
                    subject.append(tmp[1])
                else:
                    pass
        self.template = "".join(template)
        self.aln = "".join(aln)
        self.subject = "".join(subject)
        self.templ_len = len(self.template) - self.template.count("-")
        return

    def coverage_check(self):
        # simplify problem: number of gaps cant exceed 1% of total length
        no_gaps = self.subject.count("-")
        if no_gaps / len(self.subject) > 0.03:
            return False
        else:
            return True

    def create_vcf_header(self, ref):
        """Create VCF file header"""
        #ref: path, vcfname: seqname from path

        date = time.strftime("%Y%m%d")
        version = "VCFv4.2"
        source = sys.argv[0].split("/")[-1]

        # obs: can't process multiple alignment/chromosomes
        self.reference = ref
        nchrom = 1
        clens = self.templ_len

        self.vcf_content = []
        self.vcf_content.append("##fileformat={0}".format(version))
        self.vcf_content.append("##source={0}".format(source))
        self.vcf_content.append("##reference={0}".format(ref))
        self.vcf_content.append("##contig=<ID={0},length={1}>".format(ref, clens))
        self.vcf_content.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

        return

    def add_vcf_line(self, pos, alt, ref):
        """Create VCF file from alignment"""

        chrom = self.reference
        self.vcf_content.append("{0}\t{1}\t.\t{2}\t{3}\t.\tPASS\t.".format(chrom, pos+1, ref, alt))

        return

    def process_aln(self, gff_cds, gff_mask):
        aln_len = len(self.template)
        # inserts start on -1, dels on 0
        self.indel_start = None
        self.indel_type = None
        # compensate for inserts
        self.shift = 0
        mask_i = 0
        for i in range(aln_len):
            if gff_mask and not self.template[i] == "-":
                ref_pos = i - self.shift
                if ref_pos >= gff_mask[mask_i][0] and ref_pos <= gff_mask[mask_i][1]:
                    #print(i, ref_pos, self.shift)
                    if self.indel_type is not None and self.indel_type != "mask":
                        # interrupting a deletion
                        if self.indel_type == "deletion":
                            self.indel_type = "part_deletion"
                            self.process_variant(i)
                        # after an insertion
                        elif self.indel_type == "insert":
                            self.indel_type = "err_insert"
                            self.process_variant(i)
                            self.shift += i - self.indel_start - 1
                    if self.indel_type != "mask":
                        self.indel_type = "mask"
                        self.indel_start = i
                    continue

                elif ref_pos > gff_mask[mask_i][1]:
                    if mask_i + 1 != len(gff_mask):
                        mask_i +=1
                    if self.indel_type == "mask":
                        self.process_variant(i)
                        self.indel_type = None
                        self.indel_start = None

            if self.aln[i] == "_":
                if self.template[i] == "-":
                    # insertion
                    if self.indel_start is None:
                        self.indel_start = i - 1
                        self.indel_type = "insert"
                    elif self.indel_type == "deletion":
                        # not a deletion any more, switch type
                        self.process_variant(i)
                        self.indel_start = i - 1
                        self.indel_type = "insert"
                    elif self.indel_type == "mask":
                        # insert in a masked region
                        self.process_variant(i)
                        self.indel_start = i - 1
                        self.indel_type = "insert"
                    else:
                        pass

                elif self.subject[i] == "-":
                    # deletion
                    if self.indel_start is None:
                        self.indel_start = i
                        self.indel_type = "deletion"
                    elif self.indel_type == "insert":
                        # not an insert any more, switch type
                        self.process_variant(i)
                        self.shift += i - self.indel_start - 1
                        self.indel_start = i
                        self.indel_type = "deletion"
                    else:
                        pass
                else:
                    if self.indel_start is not None:
                        # indel ended on prev. pos
                        self.process_variant(i)
                        if self.indel_type == "insert":
                            self.shift += i - self.indel_start - 1
                        self.indel_type = None
                        self.indel_start = None
                    # currently mismatch
                    self.indel_type = "mismatch"
                    self.indel_start = i
                    self.process_variant(i)
                    self.indel_type = None
            else:
                # not a mismatch
                if self.indel_start is not None:
                    # indel ended on prev. pos
                    self.process_variant(i)
                    if self.indel_type == "insert":
                        self.shift += i - self.indel_start - 1
                    self.indel_type = None
                    self.indel_start = None
        if self.indel_start is not None:
            # end of string
            self.process_variant(i + 1)
        return

    def process_variant(self, pos):
        frmshft = False
        indel_len = 0
        ref_pos = self.indel_start - self.shift

        if self.indel_type == "insert":
            indel_len = pos - self.indel_start - 1
        else:
            indel_len = pos - self.indel_start

        if self.indel_type in ["insert", "deletion"]:
            frmshft = frameshift_check(ref_pos, indel_len)
            if frmshft:
                self.indel_type = "err_" + self.indel_type

        if self.indel_type == "part_deletion":
            self.indel_type = "deletion"

        if self.subject[self.indel_start].lower() == "n":
            self.indel_type = None

        if self.indel_type == "mismatch":
            self.add_vcf_line(ref_pos, self.subject[self.indel_start], self.template[self.indel_start])
        elif self.indel_type == "deletion":
            # previous nucleotide is shown
            self.add_vcf_line(ref_pos-1, self.template[self.indel_start-1] , self.template[self.indel_start-1:pos])
        elif self.indel_type == "insert":
            self.add_vcf_line(ref_pos, self.subject[self.indel_start:pos], self.template[self.indel_start])
        else:
            pass

        return

def frameshift_check(ref_pos, indel_len):
    frmshft = False
    for cds in gff_cds:
        # starts within the CDS
        if ref_pos >= cds[0] and ref_pos <= cds[1] and (indel_len % 3) != 0:
            frmshft = True
    return frmshft

def read_gff(gff_f, feature):
    feature_pos = []
    reference = ""
    with open(gff_f, "r") as gff_p:
        for line in gff_p:
            if not line.startswith("#"):
                cols = line.strip().split()
                if cols and cols[2] == feature:
                    if not reference:
                        reference = cols[0]
                    # go to 0 based pos
                    start = int(cols[3]) - 1
                    end = int(cols[4]) - 1
                    feature_pos.append([start, end])
    return feature_pos, reference

############# MAIN #############
# Parse command line options
parser = argparse.ArgumentParser(
    description='Creates VCF files from reference and query file pairs')
parser.add_argument("-i", dest="aln_file", help="Alignment file")
parser.add_argument("-r", dest="ref", help="reference gff")
parser.add_argument("-m", dest="mask", help="mask gff")
parser.add_argument("-o", dest="odir", help="write to DIR")
args = parser.parse_args()

# read CDS positions from gff file
gff_cds = []
if os.path.exists(args.ref) and os.path.getsize(args.ref) > 0:
    gff_cds, reference = read_gff(args.ref, "CDS")

# get masked positions
gff_mask = []
if os.path.exists(args.mask) and os.path.getsize(args.mask) > 0:
    gff_mask, _ = read_gff(args.mask, "mask")

if os.path.exists(args.aln_file) and os.path.getsize(args.aln_file) > 0:
    sample = SAMPLE_ALN(args.aln_file)
    sample.read_aln()
    if not sample.coverage_check():
        # check number of gaps in query
        print("#", sample.samplename, "too short")
        continue
    sample.create_vcf_header(reference)
    sample.process_aln(gff_cds, gff_mask)

    with open("{}.vcf".format(os.path.join(args.odir, sample.samplename)), "w") as of:
        print("\n".join(sample.vcf_content), file=of)
else:
    print("Error: File not found {}".format(args.aln_file), file=sys.stderr)
