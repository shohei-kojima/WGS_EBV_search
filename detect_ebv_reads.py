#!/usr/bin/env python

'''
Author: Shohei Kojima @Keio/RIKEN
Prerequisites:
    Python 3.7 or higher
    samtools (tested with v1.21)
    pysam (tested with v0.22.1)
    pandas (tested with v2.2.2)
    numpy (tested with v1.26.4)
Input file:
    This can only take BAM as an input file. This CANNOT take SAM and CRAM.
    Input BAM file must be sorted and indexed.
Example usage:
    python detect_ebv_reads.py \
    -i sample1.bam \
    -o sample1.tsv \
    -c chrEBV
Others:
    This will only use 1 thread.
'''

import os,sys,argparse,subprocess
import numpy as np
import pandas as pd
import pysam


# arguments
parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', metavar = 'str', type = str, help = 'Input paired-end BAM files containing chrEBV.')
parser.add_argument('-o', metavar = 'str', type = str, help = 'Output file name.')
parser.add_argument('-c', metavar = 'str', type = str, help = 'EBV chromosome name.', default = 'chrEBV')
args = parser.parse_args()


# check python version
py_version = sys.version_info
if (py_version[0] >= 3) and (py_version[1] >= 7):
    print('Python %d.%d.%d' % (py_version[0], py_version[1], py_version[2]))
else:
    print('Please use Python 3.7 or higher.', file = sys.stderr)
    exit(1)


def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    elif 'PATH' in os.environ:
        for path in os.environ['PATH'].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


# check samtools
if which('samtools') is None:
    print('samtools not found in PATH.')
    exit(1)


# make output directry
outdir = os.path.dirname(args.o)
if len(outdir) >= 1:
    os.makedirs(outdir, exist_ok = True)


# settings
chrEBV = args.c
REF_GENOME_LENGTH = 171_823
PHRED_THRESHOLD = 20
MIN_MAPPED_LENGTH = 145
W_REPEATS = [[12000, 15072], [15072, 18144], [18144, 21216], [21216, 24288], [24288, 27360], [27360, 30432], [30432, 33504], [33504, 35355]]
WREP_START = W_REPEATS[0][0]
WREP_END = W_REPEATS[-1][1]
N_MAX_READS = 350_000


def nt_to_i(nt):
    nt = nt.upper()
    if nt == 'A':
        return 0
    if nt == 'T':
        return 1
    if nt == 'G':
        return 2
    if nt == 'C':
        return 3
    return 4

def has_overlap(s1, e1, s2, e2):
    tmp = (e2 - s1) * (s2 - e1)
    return tmp < 0

def count_ebv_reads(f):
    samfile = pysam.AlignmentFile(f, 'rb')
    # get the number of reads mapped to chrEBV
    counts = samfile.get_index_statistics()
    chrebv_found = False
    for c in counts:
        chr = c.contig
        if chr == chrEBV:
            n_total_reads = c.mapped
            chrebv_found = True
            break
    if not chrebv_found:
        print('%s not found in %s.' % (chrEBV, f), file = sys.stderr)
        exit(1)
    # subsample if too many mapped reads found
    if n_total_reads > N_MAX_READS:
        # downsample
        subsample_bam = f + '.subsample.bam'
        subsample = N_MAX_READS / n_total_reads
        cmd = 'samtools view -bh -s %.6f -o %s %s %s' % (subsample, subsample_bam, f, chrEBV)
        subprocess.run(cmd, shell = True)
        f = subsample_bam
    # detect EBV reads
    r1_set = set()
    r2_set = set()
    inside_w_set = set()
    geno = np.zeros((REF_GENOME_LENGTH, 5), dtype = np.int32)
    samfile = pysam.AlignmentFile(f, 'rb')
    for read in samfile:
        read_name = read.query_name
        if not read.is_proper_pair:
            continue
        if read.is_qcfail:
            continue
        if read.is_secondary:
            continue
        if read.is_supplementary:
            continue
        if read.is_duplicate:
            continue
        # check mapped length
        cigartuples = read.cigartuples
        mapped_length = 0
        for op, length in cigartuples:
            if op == 0:
                mapped_length += length
        if mapped_length < MIN_MAPPED_LENGTH:
            continue
        # save read name
        if read.is_read1:
            r1_set.add(read_name)
        elif read.is_read2:
            r2_set.add(read_name)
        # w-repeat check
        ref_s = read.reference_start
        ref_e = read.reference_end
        isin_wrep = has_overlap(ref_s, ref_e, WREP_START, WREP_END)
        if isin_wrep:
            inside_w_set.add(read_name)
    ebv_reads = r1_set & r2_set
    n_inside_w = len(ebv_reads & inside_w_set)
    n_outside_w = len(ebv_reads - inside_w_set)
    # variant calling
    samfile = pysam.AlignmentFile(f, 'rb')
    for read in samfile:
        if not read.query_name in ebv_reads:
            continue
        readseq = read.query_sequence
        qual = read.query_qualities
        if qual is None:
            qual = [33] * 1024
        for readpos, refpos in read.get_aligned_pairs():
            if readpos is None:
                continue
            if refpos is None:
                continue
            q = qual[readpos]
            if q < PHRED_THRESHOLD:
                continue
            nt = readseq[readpos]
            i = nt_to_i(nt)
            geno[refpos, i] += 1
    # estimate EBV read number
    if n_total_reads > N_MAX_READS:
        n_ebv_read = int((n_outside_w + n_inside_w) * n_total_reads / N_MAX_READS)
        os.remove(subsample_bam)
    else:
        n_ebv_read = n_outside_w + n_inside_w
    # summarize information
    df = pd.DataFrame(geno, columns = ['A', 'T', 'G', 'C', 'N'])
    df.index.name = 'POS'
    df = df[ df.sum(axis = 1) >= 1 ]
    non_w_depth = n_outside_w / (171_823 - 35355 + 12000)
    w_depth = n_inside_w / (15072 - 12000)
    if non_w_depth > 0:
        w_rep_count = w_depth / non_w_depth
    else:
        w_rep_count = np.nan
    # comment line in dataframe
    attrs = ['chrEBV_reads', 'EBV_reads', 'reads_outside_W', 'reads_inside_W', 'W_copy_num_estimate', 'total_mapped_length']
    data = [n_total_reads, n_ebv_read, n_outside_w, n_inside_w, w_rep_count, df.shape[0]]
    tmp = []
    for a, d in zip(attrs, data):
        tmp.append('%s=%s' % (a, str(d)))
    comment = '#' + ','.join(tmp)
    return comment, df




# main
comment, df = count_ebv_reads(args.i)
with open(args.o, 'w') as outfile:
    outfile.write(comment + '\n')
    df.to_csv(outfile, sep = '\t')
