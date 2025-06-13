#!/usr/bin/env python

import gzip
import edlib
import os
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from argparse import ArgumentParser
from collections import defaultdict

def parse_args():
    '''Parse arguments'''
    description = '''
        Custom reference free counts.
        '''
    parser = ArgumentParser(description=description)
    parser.add_argument('-r',
                        '--reads',
                        metavar='READS',
                        type=str,
                        help='Trimmed reads in fastq format.')
    parser.add_argument('-f',
                        '--fwd_primer',
                        metavar='FWD_PRIMER',
                        type=str,
                        help='Forward primer sequence.')
    parser.add_argument('-v',
                        '--rev_primer',
                        metavar='REV_PRIMER',
                        type=str,
                        help='Reverse primer sequence.')
    parser.add_argument('-m',
                        '--mismatches',
                        metavar='MISMATCHES',
                        default=3,
                        type=int,
                        help='Number of mismatches to allow.')
    parser.add_argument('-s',
                        '--sample',
                        metavar='SAMPLE',
                        type=str,
                        help='Sample Name')
 
    return parser.parse_args()

ref_dict={}
def fetch_seq(read_seq, read_qual, fwd_primer, rev_primer, mismatches):
    fwd_result = edlib.align(fwd_primer, read_seq, mode="HW", task="path", k=mismatches)
    rev_result = edlib.align(rev_primer, read_seq, mode="HW", task="path", k=mismatches)

#    print(f"Read : {read_seq} -- {fwd_result["locations"]} -- {rev_result["locations"]} ")
    if fwd_result["locations"] and rev_result["locations"]:
        start = fwd_result["locations"][0][1] + 1
        end = rev_result["locations"][0][0]

        if start < 0 or end > len(read_seq):
            # cannot fetch sequence at correct barcode position
            return None

        return (read_seq[start:(end)])

    return None



def parse_ref(file):
   with open(file) as f:
      for line in f:
         l=line.strip().split("\t")
         ref_dict[l[1]]=[l[0],l[2]]
   
   f.close()
#   return ref_dict
   
def main():
    args = parse_args()
    seq_count = {}
    seq_num=1
    if args.reads.endswith(".gz"):
        in_handle = gzip.open(args.reads, "rt")
    else:
        in_handle = open(args.reads, "r")

#    parse_ref(args.guide)
#    print(ref_seq)

    trimmed_reads, untrimmed_reads, total_reads = 0, 0, 0
    with in_handle:
        for read_id, read_seq, read_qual in FastqGeneralIterator(in_handle):
            total_reads += 1
            trimmed_seq = fetch_seq(read_seq, read_qual, args.fwd_primer, args.rev_primer, args.mismatches)
            if trimmed_seq:
               try:
                  seq_count[trimmed_seq] = seq_count[trimmed_seq] +1
               except KeyError as e:
                      seq_count[trimmed_seq] =  1
#               continue

    # print sequence count to file
    sample_name = os.path.basename(args.reads).split(".")[0]
    print(f"Sample - {sample_name}")
    (fwd_idx, rev_idx) = sample_name.split("-") if '-' in sample_name else ("F","R")

    file_name = f"{args.sample}-{sample_name}_ref-free_count.txt"
    filt_fname = f"{args.sample}-{sample_name}_filtered_count.txt"
    ref_file=open(file_name, "w")
#    filt_file=open(filt_fname, "w")

    for k, v in seq_count.items():
       ref_file.write("%s\t%s\t%s\t%d\t%s\n" %(k,fwd_idx, rev_idx ,seq_count[k],args.sample))
       if seq_count[k] >= 25:
          with open(filt_fname, "a") as filt_file:
             filt_file.write("%s\t%s\t%s\t%d\t%s\n" %(k,fwd_idx, rev_idx ,seq_count[k],args.sample))

#    ref_file.close()
#    filt_file.close()


if __name__ == "__main__":
    main()
