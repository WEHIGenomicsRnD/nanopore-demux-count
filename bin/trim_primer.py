#!/usr/bin/env python
'''
Module      : trim_primer
Description : Locates forward and reverse primers in FASTQ reads and trims excess sequence
Copyright   : (c) WEHI Genomics R&D, 2024
License     : TBD
Maintainer  : Marek Cmero
Portability : POSIX
'''
import gzip
import edlib
import os
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from argparse import ArgumentParser

tab = str.maketrans("ACTG", "TGAC")

def rc(seq):
    # reverse complement sequence
    # https://bioinformatics.stackexchange.com/questions/3583
    return seq.translate(tab)[::-1]

def parse_args():
    '''Parse arguments'''
    description = '''
        Trim read by primer sequences.
        '''
    parser = ArgumentParser(description=description)
    parser.add_argument('reads',
                        metavar='READS',
                        type=str,
                        help='Reads in fastq format.')
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
    parser.add_argument('-b',
                        '--barcode_length',
                        metavar='BARCODE_LENGTH',
                        default=13,
                        type=int,
                        help='Barcode length.')
    parser.add_argument('--untrimmed_fastq',
                        type=str,
                        default='',
                        help='Output file where to write untrimmed reads (optional).')
 
    return parser.parse_args()

def trim_read(read_seq, read_qual, fwd_primer, rev_primer, mismatches, barcode_length, rc_read=False):
    read_seq = rc(read_seq) if rc_read else read_seq
    read_qual = read_qual[::-1] if rc_read else read_qual
    fwd_result = edlib.align(fwd_primer, read_seq, mode="HW", task="path", k=mismatches)
    rev_result = edlib.align(rev_primer, read_seq, mode="HW", task="path", k=mismatches)

    if fwd_result["locations"] and rev_result["locations"]:
        start = fwd_result["locations"][0][0] - barcode_length
        end = rev_result["locations"][0][1] + barcode_length

        if start < 0 or end > len(read_seq):
            # cannot trim at correct barcode position
            return None, None

        return (read_seq[start:(end+1)], read_qual[start:(end+1)])

    return None, None

def main():
    args = parse_args()

    if args.untrimmed_fastq.strip() != "":
        untrimmed_out = gzip.open(args.untrimmed_fastq, "wt")
    else:
        untrimmed_out = None

    if args.reads.endswith(".gz"):
        in_handle = gzip.open(args.reads, "rt")
    else:
        in_handle = open(args.reads, "r")

    trimmed_reads, untrimmed_reads, total_reads = 0, 0, 0
    with in_handle:
        for read_id, read_seq, read_qual in FastqGeneralIterator(in_handle):
            total_reads += 1
            trimmed_seq, trimmed_qal = trim_read(read_seq, read_qual, args.fwd_primer, args.rev_primer, args.mismatches, args.barcode_length)
            if trimmed_seq:
                trimmed_reads += 1
                print("@%s\n%s\n+\n%s" % (read_id, trimmed_seq, trimmed_qal))
                continue

            # try with reverse complement
            trimmed_seq, trimmed_qal = trim_read(read_seq, read_qual, args.fwd_primer, args.rev_primer, args.mismatches, args.barcode_length, rc_read=True)
            if trimmed_seq:
                trimmed_reads += 1
                print("@%s\n%s\n+\n%s" % (read_id, trimmed_seq, trimmed_qal))
            elif untrimmed_out:
                # couldn't trim
                untrimmed_reads += 1
                untrimmed_out.write("@%s\n%s\n+\n%s\n" % (read_id, read_seq, read_qual))

    # print stats to file
    sample_name = os.path.basename(args.reads).split(".")[0]
    stats_file_name = f"{sample_name}_trimming_stats.txt"
    percent_trimmed = round((trimmed_reads/total_reads)*100, 2)
    percent_untrimmed = round((untrimmed_reads/total_reads)*100, 2) if untrimmed_out else 0

    with open(stats_file_name, "w") as stats_file:
        stats_file.write(f"Trimmed reads: {trimmed_reads} ({percent_trimmed}%)\n")
        if untrimmed_out:
            stats_file.write(f"Untrimmed reads: {untrimmed_reads} ({percent_untrimmed}%)\n")
        stats_file.write("Total reads: %d\n" % total_reads)
        stats_file.write("Mismatches: %d\n" % args.mismatches)
        stats_file.write("Barcode length: %d\n" % args.barcode_length)

    if untrimmed_out:
        untrimmed_out.close()

if __name__ == "__main__":
    main()
