#!/usr/bin/env python
'''
Module      : trim_primer
Description : Locates forward and reverse primers in FASTQ reads and trims excess sequence
Copyright   : (c) WEHI Genomics R&D, 2024
License     : TBD
Maintainer  : Marek Cmero
Portability : POSIX
'''
import logging
import os
import gzip
import edlib
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from datetime import datetime
from argparse import ArgumentParser

tab = str.maketrans("ACTG", "TGAC")

def init_log(filename):
    '''
    Initialise logging
    '''
    logging.basicConfig(filename=filename,
                        level=logging.DEBUG,
                        filemode='w',
                        format='%(asctime)s %(levelname)s - %(message)s',
                        datefmt="%Y-%m-%dT%H:%M:%S%z")

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
    parser.add_argument('--quiet',
                        action='store_true',
                        help='Do not write log file.')
 
    return parser.parse_args()

def trim_read(read_id, read_seq, fwd_primer, rev_primer, mismatches, barcode_length, logfile):
    fwd_result = edlib.align(fwd_primer, read_seq, mode="HW", task="path", k=mismatches)
    rev_result = edlib.align(rev_primer, read_seq, mode="HW", task="path", k=mismatches)

    if fwd_result["locations"] and rev_result["locations"]:
        start = fwd_result["locations"][0][0] - barcode_length
        end = rev_result["locations"][0][1] + barcode_length

        if start < 0 or end > len(read_seq):
            if logfile: logging.debug("Failed to trim read %s (barcode partially missing)", read_id)
            return None

        return read_seq[start:end]

    elif logfile:
        logging.debug("Failed to trim read %s (primers not found)", read_id)

    return None

def main():
    args = parse_args()

    if args.quiet:
        logfile = None
    else:
        logfile = datetime.now().strftime('runlog_%H_%M_%d_%m_%Y.log')
        init_log(logfile)

    if args.reads.endswith(".gz"):
        in_handle = gzip.open(args.reads, "rt")
    else:
        in_handle = open(args.reads, "r")

    with in_handle:
        for read_id, read_seq, read_qual in FastqGeneralIterator(in_handle):
            trimmed_seq = trim_read(read_id, read_seq, args.fwd_primer, args.rev_primer, args.mismatches, args.barcode_length, logfile)
            if trimmed_seq:
                print("@%s\n%s\n+\n%s\n" % (read_id, trimmed_seq, read_qual))
                continue

            # try with reverse complement
            trimmed_seq = trim_read(read_id, rc(read_seq), args.fwd_primer, args.rev_primer, args.mismatches, args.barcode_length, logfile)
            if trimmed_seq:
                print("@%s\n%s\n+\n%s\n" % (read_id, trimmed_seq, read_qual))

if __name__ == "__main__":
    main()