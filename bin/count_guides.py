#!/usr/bin/env python
'''
Module      : count_guides
Description : Count guide sequences from aligned bam,
              ensuring that read spans the whole guide
Copyright   : (c) WEHI Genomics R&D, 2024
License     : TBD
Maintainer  : Marek Cmero
Portability : POSIX
'''
import os
import sys
import pysam
from Bio import SeqIO
from argparse import ArgumentParser
import collections

def parse_args():
    '''Parse arguments'''
    description = '''
        Count guide sequences

        Usage:
            count_guides.py <bam> <guide_reference>

        Outputs a count table of guide sequences from aligned bam.
        '''
    parser = ArgumentParser(description=description)
    parser.add_argument('bam',
                        type=str,
                        help='Reads in fastq format.')
    parser.add_argument('guide_reference',
                        type=str,
                        help='Guide reference fasta file.')
    parser.add_argument('sample',
                        type=str,
                        help='Sample name.')
    parser.add_argument('--lenient',
                        action='store_true',
                        help='Count partial mappings.')

    return parser.parse_args()


def main():
    args = parse_args()

    if not os.path.exists(args.bam):
        print(f'{args.bam} does not exist', file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.guide_reference):
        print(f'{args.guide_reference} does not exist', file=sys.stderr)
        sys.exit(1)

    sample_name = args.sample + "_" + os.path.basename(args.bam).split(".")[0]
    outfile=f"{sample_name}_overall.txt"
    out=open(outfile,'w')
    out.write("guide\ttype\t"+sample_name+"\n") 
    # create a reference look up for guie lengths
    guide_lens = {}
    with open(args.guide_reference) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            guide_lens[record.id] = len(record.seq)

    ref_partial=collections.defaultdict(dict)
    counts = {'unmapped': 0, 'partial_map': 0}
    for guide in guide_lens:
        counts[guide] = 0
        ref_partial[guide]['lenient'] = 0
        ref_partial[guide]['strict'] = 0

    bamfile = pysam.AlignmentFile(args.bam, 'r')
    for read in bamfile:
#        print(read)
        if read.is_unmapped:
            counts['unmapped'] += 1
            continue

        if read.reference_name not in guide_lens:
            print(f'Guide {read.reference_name} not found in reference',
                  file=sys.stderr)
            continue

        guide_len = guide_lens[read.reference_name]
        bam_len = read.reference_end - read.reference_start
        if bam_len != guide_len:
                ref_partial[read.reference_name]['lenient'] += 1
        else:
            ref_partial[read.reference_name]['strict'] += 1


        read_spans_ref = read.reference_end - read.reference_start == guide_len
        if args.lenient or read_spans_ref:
            counts[read.reference_name] += 1
        else:
            counts['partial_map'] += 1
            print(f'Guide {read.reference_name} does not span the whole guide',
                  file=sys.stderr)

    sample_name = args.sample + "_" + os.path.basename(args.bam).split(".")[0]
    print(f'guide\t{sample_name}')
    for guide, count in counts.items():
        print(f'{guide}\t{count}')

    for g,c in ref_partial.items():
        for k,v in c.items():
           print(f'{k}- {g} - {ref_partial[g][k]}',file=sys.stderr)
           out.write(g+"\t"+k+"\t"+str(ref_partial[g][k])+"\n")

if __name__ == "__main__":
    main()
