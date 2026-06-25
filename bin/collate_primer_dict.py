#!/usr/bin/python

import os
import argparse
from argparse import ArgumentParser

def parse_args():
    '''Parse arguments'''
    description = '''
        Merge Dict based count files for all primers

        Usage:
            collate_stats.py files>

        Outputs a table of merged count file
        Fasta file for each primer group
        '''
    parser = ArgumentParser(description=description)
    parser.add_argument('files',
                        nargs='+',
                        type=str,
                        help='Primer count files to merge')

    return parser.parse_args()


#GGAGCCGAATGGCTCCTCCCTGAAGAACATAGAGAAGTATCTCAGAAGTCAAAGTGATCTCACAAGCACCACCAACAACCCAGCCTTTCAGCAGCGGCTGCGACTGGGGGCCAAACGCGCTGTGAATAATGGGAGGTTACTGAAAGACGGACCGCAGTACAGGGTCAATTAATGGGAGCTTAGATGGCAAAGGGGCACCTCAGTATCCCAGTGCATTCCCATCCTCGCTCCCACCTGTCAGCCTTCTACCCCATGA  Fwd_05  Rev_08  1923    barcode03
#insert  FwdBarcodeID    RevBarcodeID    TotalForSample  AmpliconPrimer  SampleID        TotalForSequence        percent username
#CGACACTCACTTCCGCACCTTCCGCTCCCACTCCGATTACCGGCGCATCACGCGGACCAGCGCGCTCCTGGACGCCTGCGGCTTCTATTGGGGACCCCTGAGCGTGCACGGGGCGCACGAGCGGCTGCGTGCCGAGCCCGTGGGCACCTTCTTGGTGGCCGACAGTCGCCAACGGAACTGCTTCTTCGCGCTCAGCGTGAAGATGGCTTCGGGCCCC Fwd_11  Rev_03  306     CGACACTCAC      Fwd_11-Rev_03-CGACACTCAC        70      0.22875816993464052     evelyn


def main():

    args = parse_args()

    out=open("AllPrimer_stats.txt",'w')
    out.write("insert\tFwdBarcodeID\tRevBarcodeID\tAmpliconPrimer\tSampleID\tTotalForSequence\n")

    for f in args.files:
        primer = f.replace(".merged_filtered_counts.txt", "")
        out_fa=open(f"{primer}.fa",'w')
        with open(f) as fh:
            for line in fh:
               if line.startswith("Sequence"):
                  continue

               lsplit=line.split("\t")
               out_fa.write(f">{primer}:sampleID:{lsplit[1]}-{lsplit[2]}-{primer}:count:{lsplit[3]}\n")
               out_fa.write(f"{lsplit[0]}\n")
               out.write(f"{lsplit[0]}\t{lsplit[1]}\t{lsplit[2]}\t{primer}\t{lsplit[1]}-{lsplit[2]}-{primer}\t{lsplit[3]}\n")
            out_fa.close()



if __name__ == "__main__":
    main()
