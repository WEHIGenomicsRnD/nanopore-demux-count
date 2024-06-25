#!/usr/bin/env python
'''
Module      : Align sequence pairs
Description : Take two fasta files and perform pairwise-alignment
              on each sequence pair
Copyright   : (c) WEHI Genomics R&D, 2024
License     : TBD
Maintainer  : Marek Cmero
Portability : POSIX
'''

import sys
import os
import edlib
import logging
import pyfastx as fx
from datetime import datetime
from argparse import ArgumentParser

debug = False
BASE_BUFFER_LEN = 50


def parse_args():
    '''Parse arguments'''
    description = '''
        Take two fasta files and perform pairwise-alignment
        on each sequence pair.

        Usage:
            python align_sequence_pairs.py query.fasta target.fasta

        Output:
            Returns aligned sequence pairs to stdout
        '''
    parser = ArgumentParser(description=description)
    parser.add_argument('query_fasta',
                        type=str,
                        help='Input query FASTA file')
    parser.add_argument('target_fasta',
                        type=str,
                        help='Input target FASTA file')
    parser.add_argument('--debug',
                        action='store_true',
                        help='Enable debug logging')

    return parser.parse_args()


def init_logging(log_filename):
    '''
    Initiate logging
    '''
    logging.basicConfig(filename=log_filename,
                        level=logging.DEBUG,
                        filemode='w',
                        format='%(asctime)s %(levelname)s - %(message)s',
                        datefmt="%Y-%m-%dT%H:%M:%S%z")
    logging.info('Align sequence pairs started')


def main():
    '''Main'''
    args = parse_args()

    if not os.path.exists(args.query_fasta) or \
       not os.path.exists(args.target_fasta):
        sys.exit('Error: FASTA file(s) do not exist')

    if args.debug:
        global debug
        debug = True
        init_logging('rotate_{:%Y-%m-%d_%H%M}.log'.format(datetime.now()))

    queries = fx.Fasta(args.query_fasta, build_index=False)
    targets = fx.Fasta(args.target_fasta, build_index=False)

    target_dict = {}
    for target in targets:
        name, seq = target
        target_dict[name] = seq

    for query in queries:
        name, seq = query

        if name not in target_dict:
            logging.warning(f'Query {name} not found in target')
            continue

        target_seq = target_dict[name]
        result = edlib.align(seq, target_seq, mode='HW', task='path')
        result = edlib.getNiceAlignment(result, seq, target_seq)

        print('----------------------------------------')
        print(f'{name}')
        print('----------------------------------------')
        align_len = len(result['query_aligned'])
        for i in range(0, align_len, BASE_BUFFER_LEN):
            query_seq = result['query_aligned'][i:i+BASE_BUFFER_LEN]
            matched_seq = result['matched_aligned'][i:i+BASE_BUFFER_LEN]
            target_seq = result['target_aligned'][i:i+BASE_BUFFER_LEN]
            print(f'{name}_QUERY\t{query_seq}')
            print(f'{name}_MATCH\t{matched_seq}')
            print(f'{name}_TARGT\t{target_seq}\n')

    if debug:
        logging.info('Align sequence pairs finished')


if __name__ == '__main__':
    main()
