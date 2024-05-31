#!/usr/bin/env python
'''
Module      : collate_counts
Description : Collate multiple count files into one sorted table
Copyright   : (c) WEHI Genomics R&D, 2024
License     : TBD
Maintainer  : Marek Cmero
Portability : POSIX
'''
import pandas as pd
import sys
from io import StringIO
from natsort import index_natsorted, order_by_index
from argparse import ArgumentParser


def parse_args():
    '''Parse arguments'''
    description = '''
        Count guide sequences

        Usage:
            collate_counts.py <count_files>

        Outputs a table of merged counts
        '''
    parser = ArgumentParser(description=description)
    parser.add_argument('count_files',
                        nargs='+',
                        type=str,
                        help='Count files to merge')

    return parser.parse_args()


def main():
    args = parse_args()

    result_df = pd.DataFrame()
    for count_file in args.count_files:
        df = pd.read_csv(count_file, sep='\t')

        if len(result_df) == 0:
            result_df = df
        else:
            result_df = pd.merge(result_df, df, how='outer', on='guide')

    new_index = order_by_index(result_df.index, index_natsorted(result_df['guide']))
    result_df = result_df.reindex(index=new_index)

    # write to stdout
    output = StringIO()
    result_df.to_csv(output, sep='\t', index=False)
    output.seek(0)
    print(output.read(), file=sys.stdout)


if __name__ == "__main__":
    main()
