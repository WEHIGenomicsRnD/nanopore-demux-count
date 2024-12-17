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


def merge_files(cfile,result_df,rfile):
   df = pd.read_csv(cfile, sep='\t')
   if len(result_df) == 0:
      result_df = df
   else:
      onval='guide'
      if rfile:
          onval=['guide','type']
      result_df = pd.merge(result_df, df, how='outer', on=onval) 

   return result_df

def main():
    args = parse_args()

    out=open("collated_overall.txt",'w')
    result_df = pd.DataFrame()
    overall_df= pd.DataFrame()
    for count_file in args.count_files:
        if "overall" in count_file:
            overall_df=merge_files(count_file,overall_df,1)
        else:
            result_df=merge_files(count_file,result_df,0)
           

    new_index = order_by_index(result_df.index, index_natsorted(result_df['guide']))
    new_index1 = order_by_index(overall_df.index, index_natsorted(zip(overall_df.guide, overall_df.type)))
    result_df = result_df.reindex(index=new_index)
    overall_df = overall_df.reindex(index=new_index1)

    overall_df.to_csv(out, sep='\t', index=False)
    # write to stdout
    output = StringIO()
    result_df.to_csv(output, sep='\t', index=False)
    output.seek(0)
    print(output.read(), file=sys.stdout)


if __name__ == "__main__":
    main()
