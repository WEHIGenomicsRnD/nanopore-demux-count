#!/usr/bin/env python

import pandas as pd
import sys
from io import StringIO
from natsort import index_natsorted, order_by_index
from argparse import ArgumentParser


def parse_args():
    '''Parse arguments'''
    description = '''
        Merge files for dict based sequence counts

        Usage:
            merge_custom_counts.py <count_files>

        Outputs a table of merged counts
        '''
    parser = ArgumentParser(description=description)
    parser.add_argument('count_files',
                        nargs='+',
                        type=str,
                        help='Count files to merge')

    return parser.parse_args()


def merge_files(count_file,df_param):

    print(f"File - {count_file}")
    try:
       result_df = pd.read_csv(count_file, sep='\t', header=None)

       res=result_df.sort_values(by=3, ascending=False)
       if len(df_param) == 0:
           df_param= res
       else:
           df_param = pd.concat([df_param,res], ignore_index=True)
    except pd.errors.EmptyDataError:
        df_param= df_param


    return df_param

def main():

    filt_outfile=open('merged_filtered_counts.txt','w')
    all_outfile=open('merged_all_counts.txt','w')
    args = parse_args()

    filt_df = pd.DataFrame()
    raw_df = pd.DataFrame()
    for count_file in args.count_files:
       if 'unmapped' in count_file:
           continue

       if 'filtered' in count_file:
          filt_df = merge_files(count_file, filt_df)
       elif 'ref' in count_file:
          raw_df = merge_files(count_file , raw_df)

# write to file
    if len(raw_df) !=0:
       filt_df.columns = ['Sequence','Forward_idx','Reverse_idx','Count','SampleName']
       raw_df.columns = ['Sequence','Forward_idx','Reverse_idx','Count','SampleName']

       filt_df.to_csv(filt_outfile, sep='\t', index=False)
       raw_df.to_csv(all_outfile, sep='\t', index=False)

       filt_outfile.close()
       all_outfile.close()

if __name__ == "__main__":
    main()
