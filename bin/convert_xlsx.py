#!/usr/bin/env python

import pandas as pd
import sys
from argparse import ArgumentParser

def parse_args():
    '''Parse arguments'''
    description = '''
        Process input excel file for primer ans index

        Usage:
            convert_xlsx.py <excel_file>

        Outputs a primers and indexes
        '''
    parser = ArgumentParser(description=description)
    parser.add_argument('-i',
                        '--excel_file',
                        type=str,
                        help='Input excel file')

    return parser.parse_args()


def main():

    args = parse_args()

    print(f"Filename - {args.excel_file}")
    df = pd.ExcelFile(args.excel_file)
    
    if len(df.sheet_names) > 1:

       index_out=open("index.csv",'w')
       primer_out=open("primer.csv",'w')

       primers = pd.read_excel(args.excel_file, sheet_name=0)
       index =   pd.read_excel(args.excel_file, sheet_name=1)

       primers['name'] = ['primer'+str(n+1) for n in range(primers.shape[0])]
#    index['direction'] = index['direction'].apply(lambda x: 'Fwd' if x.startswith('F') else 'Rev')

       primers.to_csv(primer_out, sep=',', index=False)
       index.to_csv(index_out, sep=',', index=False)
    else:
       print(f"WARN:Excel file should have primers and index")
    


    if len(df.sheet_names) > 2:
        insert_out=open("guides.fa",'w')

        insert = pd.read_excel(args.excel_file, sheet_name=2)
        print(f"Insert- {insert}")
        for i, row in insert.iterrows():
           insert_out.write(">" + str(row['name']) + '\n')
           insert_out.write(str(row['sequence']) + '\n')
        

if __name__ == "__main__":
    main()
