#!/usr/bin/env python3
"""
Author : kyclark
Date   : 2019-06-12
Purpose: Rock the Casbah
"""

import argparse
import os
import sys
import pandas as pd
from dire import die

# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    hdrs = ('qseqid sseqid pident length mismatch gapopen qstart '
            'qend sstart send evalue bitscore').split()
    df = pd.read_csv('blast_results_E5_HPV16.tab', sep='\t', names=hdrs)
    evalue = 1
    wanted = df[df.apply(lambda x: '_' in x['sseqid'] and x['evalue'] < evalue,
                         axis=1)]

    if len(wanted) < 1:
        die('Nobody loves you')

    d = dict(zip(wanted['sseqid'], wanted['qseqid']))
    print(d)

    #d = {}
    #for i, row in wanted.iterrows():
    #    d[row['sseqid']] = row['qseqid']
    #print(d)

# --------------------------------------------------
if __name__ == '__main__':
    main()
