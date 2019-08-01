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
from collections import Counter
from dire import die


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = sys.argv[1:]
    assert args
    file = args[0]

    hdrs = ('qseqid sseqid pident length mismatch gapopen qstart '
            'qend sstart send evalue bitscore').split()
    df = pd.read_csv(file, sep='\t', names=hdrs)
    evalue = 1
    wanted = df[df.apply(lambda x: '_' in x['sseqid'] and x['evalue'] < evalue,
                         axis=1)]

    c = Counter(wanted['qseqid'])
    for qseqid, count in c.items():
        if count > 1:
            best = sorted(wanted[wanted['qseqid'] == qseqid]['evalue'])[0]
            indexes = wanted.index[(wanted['qseqid'] == qseqid)
                                   & (wanted['evalue'] > best)].tolist()
            wanted = wanted.drop(indexes)

    print(wanted)

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
