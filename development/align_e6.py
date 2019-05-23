#!/usr/bin/env python3
"""
Author : kyclark
Date   : 2019-05-23
Purpose: Rock the Casbah
"""

import argparse
import os
import re
import sys
from pprint import pprint as pp
from dire import die, warn
from Bio import AlignIO


# --------------------------------------------------
def get_args():
    """get command-line arguments"""
    parser = argparse.ArgumentParser(
        description='Argparse Python script',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('fasta', metavar='str', help='FASTA')

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Make a jazz noise here"""
    args = get_args()
    fasta = args.fasta

    if not os.path.isfile(fasta):
        die('"{}" is not a file'.format(fasta))

    alignment = AlignIO.read(fasta, 'fasta')
    print(alignment)
    alignment_length = alignment.get_alignment_length()
    conserved = {}
    dashes = {}
    seq_by_id = dict([(rec.id, str(rec.seq)) for rec in alignment])

    for i, rec in enumerate(alignment):
        seq = str(rec.seq)
        match = re.search(r'^([-]+)', seq)
        num_dashes = 0
        if match:
            dash = match.group(1)
            num_dashes = len(dash)

        dashes[rec.id] = num_dashes

    pp(dashes)

    for i in range(0, alignment_length):
        col = alignment[:, i]
        conserved[i] = col.count('M')

    #pp(conserved)
    max_num_seqs = max(conserved.values())
    print(max_num_seqs)
    seqs_at_max = list(filter(lambda t: t[1] == max_num_seqs, conserved.items()))

    if len(seqs_at_max) > 1:
        warn('Hey, look at this!')

    best_col = seqs_at_max[0][0]
    print(seqs_at_max)
    prefix = seq_by_id['unknown'][0:best_col]

    if len(prefix) * '-' == prefix:
        print('all dashes, do nothing')
    elif prefix.count('-') == 0:
        print('E6 = ', seq_by_id['unknown'][best_col:])
    else:
        print('remove {} from beginning'.format(len(seq.replace('-', ''))))


    # max_num_seqs = max(sorted([v for k, v in conserved.items()])
    # pp(pos)
    # print(max(pos))


# --------------------------------------------------
if __name__ == '__main__':
    main()
