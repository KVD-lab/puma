#!/usr/bin/env python3
"""
Author : kyclark
Date   : 2019-08-19
Purpose: Rock the Casbah
"""

import argparse
import csv
import os
import sys

# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Argparse Python script',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i',
                        '--input',
                        help='Input file',
                        metavar='FILE',
                        type=argparse.FileType('r'),
                        required=True)

    parser.add_argument('-o',
                        '--output',
                        help='Output file',
                        metavar='FILE',
                        type=str,
                        required=True)

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()

    reader = csv.DictReader(args.input, delimiter=',')
    out_fh = open(args.output, 'wt')
    num_written = 0

    for row in reader:
        if row['gene'] == 'E1':
            acc = row['accession']
            seq = row['translated seq']
            if seq.endswith('*'):
                seq = seq[:-1]

            out_fh.write('>{}\n{}\n'.format(acc, seq))
            num_written += 1

    out_fh.close()

    print(f'Done, wrote {num_written} to "{args.output}".')

# --------------------------------------------------
if __name__ == '__main__':
    main()
