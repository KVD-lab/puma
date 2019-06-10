#!/usr/bin/env python3
"""
Author : kyclark
Date   : 2019-03-29
Purpose: Rock the Casbah
"""

import argparse
import re
import os
import sys
import puma


# --------------------------------------------------
def warn(msg):
    """Print a message to STDERR"""
    print(msg, file=sys.stderr)


# --------------------------------------------------
def die(msg='Something bad happened'):
    """warn() and exit with error"""
    warn(msg)
    sys.exit(1)


# --------------------------------------------------
def get_args():
    args = sys.argv
    bin_dir = os.path.dirname(args[0])

    parser = argparse.ArgumentParser(
        description='HPC annotator'
        'information within a given papillomavirus '
        'genome.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i',
                        '--input',
                        metavar='FILE',
                        # type=argparse.FileType('r'),
                        help=('Path to a FASTA-formatted file that '
                              'contains a Papillomavirus genome.'),
                        required=True)

    parser.add_argument('-f',
                        '--format',
                        metavar='FORMAT',
                        choices=['fasta'],
                        default='fasta',
                        help='File format')

    parser.add_argument('-d',
                        '--data_dir',
                        metavar='DIR',
                        type=str,
                        default=os.path.join(bin_dir, 'data_dir'),
                        help='Directory that has all database files for use')

    parser.add_argument('-o',
                        '--outdir',
                        metavar='DIR',
                        type=str,
                        default=os.path.join(bin_dir, 'puma-out'),
                        help='Output directory')

    parser.add_argument('-g',
                        '--gff3',
                        action='store_true',
                        help='Outputs information in gff3 format')

    parser.add_argument('-c',
                        '--csv',
                        action='store_true',
                        help='Outputs information in csv format')

    parser.add_argument('-e',
                        '--evalue',
                        metavar='FLOAT',
                        type=float,
                        default=0.00001,
                        help='BLAST evalue')

    parser.add_argument('-m',
                        '--min_prot_len',
                        metavar='NUM',
                        type=int,
                        default=25,
                        help='Minimum protein length')

    parser.add_argument(
        '-s',
        '--sites',
        metavar='STR',
        type=str,
        default='ALL',
        help='Comma-separated string of L1 L2 E1 E2 E4 E5 E5_delta '
        'E5_zeta E5_epsilon E6 E7 '
        'E9 E10 '
        'E2BS E1BS URR ALL')

    parser.add_argument(
        '-D',
        '--debug_level',
        metavar='STR',
        type=str,
        default='critical',
        choices=['debug', 'info', 'warning', 'error', 'critical'],
        help='Debug level')

    parser.add_argument('-L',
                        '--log_file',
                        metavar='FILE',
                        type=str,
                        default='run.log',
                        help='Debug log file')

    #return parser.parse_known_args()
    return parser.parse_args()


# --------------------------------------------------
def main():
    """main"""
    args = vars(get_args())  # make into dict

    args['sites'] = re.split(r'\s*,\s*', args['sites'].upper())

    # input_file = args.input
    # out_dir = args.outdir
    # data_dir = args.data_dir
    # input_format = args.format.lower()
    # min_prot_len = args.min_prot_len
    # evalue = args.evalue
    # gff3 = args.gff3
    # csv = args.csv
    # blastE1E8_dir = os.path.join(out_dir, 'blastE1E8')

    # puma.run(args)
    try:
        puma.run(args)
        print("ALL GOOD")
    except Exception as e:
        die('There was a fatal error: {}'.format(e))


# --------------------------------------------------
if __name__ == '__main__':
    main()
