#!/usr/bin/env python3
"""
Author : kyclark
Date   : 2019-03-29
Purpose: Rock the Casbah
"""

import argparse
import logging
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
                        type=argparse.FileType('r'),
                        help=('Path to a FASTA-formatted file that '
                              'contains a Papillomavirus genome.'),
                        required=True)

    parser.add_argument('-f',
                        '--format',
                        metavar='FORMAT',
                        choices=['fasta', 'genbank'],
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

    args = parser.parse_args()

    if not os.path.isdir(args.data_dir):
        parser.error('--data_dir "{}" is not a directory.'.format(
            args.data_dir))

    return args


# --------------------------------------------------
def main():
    """main"""
    args = get_args()

    level = {
        'debug': logging.DEBUG,
        'info': logging.INFO,
        'warning': logging.WARNING,
        'error': logging.ERROR,
        'critical': logging.CRITICAL
    }
    logging.basicConfig(level=level[args.debug_level],
                        filename=args.log_file,
                        filemode='w')
    #try:
    puma.run(vars(args))
    print("ALL GOOD")
    #except Exception as e:
    #    die('There was a fatal error: {}'.format(e))


# --------------------------------------------------
if __name__ == '__main__':
    main()
