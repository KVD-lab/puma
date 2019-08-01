#!/usr/bin/env python3

import os
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from dire import die


def run_blastp(query, subject, outfile, evalue=1e-4):
    cmd = blastp(
        query=query,
        subject=subject,
        evalue=evalue,
        outfmt=6,
        out=outfile)

    print(cmd)
    stdout, stderr = cmd()

    #if stderr:
    #    logging.warn("STDERR = ", stderr)

    if not os.path.isfile(outfile):
        raise Exception('No BLAST output')

    if os.path.getsize(outfile) == 0:
        return 0

    return len(open(outfile).read().splitlines())

def main():
    orfs = 'orfs_E5_HPV6.fa'
    e5 = 'blast_E5.fa'
    outfile = 'e5_results.tab'

    try:
        num_hits = run_blastp(orfs, e5, outfile, evalue=1e-10)
        print('Got {} hits!'.format(num_hits))
    except Exception as e:
        die('EXCEPTION: {}'.format(e))

    print('All was good!')

if __name__ == '__main__':
    main()
