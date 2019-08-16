#!/usr/bin/env python3

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC

fa = list(SeqIO.parse('foo.fa', 'fasta'))[0]
raw = str(fa.seq)[6077:]
unambig = Seq(str(fa.seq), IUPAC.unambiguous_dna)[6077:]
ambig = Seq(str(fa.seq))[6077:]

#print(unambig, file=open('unambig.txt', 'wt'))
#print(ambig, file=open('ambig.txt', 'wt'))
#print('ambig', ambig)
#print(seq)
#print(type(seq))
#print(len(seq[6077:]))

#print(ambig.translate())
import re
print(re.findall('[^ATCG]', raw.upper(), re.IGNORECASE))
print(Seq(raw.upper(), IUPAC.unambiguous_dna).translate())
