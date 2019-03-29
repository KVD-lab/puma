import puma
from Bio import SeqIO

# --------------------------------------------------
def test_foo():
    foo = puma.foo()
    assert len(foo) == 1
    assert foo == ['foo']

# --------------------------------------------------
def test_trans_orf():
    seq = SeqIO.sequence('foo')
    res = puma.trans_orf(seq, -1, 1)
    assert res == ''
