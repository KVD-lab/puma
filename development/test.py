import puma
import os

from Bio import SeqIO
from Bio.Seq import Seq

# --------------------------------------------------
def test_trans_orf():

    new_file = os.path.join('data_dir', 'tests', 'linearized_genome.txt')
    assert os.path.isfile(new_file)
    new_genome = open(new_file).read().strip()
    assert new_genome

    res = puma.trans_orf(Seq(new_genome), 25)

    assert isinstance(res, dict)
    assert len(res) == 38
    assert 'ISVVCCVYVCMCLYVLVNIKLYVCLYVWYNKHVCMCF' in res
    assert res['ISVVCCVYVCMCLYVLVNIKLYVCLYVWYNKHVCMCF'] == 18


# --------------------------------------------------
def test_make_l1_end():
    """Test"""

    #HPV16
    l1Result = [5639, 7156]
    orig_file = os.path.join('data_dir', 'tests', 'original_genome.txt')
    assert os.path.isfile(orig_file)
    original_genome = open(orig_file).read().strip()

    new_file = os.path.join('data_dir', 'tests', 'linearized_genome.txt')
    assert os.path.isfile(new_file)
    new_genome = open(new_file).read().strip()
    assert new_genome

    res = puma.make_l1_end(l1Result, original_genome, len(original_genome))
    assert res == new_genome
