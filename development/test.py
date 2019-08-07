import puma
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import BiopythonWarning
import warnings

# --------------------------------------------------
def test_make_l1_end():
    warnings.simplefilter('ignore', BiopythonWarning)
    #HPV16
    l1Result = [5639, 7156]
    orig_file = os.path.join('data_dir', 'tests', 'original_genome.txt')
    assert os.path.isfile(orig_file)
    original_genome = open(orig_file).read().strip()

    new_file = os.path.join('data_dir', 'tests', 'linearized_genome.txt')
    assert os.path.isfile(new_file)
    new_genome = open(new_file).read().strip()
    assert new_genome

    res = puma.make_l1_end(l1Result, original_genome)
    assert res == new_genome
# --------------------------------------------------
# def test_linearize_genome():
#     # HPV16
#     warnings.simplefilter('ignore', BiopythonWarning)
#     orig_file = os.path.join('data_dir', 'tests', 'original_genome.txt')
#     assert os.path.isfile(orig_file)
#     original_genome = open(orig_file).read().strip()
#
#     args = { 'out_dir': './puma_out', 'data_dir': 'data_dir',
#         'input_format': 'genbank', 'min_prot_len': 25,
#         'e_value': 1e-05, 'for_user_dir': './puma_out/for_user',
#         'program_files_dir': './puma_out/program_files'}
#
#     new_file = os.path.join('data_dir', 'tests', 'linearized_genome.txt')
#     assert os.path.isfile(new_file)
#     new_genome = open(new_file).read().strip()
#
#     res = puma.linearize_genome(original_genome,args)
#     assert res == new_genome
# --------------------------------------------------
def test_trans_orf():
    warnings.simplefilter('ignore', BiopythonWarning)
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
def test_run_blastp():
    warnings.simplefilter('ignore', BiopythonWarning)

    orfs = os.path.join('data_dir', 'tests', 'orfs.fa')
    assert os.path.isfile(orfs)
    blast_sub = os.path.join('data_dir', 'tests', 'main_blast.fa')
    assert os.path.isfile(blast_sub)
    blast_out = os.path.join('data_dir', 'tests', 'blast_results_main.tab')
    assert os.path.isfile(blast_out)
    res = puma.run_blastp(orfs,blast_sub,blast_out)
    assert res == 6



# --------------------------------------------------
#def test_blast_main_orfs():


# --------------------------------------------------
#def test_identify_main_proteins():

# --------------------------------------------------
#def test_verify_l1():

# --------------------------------------------------
#def test_blast_e5_variants():

# --------------------------------------------------
#def test_identify_e5_variants():


# --------------------------------------------------
#def test_blast_verify_e6():


# --------------------------------------------------
#def test_parse_blast_results_verify_e6():

# --------------------------------------------------
#def test_align_verify_e6():

# --------------------------------------------------
#def test_verify_e6():

# --------------------------------------------------
#def test_find_urr():

# --------------------------------------------------
#def test_fimo_e1bs():

# --------------------------------------------------
#def test_find_e1bs():

# --------------------------------------------------
#def test_fimo_e2bs():

# --------------------------------------------------
#def test_find_e2bs():

# --------------------------------------------------
#def test_blast_splice_acceptor():

# --------------------------------------------------
#def test_locate_known_splice_acceptor():

# --------------------------------------------------
#def test_align_splice_acceptor():

# --------------------------------------------------
#def test_find_splice_acceptor():

# --------------------------------------------------
#def test_blast_spliced_e1_e8():

# --------------------------------------------------
#def test_locate_known_e1_splice_donor():

# --------------------------------------------------
#def test_align_splice_donor_e1():

# --------------------------------------------------
#def test_find_e1_e4():

# --------------------------------------------------
#def test_locate_known_e8_splice_donor():

# --------------------------------------------------
#def test_align_splice_donor_e8():

# --------------------------------------------------
#def test_find_e8():

# --------------------------------------------------
#def test_find_e8_e2():

# --------------------------------------------------
#def test_validate_args():












