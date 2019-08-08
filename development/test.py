"""
Tests are currently written with HPV16
"""

import puma
import os
import csv

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import BiopythonWarning
import warnings


# --------------------------------------------------
def test_make_l1_end():
    """Docstring"""

    warnings.simplefilter('ignore', BiopythonWarning)
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
def test_trans_orf():
    """Docstring"""

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
    """Docstring"""

    warnings.simplefilter('ignore', BiopythonWarning)
    orfs = os.path.join('data_dir', 'tests', 'orfs.fa')
    assert os.path.isfile(orfs)
    blast_sub = os.path.join('data_dir', 'tests', 'main_blast.fa')
    assert os.path.isfile(blast_sub)
    blast_out = os.path.join('data_dir', 'tests', 'blast_results_main.tab')
    assert os.path.isfile(blast_out)
    res = puma.run_blastp(orfs, blast_sub, blast_out)
    assert res == 6


# --------------------------------------------------
def test_linearize_genome():
    """Docstring"""

    warnings.simplefilter('ignore', BiopythonWarning)
    orig_file = os.path.join('data_dir', 'tests', 'original_genome.txt')
    assert os.path.isfile(orig_file)
    original_genome = open(orig_file).read().strip()
    args = {
        'out_dir': './puma_out_test',
        'data_dir': 'data_dir',
        'input_format': 'genbank',
        'min_prot_len': 25,
        'e_value': 1e-05,
        'for_user_dir': './puma_out_test/for_user',
        'program_files_dir': './puma_out_test/program_files'
    }
    new_file = os.path.join('data_dir', 'tests', 'linearized_genome.txt')
    assert os.path.isfile(new_file)
    new_genome = open(new_file).read().strip()
    res = puma.linearize_genome(Seq(original_genome), args)
    assert res == new_genome


# --------------------------------------------------
def test_verify_l1():
    """Docstring"""

    warnings.simplefilter('ignore', BiopythonWarning)

    proteins = 'proteins.csv'
    assert os.path.isfile(proteins)

    found = dict()
    with open(proteins) as fh:
        reader = csv.reader(fh)
        for row in reader:
            name = row[0]
            found[name] = [int(row[1]), int(row[2]), row[3], Seq(row[4])]

    res = puma.verify_l1(found)
    assert res == found


# --------------------------------------------------
def test_blast_main_orfs():
    """Docstring"""

    warnings.simplefilter('ignore', BiopythonWarning)
    new_file = os.path.join('data_dir', 'tests', 'linearized_genome.txt')
    assert os.path.isfile(new_file)
    new_genome = open(new_file).read().strip()
    assert new_genome
    orfs = os.path.join('data_dir', 'tests', 'orfs.fa')
    assert os.path.isfile(orfs)
    blast_sub = os.path.join('data_dir', 'tests', 'main_blast.fa')
    assert os.path.isfile(blast_sub)
    blast_out = os.path.join('data_dir', 'tests', 'blast_results_main.tab')
    assert os.path.isfile(blast_out)
    args = {
        'out_dir': './puma_out_test',
        'data_dir': 'data_dir',
        'input_format': 'genbank',
        'min_prot_len': 25,
        'e_value': 1e-05,
        'for_user_dir': './puma_out_test/for_user',
        'program_files_dir': './puma_out_test/program_files'
    }
    res_1, res_2 = puma.blast_main_orfs(Seq(new_genome), args)

    assert os.path.isfile(res_1)
    assert os.path.isfile(res_1)


# --------------------------------------------------
def test_identify_main_proteins():
    """Docstring"""

    warnings.simplefilter('ignore', BiopythonWarning)
    new_file = os.path.join('data_dir', 'tests', 'linearized_genome.txt')
    assert os.path.isfile(new_file)
    new_genome = open(new_file).read().strip()
    assert new_genome
    args = {
        'out_dir': './puma_out_test',
        'data_dir': 'data_dir',
        'input_format': 'genbank',
        'min_prot_len': 25,
        'e_value': 1e-05,
        'for_user_dir': './puma_out_test/for_user',
        'program_files_dir': './puma_out_test/program_files'
    }
    found = {
        'E7': [
            8468, 8764,
            'atgcatggagatacacctacattgcatgaatatatgttagatttgcaaccagagacaactgatctctactgttatgagcaattaaatgacagctcagaggaggaggatgaaatagatggtccagctggacaagcagaaccggacagagcccattacaatattgtaaccttttgttgcaagtgtgactctacgcttcggttgtgcgtacaaagcacacacgtagacattcgtactttggaagacctgttaatgggcacactaggaattgtgtgccccatctgttctcagaaaccataa',
            Seq('MHGDTPTLHEYMLDLQPETTDLYCYEQLNDSSEEEDEIDGPAGQAEPDRAHYNI...KP*')
        ],
        'L2': [
            4237, 5658,
            'atgcgacacaaacgttctgcaaaacgcacaaaacgtgcatcggctacccaactttataaaacatgcaaacaggcaggtacatgtccacctgacattatacctaaggttgaaggcaaaactattgctgatcaaatattacaatatggaagtatgggtgtattttttggtgggttaggaattggaacagggtcgggtacaggcggacgcactgggtatattccattgggaacaaggcctcccacagctacagatacacttgctcctgtaagaccccctttaacagtagatcctgtgggcccttctgatccttctatagtttctttagtggaagaaactagttttattgatgctggtgcaccaacatctgtaccttccattcccccagatgtatcaggatttagtattactacttcaactgataccacacctgctatattagatattaataatactgttactactgttactacacataataatcccactttcactgacccatctgtattgcagcctccaacacctgcagaaactggagggcattttacactttcatcatccactattagtacacataattatgaagaaattcctatggatacatttattgttagcacaaaccctaacacagtaactagtagcacacccataccagggtctcgcccagtggcacgcctaggattatatagtcgcacaacacaacaggttaaagttgtagaccctgcttttgtaaccactcccactaaacttattacatatgataatcctgcatatgaaggtatagatgtggataatacattatatttttctagtaatgataatagtattaatatagctccagatcctgactttttggatatagttgctttacataggccagcattaacctctaggcgtactggcattaggtacagtagaattggtaataaacaaacactacgtactcgtagtggaaaatctataggtgctaaggtacattattattatgatttaagtactattgatcctgcagaagaaatagaattacaaactataacaccttctacatatactaccacttcacatgcagcctcacctacttctattaataatggattatatgatatttatgcagatgactttattacagatacttctacaaccccggtaccatctgtaccctctacatctttatcaggttatattcctgcaaatacaacaattccttttggtggtgcatacaatattcctttagtatcaggtcctgatatacccattaatataactgaccaagctccttcattaattcctatagttccagggtctccacaatatacaattattgctgatgcaggtgacttttatttacatcctagttattacatgttacgaaaacgacgtaaacgtttaccatattttttttcagatgtctctttggctgcctag',
            Seq('MRHKRSAKRTKRASATQLYKTCKQAGTCPPDIIPKVEGKTIADQILQYGSMGVF...AA*')
        ],
        'E6': [
            7989, 8465,
            'atgcaccaaaagagaactgcaatgtttcaggacccacaggagcgacccagaaagttaccacagttatgcacagagctgcaaacaactatacatgatataatattagaatgtgtgtactgcaagcaacagttactgcgacgtgaggtatatgactttgcttttcgggatttatgcatagtatatagagatgggaatccatatgctgtatgtgataaatgtttaaagttttattctaaaattagtgagtatagacattattgttatagtttgtatggaacaacattagaacagcaatacaacaaaccgttgtgtgatttgttaattaggtgtattaactgtcaaaagccactgtgtcctgaagaaaagcaaagacatctggacaaaaagcaaagattccataatataaggggtcggtggaccggtcgatgtatgtcttgttgcagatcatcaagaacacgtagagaaacccagctgtaa',
            Seq('MHQKRTAMFQDPQERPRKLPQLCTELQTTIHDIILECVYCKQQLLRREVYDFAF...QL*')
        ],
        'E2': [
            2756, 3853,
            'atggagactctttgccaacgtttaaatgtgtgtcaggacaaaatactaacacattatgaaaatgatagtacagacctacgtgaccatatagactattggaaacacatgcgcctagaatgtgctatttattacaaggccagagaaatgggatttaaacatattaaccaccaggtggtgccaacactggctgtatcaaagaataaagcattacaagcaattgaactgcaactaacgttagaaacaatatataactcacaatatagtaatgaaaagtggacattacaagacgttagccttgaagtgtatttaactgcaccaacaggatgtataaaaaaacatggatatacagtggaagtgcagtttgatggagacatatgcaatacaatgcattatacaaactggacacatatatatatttgtgaagaagcatcagtaactgtggtagagggtcaagttgactattatggtttatattatgttcatgaaggaatacgaacatattttgtgcagtttaaagatgatgcagaaaaatatagtaaaaataaagtatgggaagttcatgcgggtggtcaggtaatattatgtcctacatctgtgtttagcagcaacgaagtatcctctcctgaaattattaggcagcacttggccaaccaccccgccgcgacccataccaaagccgtcgccttgggcaccgaagaaacacagacgactatccagcgaccaagatcagagccagacaccggaaacccctgccacaccactaagttgttgcacagagactcagtggacagtgctccaatcctcactgcatttaacagctcacacaaaggacggattaactgtaatagtaacactacacccatagtacatttaaaaggtgatgctaatactttaaaatgtttaagatatagatttaaaaagcattgtacattgtatactgcagtgtcgtctacatggcattggacaggacataatgtaaaacataaaagtgcaattgttacacttacatatgatagtgaatggcaacgtgaccaatttttgtctcaagttaaaataccaaaaactattacagtgtctactggatttatgtctatatga',
            Seq('METLCQRLNVCQDKILTHYENDSTDLRDHIDYWKHMRLECAIYYKAREMGFKHI...SI*')
        ],
        'E1': [
            8771, 9907,
            'atggctgatcctgcaggtaccaatggggaagagggtacgggatgtaatggatggttttatgtagaggctgtagtggaaaaaaaaacaggggatgctatatcagatgacgagaacgaaaatgacagtgatacaggtgaagatttggtagattttatagtaaatgataatgattatttaacacaggcagaaacagagacagcacatgcgttgtttactgcacaggaagcaaaacaacatagagatgcagtacaggttctaaaacgaaagtatttgggtagtccacttagtgatattagtggatgtgtagacaataatattagtcctagattaaaagctatatgtatagaaaaacaaagtagagctgcaaaaaggagattatttgaaagcgaagacagcgggtatggcaatactgaagtggaaactcagcagatgttacaggtagaagggcgccatgagactgaaacaccatgtagtcagtatagtggtggaagtgggggtggttgcagtcagtacagtagtggaagtgggggagagggtgttagtgaaagacacactatatgccaaacaccacttacaaatattttaaatgtactaaaaactagtaatgcaaaggcagcaatgttagcaaaatttaaagagttatacggggtgagtttttcagaattagtaagaccatttaaaagtaataaatcaacgtgttgcgattggtgtattgctgcatttggacttacacccagtatagctgacagtataaaaacactattacaacaatattgtttatatttacacattcaaagtttagcatgttcatggggaatggttgtgttactattagtaagatataaatgtggaaaaaatagagaaacaattgaaaaattgctgtctaaactattatgtgtgtctccaatgtgtatgatgatagagcctccaaaattgcgtagtacagcagcagcattatattggtataaaacaggtatatcaaatattagtgaagtgtatggagacacgccagaatggatacaaagacaaacagtattacaacatagttttaatgattgtacatttgaattatcacagatggtacaatgggcctacgataatgacatagtagacgatagtgaaattgcatataaata',
            Seq('MADPAGTNGEEGTGCNGWFYVEAVVEKKTGDAISDDENENDSDTGEDLVDFIVN...AYK')
        ],
        'L1': [
            5639, 7156,
            'atgtctctttggctgcctagtgaggccactgtctacttgcctcctgtcccagtatctaaggttgtaagcacggatgaatatgttgcacgcacaaacatatattatcatgcaggaacatccagactacttgcagttggacatccctattttcctattaaaaaacctaacaataacaaaatattagttcctaaagtatcaggattacaatacagggtatttagaatacatttacctgaccccaataagtttggttttcctgacacctcattttataatccagatacacagcggctggtttgggcctgtgtaggtgttgaggtaggtcgtggtcagccattaggtgtgggcattagtggccatcctttattaaataaattggatgacacagaaaatgctagtgcttatgcagcaaatgcaggtgtggataatagagaatgtatatctatggattacaaacaaacacaattgtgtttaattggttgcaaaccacctataggggaacactggggcaaaggatccccatgtaccaatgttgcagtaaatccaggtgattgtccaccattagagttaataaacacagttattcaggatggtgatatggttgatactggctttggtgctatggactttactacattacaggctaacaaaagtgaagttccactggatatttgtacatctatttgcaaatatccagattatattaaaatggtgtcagaaccatatggcgacagcttatttttttatttacgaagggaacaaatgtttgttagacatttatttaatagggctggtactgttggtgaaaatgtaccagacgatttatacattaaaggctctgggtctactgcaaatttagccagttcaaattattttcctacacctagtggttctatggttacctctgatgcccaaatattcaataaaccttattggttacaacgagcacagggccacaataatggcatttgttggggtaaccaactatttgttactgttgttgatactacacgcagtacaaatatgtcattatgtgctgccatatctacttcagaaactacatataaaaatactaactttaaggagtacctacgacatggggaggaatatgatttacagtttatttttcaactgtgcaaaataaccttaactgcagacgttatgacatacatacattctatgaattccactattttggaggactggaattttggtctacaacctcccccaggaggcacactagaagatacttataggtttgtaacatcccaggcaattgcttgtcaaaaacatacacctccagcacctaaagaagatccccttaaaaaatacactttttgggaagtaaatttaaaggaaaagttttctgcagacctagatcagtttcctttaggacgcaaatttttactacaagcaggattgaaggccaaaccaaaatttacattaggaaaacgaaaagctacacccaccacctcatctacctctacaactgctaaacgcaaaaaacgtaagctgtaa',
            Seq('MSLWLPSEATVYLPPVPVSKVVSTDEYVARTNIYYHAGTSRLLAVGHPYFPIKK...KL*')
        ]
    }
    res = puma.identify_main_proteins(Seq(new_genome), args)
    assert isinstance(res, dict)
    assert res.keys() == found.keys()


# --------------------------------------------------
def test_blast_e5_variants():
    """Docstring"""

    warnings.simplefilter('ignore', BiopythonWarning)
    virus = {
        'E7': [
            1312, 1608,
            'atgcatggagatacacctacattgcatgaatatatgttagatttgcaaccagagacaactgatctctactgttatgagcaattaaatgacagctcagaggaggaggatgaaatagatggtccagctggacaagcagaaccggacagagcccattacaatattgtaaccttttgttgcaagtgtgactctacgcttcggttgtgcgtacaaagcacacacgtagacattcgtactttggaagacctgttaatgggcacactaggaattgtgtgccccatctgttctcagaaaccataa',
            Seq('MHGDTPTLHEYMLDLQPETTDLYCYEQLNDSSEEEDEIDGPAGQAEPDRAHYNI...KP*')
        ],
        'E1': [
            1615, 3564,
            'atggctgatcctgcaggtaccaatggggaagagggtacgggatgtaatggatggttttatgtagaggctgtagtggaaaaaaaaacaggggatgctatatcagatgacgagaacgaaaatgacagtgatacaggtgaagatttggtagattttatagtaaatgataatgattatttaacacaggcagaaacagagacagcacatgcgttgtttactgcacaggaagcaaaacaacatagagatgcagtacaggttctaaaacgaaagtatttgggtagtccacttagtgatattagtggatgtgtagacaataatattagtcctagattaaaagctatatgtatagaaaaacaaagtagagctgcaaaaaggagattatttgaaagcgaagacagcgggtatggcaatactgaagtggaaactcagcagatgttacaggtagaagggcgccatgagactgaaacaccatgtagtcagtatagtggtggaagtgggggtggttgcagtcagtacagtagtggaagtgggggagagggtgttagtgaaagacacactatatgccaaacaccacttacaaatattttaaatgtactaaaaactagtaatgcaaaggcagcaatgttagcaaaatttaaagagttatacggggtgagtttttcagaattagtaagaccatttaaaagtaataaatcaacgtgttgcgattggtgtattgctgcatttggacttacacccagtatagctgacagtataaaaacactattacaacaatattgtttatatttacacattcaaagtttagcatgttcatggggaatggttgtgttactattagtaagatataaatgtggaaaaaatagagaaacaattgaaaaattgctgtctaaactattatgtgtgtctccaatgtgtatgatgatagagcctccaaaattgcgtagtacagcagcagcattatattggtataaaacaggtatatcaaatattagtgaagtgtatggagacacgccagaatggatacaaagacaaacagtattacaacatagttttaatgattgtacatttgaattatcacagatggtacaatgggcctacgataatgacatagtagacgatagtgaaattgcatataaatatgcacaattggcagacactaatagtaatgcaagtgcctttctaaaaagtaattcacaggcaaaaattgtaaaggattgtgcaacaatgtgtagacattataaacgagcagaaaaaaaacaaatgagtatgagtcaatggataaaatatagatgtgatagggtagatgatggaggtgattggaagcaaattgttatgtttttaaggtatcaaggtgtagagtttatgtcatttttaactgcattaaaaagatttttgcaaggcatacctaaaaaaaattgcatattactatatggtgcagctaacacaggtaaatcattatttggtatgagtttaatgaaatttctgcaagggtctgtaatatgttttgtaaattctaaaagccatttttggttacaaccattagcagatgccaaaataggtatgttagatgatgctacagtgccctgttggaactacatagatgacaatttaagaaatgcattggatggaaatttagtttctatggatgtaaagcatagaccattggtacaactaaaatgccctccattattaattacatctaacattaatgctggtacagattctaggtggccttatttacataatagattggtggtgtttacatttcctaatgagtttccatttgacgaaaacggaaatccagtgtatgagcttaatgataagaactggaaatcctttttctcaaggacgtggtccagattaagtttgcacgaggacgaggacaaggaaaacgatggagactctttgccaacgtttaaatgtgtgtcaggacaaaatactaacacattatga',
            Seq('MADPAGTNGEEGTGCNGWFYVEAVVEKKTGDAISDDENENDSDTGEDLVDFIVN...TL*')
        ],
        'L2': [
            4987, 6408,
            'atgcgacacaaacgttctgcaaaacgcacaaaacgtgcatcggctacccaactttataaaacatgcaaacaggcaggtacatgtccacctgacattatacctaaggttgaaggcaaaactattgctgatcaaatattacaatatggaagtatgggtgtattttttggtgggttaggaattggaacagggtcgggtacaggcggacgcactgggtatattccattgggaacaaggcctcccacagctacagatacacttgctcctgtaagaccccctttaacagtagatcctgtgggcccttctgatccttctatagtttctttagtggaagaaactagttttattgatgctggtgcaccaacatctgtaccttccattcccccagatgtatcaggatttagtattactacttcaactgataccacacctgctatattagatattaataatactgttactactgttactacacataataatcccactttcactgacccatctgtattgcagcctccaacacctgcagaaactggagggcattttacactttcatcatccactattagtacacataattatgaagaaattcctatggatacatttattgttagcacaaaccctaacacagtaactagtagcacacccataccagggtctcgcccagtggcacgcctaggattatatagtcgcacaacacaacaggttaaagttgtagaccctgcttttgtaaccactcccactaaacttattacatatgataatcctgcatatgaaggtatagatgtggataatacattatatttttctagtaatgataatagtattaatatagctccagatcctgactttttggatatagttgctttacataggccagcattaacctctaggcgtactggcattaggtacagtagaattggtaataaacaaacactacgtactcgtagtggaaaatctataggtgctaaggtacattattattatgatttaagtactattgatcctgcagaagaaatagaattacaaactataacaccttctacatatactaccacttcacatgcagcctcacctacttctattaataatggattatatgatatttatgcagatgactttattacagatacttctacaaccccggtaccatctgtaccctctacatctttatcaggttatattcctgcaaatacaacaattccttttggtggtgcatacaatattcctttagtatcaggtcctgatatacccattaatataactgaccaagctccttcattaattcctatagttccagggtctccacaatatacaattattgctgatgcaggtgacttttatttacatcctagttattacatgttacgaaaacgacgtaaacgtttaccatattttttttcagatgtctctttggctgcctag',
            Seq('MRHKRSAKRTKRASATQLYKTCKQAGTCPPDIIPKVEGKTIADQILQYGSMGVF...AA*')
        ],
        'E6': [
            833, 1309,
            'atgcaccaaaagagaactgcaatgtttcaggacccacaggagcgacccagaaagttaccacagttatgcacagagctgcaaacaactatacatgatataatattagaatgtgtgtactgcaagcaacagttactgcgacgtgaggtatatgactttgcttttcgggatttatgcatagtatatagagatgggaatccatatgctgtatgtgataaatgtttaaagttttattctaaaattagtgagtatagacattattgttatagtttgtatggaacaacattagaacagcaatacaacaaaccgttgtgtgatttgttaattaggtgtattaactgtcaaaagccactgtgtcctgaagaaaagcaaagacatctggacaaaaagcaaagattccataatataaggggtcggtggaccggtcgatgtatgtcttgttgcagatcatcaagaacacgtagagaaacccagctgtaa',
            Seq('MHQKRTAMFQDPQERPRKLPQLCTELQTTIHDIILECVYCKQQLLRREVYDFAF...QL*')
        ],
        'E2': [
            3506, 4603,
            'atggagactctttgccaacgtttaaatgtgtgtcaggacaaaatactaacacattatgaaaatgatagtacagacctacgtgaccatatagactattggaaacacatgcgcctagaatgtgctatttattacaaggccagagaaatgggatttaaacatattaaccaccaggtggtgccaacactggctgtatcaaagaataaagcattacaagcaattgaactgcaactaacgttagaaacaatatataactcacaatatagtaatgaaaagtggacattacaagacgttagccttgaagtgtatttaactgcaccaacaggatgtataaaaaaacatggatatacagtggaagtgcagtttgatggagacatatgcaatacaatgcattatacaaactggacacatatatatatttgtgaagaagcatcagtaactgtggtagagggtcaagttgactattatggtttatattatgttcatgaaggaatacgaacatattttgtgcagtttaaagatgatgcagaaaaatatagtaaaaataaagtatgggaagttcatgcgggtggtcaggtaatattatgtcctacatctgtgtttagcagcaacgaagtatcctctcctgaaattattaggcagcacttggccaaccaccccgccgcgacccataccaaagccgtcgccttgggcaccgaagaaacacagacgactatccagcgaccaagatcagagccagacaccggaaacccctgccacaccactaagttgttgcacagagactcagtggacagtgctccaatcctcactgcatttaacagctcacacaaaggacggattaactgtaatagtaacactacacccatagtacatttaaaaggtgatgctaatactttaaaatgtttaagatatagatttaaaaagcattgtacattgtatactgcagtgtcgtctacatggcattggacaggacataatgtaaaacataaaagtgcaattgttacacttacatatgatagtgaatggcaacgtgaccaatttttgtctcaagttaaaataccaaaaactattacagtgtctactggatttatgtctatatga',
            Seq('METLCQRLNVCQDKILTHYENDSTDLRDHIDYWKHMRLECAIYYKAREMGFKHI...SI*')
        ],
        'L1': [
            6389, 7906,
            'atgtctctttggctgcctagtgaggccactgtctacttgcctcctgtcccagtatctaaggttgtaagcacggatgaatatgttgcacgcacaaacatatattatcatgcaggaacatccagactacttgcagttggacatccctattttcctattaaaaaacctaacaataacaaaatattagttcctaaagtatcaggattacaatacagggtatttagaatacatttacctgaccccaataagtttggttttcctgacacctcattttataatccagatacacagcggctggtttgggcctgtgtaggtgttgaggtaggtcgtggtcagccattaggtgtgggcattagtggccatcctttattaaataaattggatgacacagaaaatgctagtgcttatgcagcaaatgcaggtgtggataatagagaatgtatatctatggattacaaacaaacacaattgtgtttaattggttgcaaaccacctataggggaacactggggcaaaggatccccatgtaccaatgttgcagtaaatccaggtgattgtccaccattagagttaataaacacagttattcaggatggtgatatggttgatactggctttggtgctatggactttactacattacaggctaacaaaagtgaagttccactggatatttgtacatctatttgcaaatatccagattatattaaaatggtgtcagaaccatatggcgacagcttatttttttatttacgaagggaacaaatgtttgttagacatttatttaatagggctggtactgttggtgaaaatgtaccagacgatttatacattaaaggctctgggtctactgcaaatttagccagttcaaattattttcctacacctagtggttctatggttacctctgatgcccaaatattcaataaaccttattggttacaacgagcacagggccacaataatggcatttgttggggtaaccaactatttgttactgttgttgatactacacgcagtacaaatatgtcattatgtgctgccatatctacttcagaaactacatataaaaatactaactttaaggagtacctacgacatggggaggaatatgatttacagtttatttttcaactgtgcaaaataaccttaactgcagacgttatgacatacatacattctatgaattccactattttggaggactggaattttggtctacaacctcccccaggaggcacactagaagatacttataggtttgtaacatcccaggcaattgcttgtcaaaaacatacacctccagcacctaaagaagatccccttaaaaaatacactttttgggaagtaaatttaaaggaaaagttttctgcagacctagatcagtttcctttaggacgcaaatttttactacaagcaggattgaaggccaaaccaaaatttacattaggaaaacgaaaagctacacccaccacctcatctacctctacaactgctaaacgcaaaaaacgtaagctgtaa',
            Seq('MSLWLPSEATVYLPPVPVSKVVSTDEYVARTNIYYHAGTSRLLAVGHPYFPIKK...KL*')
        ],
        'name':
        'Human papillomavirus 16 (HPV16)',
        'accession':
        'HPV16REF',
        'genome':
        'gtattgtatgtatgttgaattagtgttgtttgttgtgtatatgtttgtatgtgcttgtatgtgcttgtaaatattaagttgtatgtgtgtttgtatgtatggtataataaacacgtgtgtatgtgtttttaaatgcttgtgtaactattgtgtcatgcaacataaataaacttattgtttcaacacctactaattgtgttgtggttattcattgtatataaactatatttgctacatcctgtttttgttttatatatactatattttgtagcgccagcggccattttgtagcttcaaccgaattcggttgcatgctttttggcacaaaatgtgtttttttaaatagttctatgtcagcaactatggtttaaacttgtacgtttcctgcttgccatgcgtgccaaatccctgttttcctgacctgcactgcttgccaaccattccattgttttttacactgcactatgtgcaactactgaatcactatgtacattgtgtcatataaaataaatcactatgcgccaacgccttacataccgctgttaggcacatatttttggcttgttttaactaacctaattgcatatttggcataaggtttaaacttctaaggccaactaaatgtcaccctagttcatacatgaactgtgtaaaggttagtcatacattgttcatttgtaaaactgcacatgggtgtgtgcaaaccgttttgggttacacatttacaagcaacttatataataatactaaactacaataattcatgtataaaactaagggcgtaaccgaaatcggttgaaccgaaaccggttagtataaaagcagacattttatgcaccaaaagagaactgcaatgtttcaggacccacaggagcgacccagaaagttaccacagttatgcacagagctgcaaacaactatacatgatataatattagaatgtgtgtactgcaagcaacagttactgcgacgtgaggtatatgactttgcttttcgggatttatgcatagtatatagagatgggaatccatatgctgtatgtgataaatgtttaaagttttattctaaaattagtgagtatagacattattgttatagtttgtatggaacaacattagaacagcaatacaacaaaccgttgtgtgatttgttaattaggtgtattaactgtcaaaagccactgtgtcctgaagaaaagcaaagacatctggacaaaaagcaaagattccataatataaggggtcggtggaccggtcgatgtatgtcttgttgcagatcatcaagaacacgtagagaaacccagctgtaatcatgcatggagatacacctacattgcatgaatatatgttagatttgcaaccagagacaactgatctctactgttatgagcaattaaatgacagctcagaggaggaggatgaaatagatggtccagctggacaagcagaaccggacagagcccattacaatattgtaaccttttgttgcaagtgtgactctacgcttcggttgtgcgtacaaagcacacacgtagacattcgtactttggaagacctgttaatgggcacactaggaattgtgtgccccatctgttctcagaaaccataatctaccatggctgatcctgcaggtaccaatggggaagagggtacgggatgtaatggatggttttatgtagaggctgtagtggaaaaaaaaacaggggatgctatatcagatgacgagaacgaaaatgacagtgatacaggtgaagatttggtagattttatagtaaatgataatgattatttaacacaggcagaaacagagacagcacatgcgttgtttactgcacaggaagcaaaacaacatagagatgcagtacaggttctaaaacgaaagtatttgggtagtccacttagtgatattagtggatgtgtagacaataatattagtcctagattaaaagctatatgtatagaaaaacaaagtagagctgcaaaaaggagattatttgaaagcgaagacagcgggtatggcaatactgaagtggaaactcagcagatgttacaggtagaagggcgccatgagactgaaacaccatgtagtcagtatagtggtggaagtgggggtggttgcagtcagtacagtagtggaagtgggggagagggtgttagtgaaagacacactatatgccaaacaccacttacaaatattttaaatgtactaaaaactagtaatgcaaaggcagcaatgttagcaaaatttaaagagttatacggggtgagtttttcagaattagtaagaccatttaaaagtaataaatcaacgtgttgcgattggtgtattgctgcatttggacttacacccagtatagctgacagtataaaaacactattacaacaatattgtttatatttacacattcaaagtttagcatgttcatggggaatggttgtgttactattagtaagatataaatgtggaaaaaatagagaaacaattgaaaaattgctgtctaaactattatgtgtgtctccaatgtgtatgatgatagagcctccaaaattgcgtagtacagcagcagcattatattggtataaaacaggtatatcaaatattagtgaagtgtatggagacacgccagaatggatacaaagacaaacagtattacaacatagttttaatgattgtacatttgaattatcacagatggtacaatgggcctacgataatgacatagtagacgatagtgaaattgcatataaatatgcacaattggcagacactaatagtaatgcaagtgcctttctaaaaagtaattcacaggcaaaaattgtaaaggattgtgcaacaatgtgtagacattataaacgagcagaaaaaaaacaaatgagtatgagtcaatggataaaatatagatgtgatagggtagatgatggaggtgattggaagcaaattgttatgtttttaaggtatcaaggtgtagagtttatgtcatttttaactgcattaaaaagatttttgcaaggcatacctaaaaaaaattgcatattactatatggtgcagctaacacaggtaaatcattatttggtatgagtttaatgaaatttctgcaagggtctgtaatatgttttgtaaattctaaaagccatttttggttacaaccattagcagatgccaaaataggtatgttagatgatgctacagtgccctgttggaactacatagatgacaatttaagaaatgcattggatggaaatttagtttctatggatgtaaagcatagaccattggtacaactaaaatgccctccattattaattacatctaacattaatgctggtacagattctaggtggccttatttacataatagattggtggtgtttacatttcctaatgagtttccatttgacgaaaacggaaatccagtgtatgagcttaatgataagaactggaaatcctttttctcaaggacgtggtccagattaagtttgcacgaggacgaggacaaggaaaacgatggagactctttgccaacgtttaaatgtgtgtcaggacaaaatactaacacattatgaaaatgatagtacagacctacgtgaccatatagactattggaaacacatgcgcctagaatgtgctatttattacaaggccagagaaatgggatttaaacatattaaccaccaggtggtgccaacactggctgtatcaaagaataaagcattacaagcaattgaactgcaactaacgttagaaacaatatataactcacaatatagtaatgaaaagtggacattacaagacgttagccttgaagtgtatttaactgcaccaacaggatgtataaaaaaacatggatatacagtggaagtgcagtttgatggagacatatgcaatacaatgcattatacaaactggacacatatatatatttgtgaagaagcatcagtaactgtggtagagggtcaagttgactattatggtttatattatgttcatgaaggaatacgaacatattttgtgcagtttaaagatgatgcagaaaaatatagtaaaaataaagtatgggaagttcatgcgggtggtcaggtaatattatgtcctacatctgtgtttagcagcaacgaagtatcctctcctgaaattattaggcagcacttggccaaccaccccgccgcgacccataccaaagccgtcgccttgggcaccgaagaaacacagacgactatccagcgaccaagatcagagccagacaccggaaacccctgccacaccactaagttgttgcacagagactcagtggacagtgctccaatcctcactgcatttaacagctcacacaaaggacggattaactgtaatagtaacactacacccatagtacatttaaaaggtgatgctaatactttaaaatgtttaagatatagatttaaaaagcattgtacattgtatactgcagtgtcgtctacatggcattggacaggacataatgtaaaacataaaagtgcaattgttacacttacatatgatagtgaatggcaacgtgaccaatttttgtctcaagttaaaataccaaaaactattacagtgtctactggatttatgtctatatgacaaatcttgatactgcatccacaacattactggcgtgctttttgctttgcttttgtgtgcttttgtgtgtctgcctattaatacgtccgctgcttttgtctgtgtctacatacacatcattaataatattggtattactattgtggataacagcagcctctgcgtttaggtgttttattgtatatattatatttgtttatataccattatttttaatacatacacatgcacgctttttaattacataatgtatatgtacataatgtaattgttacatataattgttgtataccataacttactattttttcttttttattttcatatataatttttttttttgtttgtttgtttgttttttaataaactgttattacttaacaatgcgacacaaacgttctgcaaaacgcacaaaacgtgcatcggctacccaactttataaaacatgcaaacaggcaggtacatgtccacctgacattatacctaaggttgaaggcaaaactattgctgatcaaatattacaatatggaagtatgggtgtattttttggtgggttaggaattggaacagggtcgggtacaggcggacgcactgggtatattccattgggaacaaggcctcccacagctacagatacacttgctcctgtaagaccccctttaacagtagatcctgtgggcccttctgatccttctatagtttctttagtggaagaaactagttttattgatgctggtgcaccaacatctgtaccttccattcccccagatgtatcaggatttagtattactacttcaactgataccacacctgctatattagatattaataatactgttactactgttactacacataataatcccactttcactgacccatctgtattgcagcctccaacacctgcagaaactggagggcattttacactttcatcatccactattagtacacataattatgaagaaattcctatggatacatttattgttagcacaaaccctaacacagtaactagtagcacacccataccagggtctcgcccagtggcacgcctaggattatatagtcgcacaacacaacaggttaaagttgtagaccctgcttttgtaaccactcccactaaacttattacatatgataatcctgcatatgaaggtatagatgtggataatacattatatttttctagtaatgataatagtattaatatagctccagatcctgactttttggatatagttgctttacataggccagcattaacctctaggcgtactggcattaggtacagtagaattggtaataaacaaacactacgtactcgtagtggaaaatctataggtgctaaggtacattattattatgatttaagtactattgatcctgcagaagaaatagaattacaaactataacaccttctacatatactaccacttcacatgcagcctcacctacttctattaataatggattatatgatatttatgcagatgactttattacagatacttctacaaccccggtaccatctgtaccctctacatctttatcaggttatattcctgcaaatacaacaattccttttggtggtgcatacaatattcctttagtatcaggtcctgatatacccattaatataactgaccaagctccttcattaattcctatagttccagggtctccacaatatacaattattgctgatgcaggtgacttttatttacatcctagttattacatgttacgaaaacgacgtaaacgtttaccatattttttttcagatgtctctttggctgcctagtgaggccactgtctacttgcctcctgtcccagtatctaaggttgtaagcacggatgaatatgttgcacgcacaaacatatattatcatgcaggaacatccagactacttgcagttggacatccctattttcctattaaaaaacctaacaataacaaaatattagttcctaaagtatcaggattacaatacagggtatttagaatacatttacctgaccccaataagtttggttttcctgacacctcattttataatccagatacacagcggctggtttgggcctgtgtaggtgttgaggtaggtcgtggtcagccattaggtgtgggcattagtggccatcctttattaaataaattggatgacacagaaaatgctagtgcttatgcagcaaatgcaggtgtggataatagagaatgtatatctatggattacaaacaaacacaattgtgtttaattggttgcaaaccacctataggggaacactggggcaaaggatccccatgtaccaatgttgcagtaaatccaggtgattgtccaccattagagttaataaacacagttattcaggatggtgatatggttgatactggctttggtgctatggactttactacattacaggctaacaaaagtgaagttccactggatatttgtacatctatttgcaaatatccagattatattaaaatggtgtcagaaccatatggcgacagcttatttttttatttacgaagggaacaaatgtttgttagacatttatttaatagggctggtactgttggtgaaaatgtaccagacgatttatacattaaaggctctgggtctactgcaaatttagccagttcaaattattttcctacacctagtggttctatggttacctctgatgcccaaatattcaataaaccttattggttacaacgagcacagggccacaataatggcatttgttggggtaaccaactatttgttactgttgttgatactacacgcagtacaaatatgtcattatgtgctgccatatctacttcagaaactacatataaaaatactaactttaaggagtacctacgacatggggaggaatatgatttacagtttatttttcaactgtgcaaaataaccttaactgcagacgttatgacatacatacattctatgaattccactattttggaggactggaattttggtctacaacctcccccaggaggcacactagaagatacttataggtttgtaacatcccaggcaattgcttgtcaaaaacatacacctccagcacctaaagaagatccccttaaaaaatacactttttgggaagtaaatttaaaggaaaagttttctgcagacctagatcagtttcctttaggacgcaaatttttactacaagcaggattgaaggccaaaccaaaatttacattaggaaaacgaaaagctacacccaccacctcatctacctctacaactgctaaacgcaaaaaacgtaagctgtaa'
    }
    args = {
        'out_dir': './puma_out_test',
        'data_dir': 'data_dir',
        'input_format': 'genbank',
        'min_prot_len': 25,
        'e_value': 1e-05,
        'for_user_dir': './puma_out_test/for_user',
        'program_files_dir': './puma_out_test/program_files'
    }
    e5_file = os.path.join('data_dir', 'tests', 'e5_seq.txt')
    assert os.path.isfile(e5_file)
    e5_seq = str(open(e5_file).read().strip())
    orfs = os.path.join('data_dir', 'tests', 'orfs_E5.fa')
    assert os.path.isfile(orfs)
    blast_out = os.path.join('data_dir', 'tests', 'blast_results_E5.tab')
    assert os.path.isfile(blast_out)
    res_1, res_2, res_3 = puma.blast_e5_variants(virus, args)
    assert os.path.isfile(res_1)
    assert os.path.isfile(res_1)
    assert res_3 == e5_seq


# --------------------------------------------------
def test_identify_e5_variants():
    """Docstring"""

    warnings.simplefilter('ignore', BiopythonWarning)
    virus = {
        'E7': [
            1312, 1608,
            'atgcatggagatacacctacattgcatgaatatatgttagatttgcaaccagagacaactgatctctactgttatgagcaattaaatgacagctcagaggaggaggatgaaatagatggtccagctggacaagcagaaccggacagagcccattacaatattgtaaccttttgttgcaagtgtgactctacgcttcggttgtgcgtacaaagcacacacgtagacattcgtactttggaagacctgttaatgggcacactaggaattgtgtgccccatctgttctcagaaaccataa',
            Seq('MHGDTPTLHEYMLDLQPETTDLYCYEQLNDSSEEEDEIDGPAGQAEPDRAHYNI...KP*')
        ],
        'E1': [
            1615, 3564,
            'atggctgatcctgcaggtaccaatggggaagagggtacgggatgtaatggatggttttatgtagaggctgtagtggaaaaaaaaacaggggatgctatatcagatgacgagaacgaaaatgacagtgatacaggtgaagatttggtagattttatagtaaatgataatgattatttaacacaggcagaaacagagacagcacatgcgttgtttactgcacaggaagcaaaacaacatagagatgcagtacaggttctaaaacgaaagtatttgggtagtccacttagtgatattagtggatgtgtagacaataatattagtcctagattaaaagctatatgtatagaaaaacaaagtagagctgcaaaaaggagattatttgaaagcgaagacagcgggtatggcaatactgaagtggaaactcagcagatgttacaggtagaagggcgccatgagactgaaacaccatgtagtcagtatagtggtggaagtgggggtggttgcagtcagtacagtagtggaagtgggggagagggtgttagtgaaagacacactatatgccaaacaccacttacaaatattttaaatgtactaaaaactagtaatgcaaaggcagcaatgttagcaaaatttaaagagttatacggggtgagtttttcagaattagtaagaccatttaaaagtaataaatcaacgtgttgcgattggtgtattgctgcatttggacttacacccagtatagctgacagtataaaaacactattacaacaatattgtttatatttacacattcaaagtttagcatgttcatggggaatggttgtgttactattagtaagatataaatgtggaaaaaatagagaaacaattgaaaaattgctgtctaaactattatgtgtgtctccaatgtgtatgatgatagagcctccaaaattgcgtagtacagcagcagcattatattggtataaaacaggtatatcaaatattagtgaagtgtatggagacacgccagaatggatacaaagacaaacagtattacaacatagttttaatgattgtacatttgaattatcacagatggtacaatgggcctacgataatgacatagtagacgatagtgaaattgcatataaatatgcacaattggcagacactaatagtaatgcaagtgcctttctaaaaagtaattcacaggcaaaaattgtaaaggattgtgcaacaatgtgtagacattataaacgagcagaaaaaaaacaaatgagtatgagtcaatggataaaatatagatgtgatagggtagatgatggaggtgattggaagcaaattgttatgtttttaaggtatcaaggtgtagagtttatgtcatttttaactgcattaaaaagatttttgcaaggcatacctaaaaaaaattgcatattactatatggtgcagctaacacaggtaaatcattatttggtatgagtttaatgaaatttctgcaagggtctgtaatatgttttgtaaattctaaaagccatttttggttacaaccattagcagatgccaaaataggtatgttagatgatgctacagtgccctgttggaactacatagatgacaatttaagaaatgcattggatggaaatttagtttctatggatgtaaagcatagaccattggtacaactaaaatgccctccattattaattacatctaacattaatgctggtacagattctaggtggccttatttacataatagattggtggtgtttacatttcctaatgagtttccatttgacgaaaacggaaatccagtgtatgagcttaatgataagaactggaaatcctttttctcaaggacgtggtccagattaagtttgcacgaggacgaggacaaggaaaacgatggagactctttgccaacgtttaaatgtgtgtcaggacaaaatactaacacattatga',
            Seq('MADPAGTNGEEGTGCNGWFYVEAVVEKKTGDAISDDENENDSDTGEDLVDFIVN...TL*')
        ],
        'L2': [
            4987, 6408,
            'atgcgacacaaacgttctgcaaaacgcacaaaacgtgcatcggctacccaactttataaaacatgcaaacaggcaggtacatgtccacctgacattatacctaaggttgaaggcaaaactattgctgatcaaatattacaatatggaagtatgggtgtattttttggtgggttaggaattggaacagggtcgggtacaggcggacgcactgggtatattccattgggaacaaggcctcccacagctacagatacacttgctcctgtaagaccccctttaacagtagatcctgtgggcccttctgatccttctatagtttctttagtggaagaaactagttttattgatgctggtgcaccaacatctgtaccttccattcccccagatgtatcaggatttagtattactacttcaactgataccacacctgctatattagatattaataatactgttactactgttactacacataataatcccactttcactgacccatctgtattgcagcctccaacacctgcagaaactggagggcattttacactttcatcatccactattagtacacataattatgaagaaattcctatggatacatttattgttagcacaaaccctaacacagtaactagtagcacacccataccagggtctcgcccagtggcacgcctaggattatatagtcgcacaacacaacaggttaaagttgtagaccctgcttttgtaaccactcccactaaacttattacatatgataatcctgcatatgaaggtatagatgtggataatacattatatttttctagtaatgataatagtattaatatagctccagatcctgactttttggatatagttgctttacataggccagcattaacctctaggcgtactggcattaggtacagtagaattggtaataaacaaacactacgtactcgtagtggaaaatctataggtgctaaggtacattattattatgatttaagtactattgatcctgcagaagaaatagaattacaaactataacaccttctacatatactaccacttcacatgcagcctcacctacttctattaataatggattatatgatatttatgcagatgactttattacagatacttctacaaccccggtaccatctgtaccctctacatctttatcaggttatattcctgcaaatacaacaattccttttggtggtgcatacaatattcctttagtatcaggtcctgatatacccattaatataactgaccaagctccttcattaattcctatagttccagggtctccacaatatacaattattgctgatgcaggtgacttttatttacatcctagttattacatgttacgaaaacgacgtaaacgtttaccatattttttttcagatgtctctttggctgcctag',
            Seq('MRHKRSAKRTKRASATQLYKTCKQAGTCPPDIIPKVEGKTIADQILQYGSMGVF...AA*')
        ],
        'E6': [
            833, 1309,
            'atgcaccaaaagagaactgcaatgtttcaggacccacaggagcgacccagaaagttaccacagttatgcacagagctgcaaacaactatacatgatataatattagaatgtgtgtactgcaagcaacagttactgcgacgtgaggtatatgactttgcttttcgggatttatgcatagtatatagagatgggaatccatatgctgtatgtgataaatgtttaaagttttattctaaaattagtgagtatagacattattgttatagtttgtatggaacaacattagaacagcaatacaacaaaccgttgtgtgatttgttaattaggtgtattaactgtcaaaagccactgtgtcctgaagaaaagcaaagacatctggacaaaaagcaaagattccataatataaggggtcggtggaccggtcgatgtatgtcttgttgcagatcatcaagaacacgtagagaaacccagctgtaa',
            Seq('MHQKRTAMFQDPQERPRKLPQLCTELQTTIHDIILECVYCKQQLLRREVYDFAF...QL*')
        ],
        'E2': [
            3506, 4603,
            'atggagactctttgccaacgtttaaatgtgtgtcaggacaaaatactaacacattatgaaaatgatagtacagacctacgtgaccatatagactattggaaacacatgcgcctagaatgtgctatttattacaaggccagagaaatgggatttaaacatattaaccaccaggtggtgccaacactggctgtatcaaagaataaagcattacaagcaattgaactgcaactaacgttagaaacaatatataactcacaatatagtaatgaaaagtggacattacaagacgttagccttgaagtgtatttaactgcaccaacaggatgtataaaaaaacatggatatacagtggaagtgcagtttgatggagacatatgcaatacaatgcattatacaaactggacacatatatatatttgtgaagaagcatcagtaactgtggtagagggtcaagttgactattatggtttatattatgttcatgaaggaatacgaacatattttgtgcagtttaaagatgatgcagaaaaatatagtaaaaataaagtatgggaagttcatgcgggtggtcaggtaatattatgtcctacatctgtgtttagcagcaacgaagtatcctctcctgaaattattaggcagcacttggccaaccaccccgccgcgacccataccaaagccgtcgccttgggcaccgaagaaacacagacgactatccagcgaccaagatcagagccagacaccggaaacccctgccacaccactaagttgttgcacagagactcagtggacagtgctccaatcctcactgcatttaacagctcacacaaaggacggattaactgtaatagtaacactacacccatagtacatttaaaaggtgatgctaatactttaaaatgtttaagatatagatttaaaaagcattgtacattgtatactgcagtgtcgtctacatggcattggacaggacataatgtaaaacataaaagtgcaattgttacacttacatatgatagtgaatggcaacgtgaccaatttttgtctcaagttaaaataccaaaaactattacagtgtctactggatttatgtctatatga',
            Seq('METLCQRLNVCQDKILTHYENDSTDLRDHIDYWKHMRLECAIYYKAREMGFKHI...SI*')
        ],
        'L1': [
            6389, 7906,
            'atgtctctttggctgcctagtgaggccactgtctacttgcctcctgtcccagtatctaaggttgtaagcacggatgaatatgttgcacgcacaaacatatattatcatgcaggaacatccagactacttgcagttggacatccctattttcctattaaaaaacctaacaataacaaaatattagttcctaaagtatcaggattacaatacagggtatttagaatacatttacctgaccccaataagtttggttttcctgacacctcattttataatccagatacacagcggctggtttgggcctgtgtaggtgttgaggtaggtcgtggtcagccattaggtgtgggcattagtggccatcctttattaaataaattggatgacacagaaaatgctagtgcttatgcagcaaatgcaggtgtggataatagagaatgtatatctatggattacaaacaaacacaattgtgtttaattggttgcaaaccacctataggggaacactggggcaaaggatccccatgtaccaatgttgcagtaaatccaggtgattgtccaccattagagttaataaacacagttattcaggatggtgatatggttgatactggctttggtgctatggactttactacattacaggctaacaaaagtgaagttccactggatatttgtacatctatttgcaaatatccagattatattaaaatggtgtcagaaccatatggcgacagcttatttttttatttacgaagggaacaaatgtttgttagacatttatttaatagggctggtactgttggtgaaaatgtaccagacgatttatacattaaaggctctgggtctactgcaaatttagccagttcaaattattttcctacacctagtggttctatggttacctctgatgcccaaatattcaataaaccttattggttacaacgagcacagggccacaataatggcatttgttggggtaaccaactatttgttactgttgttgatactacacgcagtacaaatatgtcattatgtgctgccatatctacttcagaaactacatataaaaatactaactttaaggagtacctacgacatggggaggaatatgatttacagtttatttttcaactgtgcaaaataaccttaactgcagacgttatgacatacatacattctatgaattccactattttggaggactggaattttggtctacaacctcccccaggaggcacactagaagatacttataggtttgtaacatcccaggcaattgcttgtcaaaaacatacacctccagcacctaaagaagatccccttaaaaaatacactttttgggaagtaaatttaaaggaaaagttttctgcagacctagatcagtttcctttaggacgcaaatttttactacaagcaggattgaaggccaaaccaaaatttacattaggaaaacgaaaagctacacccaccacctcatctacctctacaactgctaaacgcaaaaaacgtaagctgtaa',
            Seq('MSLWLPSEATVYLPPVPVSKVVSTDEYVARTNIYYHAGTSRLLAVGHPYFPIKK...KL*')
        ],
        'name':
        'Human papillomavirus 16 (HPV16)',
        'accession':
        'HPV16REF',
        'genome':
        'gtattgtatgtatgttgaattagtgttgtttgttgtgtatatgtttgtatgtgcttgtatgtgcttgtaaatattaagttgtatgtgtgtttgtatgtatggtataataaacacgtgtgtatgtgtttttaaatgcttgtgtaactattgtgtcatgcaacataaataaacttattgtttcaacacctactaattgtgttgtggttattcattgtatataaactatatttgctacatcctgtttttgttttatatatactatattttgtagcgccagcggccattttgtagcttcaaccgaattcggttgcatgctttttggcacaaaatgtgtttttttaaatagttctatgtcagcaactatggtttaaacttgtacgtttcctgcttgccatgcgtgccaaatccctgttttcctgacctgcactgcttgccaaccattccattgttttttacactgcactatgtgcaactactgaatcactatgtacattgtgtcatataaaataaatcactatgcgccaacgccttacataccgctgttaggcacatatttttggcttgttttaactaacctaattgcatatttggcataaggtttaaacttctaaggccaactaaatgtcaccctagttcatacatgaactgtgtaaaggttagtcatacattgttcatttgtaaaactgcacatgggtgtgtgcaaaccgttttgggttacacatttacaagcaacttatataataatactaaactacaataattcatgtataaaactaagggcgtaaccgaaatcggttgaaccgaaaccggttagtataaaagcagacattttatgcaccaaaagagaactgcaatgtttcaggacccacaggagcgacccagaaagttaccacagttatgcacagagctgcaaacaactatacatgatataatattagaatgtgtgtactgcaagcaacagttactgcgacgtgaggtatatgactttgcttttcgggatttatgcatagtatatagagatgggaatccatatgctgtatgtgataaatgtttaaagttttattctaaaattagtgagtatagacattattgttatagtttgtatggaacaacattagaacagcaatacaacaaaccgttgtgtgatttgttaattaggtgtattaactgtcaaaagccactgtgtcctgaagaaaagcaaagacatctggacaaaaagcaaagattccataatataaggggtcggtggaccggtcgatgtatgtcttgttgcagatcatcaagaacacgtagagaaacccagctgtaatcatgcatggagatacacctacattgcatgaatatatgttagatttgcaaccagagacaactgatctctactgttatgagcaattaaatgacagctcagaggaggaggatgaaatagatggtccagctggacaagcagaaccggacagagcccattacaatattgtaaccttttgttgcaagtgtgactctacgcttcggttgtgcgtacaaagcacacacgtagacattcgtactttggaagacctgttaatgggcacactaggaattgtgtgccccatctgttctcagaaaccataatctaccatggctgatcctgcaggtaccaatggggaagagggtacgggatgtaatggatggttttatgtagaggctgtagtggaaaaaaaaacaggggatgctatatcagatgacgagaacgaaaatgacagtgatacaggtgaagatttggtagattttatagtaaatgataatgattatttaacacaggcagaaacagagacagcacatgcgttgtttactgcacaggaagcaaaacaacatagagatgcagtacaggttctaaaacgaaagtatttgggtagtccacttagtgatattagtggatgtgtagacaataatattagtcctagattaaaagctatatgtatagaaaaacaaagtagagctgcaaaaaggagattatttgaaagcgaagacagcgggtatggcaatactgaagtggaaactcagcagatgttacaggtagaagggcgccatgagactgaaacaccatgtagtcagtatagtggtggaagtgggggtggttgcagtcagtacagtagtggaagtgggggagagggtgttagtgaaagacacactatatgccaaacaccacttacaaatattttaaatgtactaaaaactagtaatgcaaaggcagcaatgttagcaaaatttaaagagttatacggggtgagtttttcagaattagtaagaccatttaaaagtaataaatcaacgtgttgcgattggtgtattgctgcatttggacttacacccagtatagctgacagtataaaaacactattacaacaatattgtttatatttacacattcaaagtttagcatgttcatggggaatggttgtgttactattagtaagatataaatgtggaaaaaatagagaaacaattgaaaaattgctgtctaaactattatgtgtgtctccaatgtgtatgatgatagagcctccaaaattgcgtagtacagcagcagcattatattggtataaaacaggtatatcaaatattagtgaagtgtatggagacacgccagaatggatacaaagacaaacagtattacaacatagttttaatgattgtacatttgaattatcacagatggtacaatgggcctacgataatgacatagtagacgatagtgaaattgcatataaatatgcacaattggcagacactaatagtaatgcaagtgcctttctaaaaagtaattcacaggcaaaaattgtaaaggattgtgcaacaatgtgtagacattataaacgagcagaaaaaaaacaaatgagtatgagtcaatggataaaatatagatgtgatagggtagatgatggaggtgattggaagcaaattgttatgtttttaaggtatcaaggtgtagagtttatgtcatttttaactgcattaaaaagatttttgcaaggcatacctaaaaaaaattgcatattactatatggtgcagctaacacaggtaaatcattatttggtatgagtttaatgaaatttctgcaagggtctgtaatatgttttgtaaattctaaaagccatttttggttacaaccattagcagatgccaaaataggtatgttagatgatgctacagtgccctgttggaactacatagatgacaatttaagaaatgcattggatggaaatttagtttctatggatgtaaagcatagaccattggtacaactaaaatgccctccattattaattacatctaacattaatgctggtacagattctaggtggccttatttacataatagattggtggtgtttacatttcctaatgagtttccatttgacgaaaacggaaatccagtgtatgagcttaatgataagaactggaaatcctttttctcaaggacgtggtccagattaagtttgcacgaggacgaggacaaggaaaacgatggagactctttgccaacgtttaaatgtgtgtcaggacaaaatactaacacattatgaaaatgatagtacagacctacgtgaccatatagactattggaaacacatgcgcctagaatgtgctatttattacaaggccagagaaatgggatttaaacatattaaccaccaggtggtgccaacactggctgtatcaaagaataaagcattacaagcaattgaactgcaactaacgttagaaacaatatataactcacaatatagtaatgaaaagtggacattacaagacgttagccttgaagtgtatttaactgcaccaacaggatgtataaaaaaacatggatatacagtggaagtgcagtttgatggagacatatgcaatacaatgcattatacaaactggacacatatatatatttgtgaagaagcatcagtaactgtggtagagggtcaagttgactattatggtttatattatgttcatgaaggaatacgaacatattttgtgcagtttaaagatgatgcagaaaaatatagtaaaaataaagtatgggaagttcatgcgggtggtcaggtaatattatgtcctacatctgtgtttagcagcaacgaagtatcctctcctgaaattattaggcagcacttggccaaccaccccgccgcgacccataccaaagccgtcgccttgggcaccgaagaaacacagacgactatccagcgaccaagatcagagccagacaccggaaacccctgccacaccactaagttgttgcacagagactcagtggacagtgctccaatcctcactgcatttaacagctcacacaaaggacggattaactgtaatagtaacactacacccatagtacatttaaaaggtgatgctaatactttaaaatgtttaagatatagatttaaaaagcattgtacattgtatactgcagtgtcgtctacatggcattggacaggacataatgtaaaacataaaagtgcaattgttacacttacatatgatagtgaatggcaacgtgaccaatttttgtctcaagttaaaataccaaaaactattacagtgtctactggatttatgtctatatgacaaatcttgatactgcatccacaacattactggcgtgctttttgctttgcttttgtgtgcttttgtgtgtctgcctattaatacgtccgctgcttttgtctgtgtctacatacacatcattaataatattggtattactattgtggataacagcagcctctgcgtttaggtgttttattgtatatattatatttgtttatataccattatttttaatacatacacatgcacgctttttaattacataatgtatatgtacataatgtaattgttacatataattgttgtataccataacttactattttttcttttttattttcatatataatttttttttttgtttgtttgtttgttttttaataaactgttattacttaacaatgcgacacaaacgttctgcaaaacgcacaaaacgtgcatcggctacccaactttataaaacatgcaaacaggcaggtacatgtccacctgacattatacctaaggttgaaggcaaaactattgctgatcaaatattacaatatggaagtatgggtgtattttttggtgggttaggaattggaacagggtcgggtacaggcggacgcactgggtatattccattgggaacaaggcctcccacagctacagatacacttgctcctgtaagaccccctttaacagtagatcctgtgggcccttctgatccttctatagtttctttagtggaagaaactagttttattgatgctggtgcaccaacatctgtaccttccattcccccagatgtatcaggatttagtattactacttcaactgataccacacctgctatattagatattaataatactgttactactgttactacacataataatcccactttcactgacccatctgtattgcagcctccaacacctgcagaaactggagggcattttacactttcatcatccactattagtacacataattatgaagaaattcctatggatacatttattgttagcacaaaccctaacacagtaactagtagcacacccataccagggtctcgcccagtggcacgcctaggattatatagtcgcacaacacaacaggttaaagttgtagaccctgcttttgtaaccactcccactaaacttattacatatgataatcctgcatatgaaggtatagatgtggataatacattatatttttctagtaatgataatagtattaatatagctccagatcctgactttttggatatagttgctttacataggccagcattaacctctaggcgtactggcattaggtacagtagaattggtaataaacaaacactacgtactcgtagtggaaaatctataggtgctaaggtacattattattatgatttaagtactattgatcctgcagaagaaatagaattacaaactataacaccttctacatatactaccacttcacatgcagcctcacctacttctattaataatggattatatgatatttatgcagatgactttattacagatacttctacaaccccggtaccatctgtaccctctacatctttatcaggttatattcctgcaaatacaacaattccttttggtggtgcatacaatattcctttagtatcaggtcctgatatacccattaatataactgaccaagctccttcattaattcctatagttccagggtctccacaatatacaattattgctgatgcaggtgacttttatttacatcctagttattacatgttacgaaaacgacgtaaacgtttaccatattttttttcagatgtctctttggctgcctagtgaggccactgtctacttgcctcctgtcccagtatctaaggttgtaagcacggatgaatatgttgcacgcacaaacatatattatcatgcaggaacatccagactacttgcagttggacatccctattttcctattaaaaaacctaacaataacaaaatattagttcctaaagtatcaggattacaatacagggtatttagaatacatttacctgaccccaataagtttggttttcctgacacctcattttataatccagatacacagcggctggtttgggcctgtgtaggtgttgaggtaggtcgtggtcagccattaggtgtgggcattagtggccatcctttattaaataaattggatgacacagaaaatgctagtgcttatgcagcaaatgcaggtgtggataatagagaatgtatatctatggattacaaacaaacacaattgtgtttaattggttgcaaaccacctataggggaacactggggcaaaggatccccatgtaccaatgttgcagtaaatccaggtgattgtccaccattagagttaataaacacagttattcaggatggtgatatggttgatactggctttggtgctatggactttactacattacaggctaacaaaagtgaagttccactggatatttgtacatctatttgcaaatatccagattatattaaaatggtgtcagaaccatatggcgacagcttatttttttatttacgaagggaacaaatgtttgttagacatttatttaatagggctggtactgttggtgaaaatgtaccagacgatttatacattaaaggctctgggtctactgcaaatttagccagttcaaattattttcctacacctagtggttctatggttacctctgatgcccaaatattcaataaaccttattggttacaacgagcacagggccacaataatggcatttgttggggtaaccaactatttgttactgttgttgatactacacgcagtacaaatatgtcattatgtgctgccatatctacttcagaaactacatataaaaatactaactttaaggagtacctacgacatggggaggaatatgatttacagtttatttttcaactgtgcaaaataaccttaactgcagacgttatgacatacatacattctatgaattccactattttggaggactggaattttggtctacaacctcccccaggaggcacactagaagatacttataggtttgtaacatcccaggcaattgcttgtcaaaaacatacacctccagcacctaaagaagatccccttaaaaaatacactttttgggaagtaaatttaaaggaaaagttttctgcagacctagatcagtttcctttaggacgcaaatttttactacaagcaggattgaaggccaaaccaaaatttacattaggaaaacgaaaagctacacccaccacctcatctacctctacaactgctaaacgcaaaaaacgtaagctgtaa'
    }
    args = {
        'out_dir': './puma_out_test',
        'data_dir': 'data_dir',
        'input_format': 'genbank',
        'min_prot_len': 25,
        'e_value': 1e-05,
        'for_user_dir': './puma_out_test/for_user',
        'program_files_dir': './puma_out_test/program_files'
    }
    found = {
        'E5_ALPHA': [
            4600, 4851,
            'atgacaaatcttgatactgcatccacaacattactggcgtgctttttgctttgcttttgtgtgcttttgtgtgtctgcctattaatacgtccgctgcttttgtctgtgtctacatacacatcattaataatattggtattactattgtggataacagcagcctctgcgtttaggtgttttattgtatatattatatttgtttatataccattatttttaatacatacacatgcacgctttttaattacataa',
            Seq('MTNLDTASTTLLACFLLCFCVLLCVCLLIRPLLLSVSTYTSLIILVLLLWITAA...IT*')
        ]
    }
    res = puma.identify_e5_variants(virus, args)
    assert res.keys() == found.keys()


# --------------------------------------------------
def test_blast_verify_e6():
    """Docstring"""

    warnings.simplefilter('ignore', BiopythonWarning)
    virus = {
        'E7': [
            1312, 1608,
            'atgcatggagatacacctacattgcatgaatatatgttagatttgcaaccagagacaactgatctctactgttatgagcaattaaatgacagctcagaggaggaggatgaaatagatggtccagctggacaagcagaaccggacagagcccattacaatattgtaaccttttgttgcaagtgtgactctacgcttcggttgtgcgtacaaagcacacacgtagacattcgtactttggaagacctgttaatgggcacactaggaattgtgtgccccatctgttctcagaaaccataa',
            Seq('MHGDTPTLHEYMLDLQPETTDLYCYEQLNDSSEEEDEIDGPAGQAEPDRAHYNI...KP*')
        ],
        'E1': [
            1615, 3564,
            'atggctgatcctgcaggtaccaatggggaagagggtacgggatgtaatggatggttttatgtagaggctgtagtggaaaaaaaaacaggggatgctatatcagatgacgagaacgaaaatgacagtgatacaggtgaagatttggtagattttatagtaaatgataatgattatttaacacaggcagaaacagagacagcacatgcgttgtttactgcacaggaagcaaaacaacatagagatgcagtacaggttctaaaacgaaagtatttgggtagtccacttagtgatattagtggatgtgtagacaataatattagtcctagattaaaagctatatgtatagaaaaacaaagtagagctgcaaaaaggagattatttgaaagcgaagacagcgggtatggcaatactgaagtggaaactcagcagatgttacaggtagaagggcgccatgagactgaaacaccatgtagtcagtatagtggtggaagtgggggtggttgcagtcagtacagtagtggaagtgggggagagggtgttagtgaaagacacactatatgccaaacaccacttacaaatattttaaatgtactaaaaactagtaatgcaaaggcagcaatgttagcaaaatttaaagagttatacggggtgagtttttcagaattagtaagaccatttaaaagtaataaatcaacgtgttgcgattggtgtattgctgcatttggacttacacccagtatagctgacagtataaaaacactattacaacaatattgtttatatttacacattcaaagtttagcatgttcatggggaatggttgtgttactattagtaagatataaatgtggaaaaaatagagaaacaattgaaaaattgctgtctaaactattatgtgtgtctccaatgtgtatgatgatagagcctccaaaattgcgtagtacagcagcagcattatattggtataaaacaggtatatcaaatattagtgaagtgtatggagacacgccagaatggatacaaagacaaacagtattacaacatagttttaatgattgtacatttgaattatcacagatggtacaatgggcctacgataatgacatagtagacgatagtgaaattgcatataaatatgcacaattggcagacactaatagtaatgcaagtgcctttctaaaaagtaattcacaggcaaaaattgtaaaggattgtgcaacaatgtgtagacattataaacgagcagaaaaaaaacaaatgagtatgagtcaatggataaaatatagatgtgatagggtagatgatggaggtgattggaagcaaattgttatgtttttaaggtatcaaggtgtagagtttatgtcatttttaactgcattaaaaagatttttgcaaggcatacctaaaaaaaattgcatattactatatggtgcagctaacacaggtaaatcattatttggtatgagtttaatgaaatttctgcaagggtctgtaatatgttttgtaaattctaaaagccatttttggttacaaccattagcagatgccaaaataggtatgttagatgatgctacagtgccctgttggaactacatagatgacaatttaagaaatgcattggatggaaatttagtttctatggatgtaaagcatagaccattggtacaactaaaatgccctccattattaattacatctaacattaatgctggtacagattctaggtggccttatttacataatagattggtggtgtttacatttcctaatgagtttccatttgacgaaaacggaaatccagtgtatgagcttaatgataagaactggaaatcctttttctcaaggacgtggtccagattaagtttgcacgaggacgaggacaaggaaaacgatggagactctttgccaacgtttaaatgtgtgtcaggacaaaatactaacacattatga',
            Seq('MADPAGTNGEEGTGCNGWFYVEAVVEKKTGDAISDDENENDSDTGEDLVDFIVN...TL*')
        ],
        'L2': [
            4987, 6408,
            'atgcgacacaaacgttctgcaaaacgcacaaaacgtgcatcggctacccaactttataaaacatgcaaacaggcaggtacatgtccacctgacattatacctaaggttgaaggcaaaactattgctgatcaaatattacaatatggaagtatgggtgtattttttggtgggttaggaattggaacagggtcgggtacaggcggacgcactgggtatattccattgggaacaaggcctcccacagctacagatacacttgctcctgtaagaccccctttaacagtagatcctgtgggcccttctgatccttctatagtttctttagtggaagaaactagttttattgatgctggtgcaccaacatctgtaccttccattcccccagatgtatcaggatttagtattactacttcaactgataccacacctgctatattagatattaataatactgttactactgttactacacataataatcccactttcactgacccatctgtattgcagcctccaacacctgcagaaactggagggcattttacactttcatcatccactattagtacacataattatgaagaaattcctatggatacatttattgttagcacaaaccctaacacagtaactagtagcacacccataccagggtctcgcccagtggcacgcctaggattatatagtcgcacaacacaacaggttaaagttgtagaccctgcttttgtaaccactcccactaaacttattacatatgataatcctgcatatgaaggtatagatgtggataatacattatatttttctagtaatgataatagtattaatatagctccagatcctgactttttggatatagttgctttacataggccagcattaacctctaggcgtactggcattaggtacagtagaattggtaataaacaaacactacgtactcgtagtggaaaatctataggtgctaaggtacattattattatgatttaagtactattgatcctgcagaagaaatagaattacaaactataacaccttctacatatactaccacttcacatgcagcctcacctacttctattaataatggattatatgatatttatgcagatgactttattacagatacttctacaaccccggtaccatctgtaccctctacatctttatcaggttatattcctgcaaatacaacaattccttttggtggtgcatacaatattcctttagtatcaggtcctgatatacccattaatataactgaccaagctccttcattaattcctatagttccagggtctccacaatatacaattattgctgatgcaggtgacttttatttacatcctagttattacatgttacgaaaacgacgtaaacgtttaccatattttttttcagatgtctctttggctgcctag',
            Seq('MRHKRSAKRTKRASATQLYKTCKQAGTCPPDIIPKVEGKTIADQILQYGSMGVF...AA*')
        ],
        'E6': [
            833, 1309,
            'atgcaccaaaagagaactgcaatgtttcaggacccacaggagcgacccagaaagttaccacagttatgcacagagctgcaaacaactatacatgatataatattagaatgtgtgtactgcaagcaacagttactgcgacgtgaggtatatgactttgcttttcgggatttatgcatagtatatagagatgggaatccatatgctgtatgtgataaatgtttaaagttttattctaaaattagtgagtatagacattattgttatagtttgtatggaacaacattagaacagcaatacaacaaaccgttgtgtgatttgttaattaggtgtattaactgtcaaaagccactgtgtcctgaagaaaagcaaagacatctggacaaaaagcaaagattccataatataaggggtcggtggaccggtcgatgtatgtcttgttgcagatcatcaagaacacgtagagaaacccagctgtaa',
            Seq('MHQKRTAMFQDPQERPRKLPQLCTELQTTIHDIILECVYCKQQLLRREVYDFAF...QL*')
        ],
        'E2': [
            3506, 4603,
            'atggagactctttgccaacgtttaaatgtgtgtcaggacaaaatactaacacattatgaaaatgatagtacagacctacgtgaccatatagactattggaaacacatgcgcctagaatgtgctatttattacaaggccagagaaatgggatttaaacatattaaccaccaggtggtgccaacactggctgtatcaaagaataaagcattacaagcaattgaactgcaactaacgttagaaacaatatataactcacaatatagtaatgaaaagtggacattacaagacgttagccttgaagtgtatttaactgcaccaacaggatgtataaaaaaacatggatatacagtggaagtgcagtttgatggagacatatgcaatacaatgcattatacaaactggacacatatatatatttgtgaagaagcatcagtaactgtggtagagggtcaagttgactattatggtttatattatgttcatgaaggaatacgaacatattttgtgcagtttaaagatgatgcagaaaaatatagtaaaaataaagtatgggaagttcatgcgggtggtcaggtaatattatgtcctacatctgtgtttagcagcaacgaagtatcctctcctgaaattattaggcagcacttggccaaccaccccgccgcgacccataccaaagccgtcgccttgggcaccgaagaaacacagacgactatccagcgaccaagatcagagccagacaccggaaacccctgccacaccactaagttgttgcacagagactcagtggacagtgctccaatcctcactgcatttaacagctcacacaaaggacggattaactgtaatagtaacactacacccatagtacatttaaaaggtgatgctaatactttaaaatgtttaagatatagatttaaaaagcattgtacattgtatactgcagtgtcgtctacatggcattggacaggacataatgtaaaacataaaagtgcaattgttacacttacatatgatagtgaatggcaacgtgaccaatttttgtctcaagttaaaataccaaaaactattacagtgtctactggatttatgtctatatga',
            Seq('METLCQRLNVCQDKILTHYENDSTDLRDHIDYWKHMRLECAIYYKAREMGFKHI...SI*')
        ],
        'L1': [
            6389, 7906,
            'atgtctctttggctgcctagtgaggccactgtctacttgcctcctgtcccagtatctaaggttgtaagcacggatgaatatgttgcacgcacaaacatatattatcatgcaggaacatccagactacttgcagttggacatccctattttcctattaaaaaacctaacaataacaaaatattagttcctaaagtatcaggattacaatacagggtatttagaatacatttacctgaccccaataagtttggttttcctgacacctcattttataatccagatacacagcggctggtttgggcctgtgtaggtgttgaggtaggtcgtggtcagccattaggtgtgggcattagtggccatcctttattaaataaattggatgacacagaaaatgctagtgcttatgcagcaaatgcaggtgtggataatagagaatgtatatctatggattacaaacaaacacaattgtgtttaattggttgcaaaccacctataggggaacactggggcaaaggatccccatgtaccaatgttgcagtaaatccaggtgattgtccaccattagagttaataaacacagttattcaggatggtgatatggttgatactggctttggtgctatggactttactacattacaggctaacaaaagtgaagttccactggatatttgtacatctatttgcaaatatccagattatattaaaatggtgtcagaaccatatggcgacagcttatttttttatttacgaagggaacaaatgtttgttagacatttatttaatagggctggtactgttggtgaaaatgtaccagacgatttatacattaaaggctctgggtctactgcaaatttagccagttcaaattattttcctacacctagtggttctatggttacctctgatgcccaaatattcaataaaccttattggttacaacgagcacagggccacaataatggcatttgttggggtaaccaactatttgttactgttgttgatactacacgcagtacaaatatgtcattatgtgctgccatatctacttcagaaactacatataaaaatactaactttaaggagtacctacgacatggggaggaatatgatttacagtttatttttcaactgtgcaaaataaccttaactgcagacgttatgacatacatacattctatgaattccactattttggaggactggaattttggtctacaacctcccccaggaggcacactagaagatacttataggtttgtaacatcccaggcaattgcttgtcaaaaacatacacctccagcacctaaagaagatccccttaaaaaatacactttttgggaagtaaatttaaaggaaaagttttctgcagacctagatcagtttcctttaggacgcaaatttttactacaagcaggattgaaggccaaaccaaaatttacattaggaaaacgaaaagctacacccaccacctcatctacctctacaactgctaaacgcaaaaaacgtaagctgtaa',
            Seq('MSLWLPSEATVYLPPVPVSKVVSTDEYVARTNIYYHAGTSRLLAVGHPYFPIKK...KL*')
        ],
        'name':
        'Human papillomavirus 16 (HPV16)',
        'accession':
        'HPV16REF',
        'genome':
        'gtattgtatgtatgttgaattagtgttgtttgttgtgtatatgtttgtatgtgcttgtatgtgcttgtaaatattaagttgtatgtgtgtttgtatgtatggtataataaacacgtgtgtatgtgtttttaaatgcttgtgtaactattgtgtcatgcaacataaataaacttattgtttcaacacctactaattgtgttgtggttattcattgtatataaactatatttgctacatcctgtttttgttttatatatactatattttgtagcgccagcggccattttgtagcttcaaccgaattcggttgcatgctttttggcacaaaatgtgtttttttaaatagttctatgtcagcaactatggtttaaacttgtacgtttcctgcttgccatgcgtgccaaatccctgttttcctgacctgcactgcttgccaaccattccattgttttttacactgcactatgtgcaactactgaatcactatgtacattgtgtcatataaaataaatcactatgcgccaacgccttacataccgctgttaggcacatatttttggcttgttttaactaacctaattgcatatttggcataaggtttaaacttctaaggccaactaaatgtcaccctagttcatacatgaactgtgtaaaggttagtcatacattgttcatttgtaaaactgcacatgggtgtgtgcaaaccgttttgggttacacatttacaagcaacttatataataatactaaactacaataattcatgtataaaactaagggcgtaaccgaaatcggttgaaccgaaaccggttagtataaaagcagacattttatgcaccaaaagagaactgcaatgtttcaggacccacaggagcgacccagaaagttaccacagttatgcacagagctgcaaacaactatacatgatataatattagaatgtgtgtactgcaagcaacagttactgcgacgtgaggtatatgactttgcttttcgggatttatgcatagtatatagagatgggaatccatatgctgtatgtgataaatgtttaaagttttattctaaaattagtgagtatagacattattgttatagtttgtatggaacaacattagaacagcaatacaacaaaccgttgtgtgatttgttaattaggtgtattaactgtcaaaagccactgtgtcctgaagaaaagcaaagacatctggacaaaaagcaaagattccataatataaggggtcggtggaccggtcgatgtatgtcttgttgcagatcatcaagaacacgtagagaaacccagctgtaatcatgcatggagatacacctacattgcatgaatatatgttagatttgcaaccagagacaactgatctctactgttatgagcaattaaatgacagctcagaggaggaggatgaaatagatggtccagctggacaagcagaaccggacagagcccattacaatattgtaaccttttgttgcaagtgtgactctacgcttcggttgtgcgtacaaagcacacacgtagacattcgtactttggaagacctgttaatgggcacactaggaattgtgtgccccatctgttctcagaaaccataatctaccatggctgatcctgcaggtaccaatggggaagagggtacgggatgtaatggatggttttatgtagaggctgtagtggaaaaaaaaacaggggatgctatatcagatgacgagaacgaaaatgacagtgatacaggtgaagatttggtagattttatagtaaatgataatgattatttaacacaggcagaaacagagacagcacatgcgttgtttactgcacaggaagcaaaacaacatagagatgcagtacaggttctaaaacgaaagtatttgggtagtccacttagtgatattagtggatgtgtagacaataatattagtcctagattaaaagctatatgtatagaaaaacaaagtagagctgcaaaaaggagattatttgaaagcgaagacagcgggtatggcaatactgaagtggaaactcagcagatgttacaggtagaagggcgccatgagactgaaacaccatgtagtcagtatagtggtggaagtgggggtggttgcagtcagtacagtagtggaagtgggggagagggtgttagtgaaagacacactatatgccaaacaccacttacaaatattttaaatgtactaaaaactagtaatgcaaaggcagcaatgttagcaaaatttaaagagttatacggggtgagtttttcagaattagtaagaccatttaaaagtaataaatcaacgtgttgcgattggtgtattgctgcatttggacttacacccagtatagctgacagtataaaaacactattacaacaatattgtttatatttacacattcaaagtttagcatgttcatggggaatggttgtgttactattagtaagatataaatgtggaaaaaatagagaaacaattgaaaaattgctgtctaaactattatgtgtgtctccaatgtgtatgatgatagagcctccaaaattgcgtagtacagcagcagcattatattggtataaaacaggtatatcaaatattagtgaagtgtatggagacacgccagaatggatacaaagacaaacagtattacaacatagttttaatgattgtacatttgaattatcacagatggtacaatgggcctacgataatgacatagtagacgatagtgaaattgcatataaatatgcacaattggcagacactaatagtaatgcaagtgcctttctaaaaagtaattcacaggcaaaaattgtaaaggattgtgcaacaatgtgtagacattataaacgagcagaaaaaaaacaaatgagtatgagtcaatggataaaatatagatgtgatagggtagatgatggaggtgattggaagcaaattgttatgtttttaaggtatcaaggtgtagagtttatgtcatttttaactgcattaaaaagatttttgcaaggcatacctaaaaaaaattgcatattactatatggtgcagctaacacaggtaaatcattatttggtatgagtttaatgaaatttctgcaagggtctgtaatatgttttgtaaattctaaaagccatttttggttacaaccattagcagatgccaaaataggtatgttagatgatgctacagtgccctgttggaactacatagatgacaatttaagaaatgcattggatggaaatttagtttctatggatgtaaagcatagaccattggtacaactaaaatgccctccattattaattacatctaacattaatgctggtacagattctaggtggccttatttacataatagattggtggtgtttacatttcctaatgagtttccatttgacgaaaacggaaatccagtgtatgagcttaatgataagaactggaaatcctttttctcaaggacgtggtccagattaagtttgcacgaggacgaggacaaggaaaacgatggagactctttgccaacgtttaaatgtgtgtcaggacaaaatactaacacattatgaaaatgatagtacagacctacgtgaccatatagactattggaaacacatgcgcctagaatgtgctatttattacaaggccagagaaatgggatttaaacatattaaccaccaggtggtgccaacactggctgtatcaaagaataaagcattacaagcaattgaactgcaactaacgttagaaacaatatataactcacaatatagtaatgaaaagtggacattacaagacgttagccttgaagtgtatttaactgcaccaacaggatgtataaaaaaacatggatatacagtggaagtgcagtttgatggagacatatgcaatacaatgcattatacaaactggacacatatatatatttgtgaagaagcatcagtaactgtggtagagggtcaagttgactattatggtttatattatgttcatgaaggaatacgaacatattttgtgcagtttaaagatgatgcagaaaaatatagtaaaaataaagtatgggaagttcatgcgggtggtcaggtaatattatgtcctacatctgtgtttagcagcaacgaagtatcctctcctgaaattattaggcagcacttggccaaccaccccgccgcgacccataccaaagccgtcgccttgggcaccgaagaaacacagacgactatccagcgaccaagatcagagccagacaccggaaacccctgccacaccactaagttgttgcacagagactcagtggacagtgctccaatcctcactgcatttaacagctcacacaaaggacggattaactgtaatagtaacactacacccatagtacatttaaaaggtgatgctaatactttaaaatgtttaagatatagatttaaaaagcattgtacattgtatactgcagtgtcgtctacatggcattggacaggacataatgtaaaacataaaagtgcaattgttacacttacatatgatagtgaatggcaacgtgaccaatttttgtctcaagttaaaataccaaaaactattacagtgtctactggatttatgtctatatgacaaatcttgatactgcatccacaacattactggcgtgctttttgctttgcttttgtgtgcttttgtgtgtctgcctattaatacgtccgctgcttttgtctgtgtctacatacacatcattaataatattggtattactattgtggataacagcagcctctgcgtttaggtgttttattgtatatattatatttgtttatataccattatttttaatacatacacatgcacgctttttaattacataatgtatatgtacataatgtaattgttacatataattgttgtataccataacttactattttttcttttttattttcatatataatttttttttttgtttgtttgtttgttttttaataaactgttattacttaacaatgcgacacaaacgttctgcaaaacgcacaaaacgtgcatcggctacccaactttataaaacatgcaaacaggcaggtacatgtccacctgacattatacctaaggttgaaggcaaaactattgctgatcaaatattacaatatggaagtatgggtgtattttttggtgggttaggaattggaacagggtcgggtacaggcggacgcactgggtatattccattgggaacaaggcctcccacagctacagatacacttgctcctgtaagaccccctttaacagtagatcctgtgggcccttctgatccttctatagtttctttagtggaagaaactagttttattgatgctggtgcaccaacatctgtaccttccattcccccagatgtatcaggatttagtattactacttcaactgataccacacctgctatattagatattaataatactgttactactgttactacacataataatcccactttcactgacccatctgtattgcagcctccaacacctgcagaaactggagggcattttacactttcatcatccactattagtacacataattatgaagaaattcctatggatacatttattgttagcacaaaccctaacacagtaactagtagcacacccataccagggtctcgcccagtggcacgcctaggattatatagtcgcacaacacaacaggttaaagttgtagaccctgcttttgtaaccactcccactaaacttattacatatgataatcctgcatatgaaggtatagatgtggataatacattatatttttctagtaatgataatagtattaatatagctccagatcctgactttttggatatagttgctttacataggccagcattaacctctaggcgtactggcattaggtacagtagaattggtaataaacaaacactacgtactcgtagtggaaaatctataggtgctaaggtacattattattatgatttaagtactattgatcctgcagaagaaatagaattacaaactataacaccttctacatatactaccacttcacatgcagcctcacctacttctattaataatggattatatgatatttatgcagatgactttattacagatacttctacaaccccggtaccatctgtaccctctacatctttatcaggttatattcctgcaaatacaacaattccttttggtggtgcatacaatattcctttagtatcaggtcctgatatacccattaatataactgaccaagctccttcattaattcctatagttccagggtctccacaatatacaattattgctgatgcaggtgacttttatttacatcctagttattacatgttacgaaaacgacgtaaacgtttaccatattttttttcagatgtctctttggctgcctagtgaggccactgtctacttgcctcctgtcccagtatctaaggttgtaagcacggatgaatatgttgcacgcacaaacatatattatcatgcaggaacatccagactacttgcagttggacatccctattttcctattaaaaaacctaacaataacaaaatattagttcctaaagtatcaggattacaatacagggtatttagaatacatttacctgaccccaataagtttggttttcctgacacctcattttataatccagatacacagcggctggtttgggcctgtgtaggtgttgaggtaggtcgtggtcagccattaggtgtgggcattagtggccatcctttattaaataaattggatgacacagaaaatgctagtgcttatgcagcaaatgcaggtgtggataatagagaatgtatatctatggattacaaacaaacacaattgtgtttaattggttgcaaaccacctataggggaacactggggcaaaggatccccatgtaccaatgttgcagtaaatccaggtgattgtccaccattagagttaataaacacagttattcaggatggtgatatggttgatactggctttggtgctatggactttactacattacaggctaacaaaagtgaagttccactggatatttgtacatctatttgcaaatatccagattatattaaaatggtgtcagaaccatatggcgacagcttatttttttatttacgaagggaacaaatgtttgttagacatttatttaatagggctggtactgttggtgaaaatgtaccagacgatttatacattaaaggctctgggtctactgcaaatttagccagttcaaattattttcctacacctagtggttctatggttacctctgatgcccaaatattcaataaaccttattggttacaacgagcacagggccacaataatggcatttgttggggtaaccaactatttgttactgttgttgatactacacgcagtacaaatatgtcattatgtgctgccatatctacttcagaaactacatataaaaatactaactttaaggagtacctacgacatggggaggaatatgatttacagtttatttttcaactgtgcaaaataaccttaactgcagacgttatgacatacatacattctatgaattccactattttggaggactggaattttggtctacaacctcccccaggaggcacactagaagatacttataggtttgtaacatcccaggcaattgcttgtcaaaaacatacacctccagcacctaaagaagatccccttaaaaaatacactttttgggaagtaaatttaaaggaaaagttttctgcagacctagatcagtttcctttaggacgcaaatttttactacaagcaggattgaaggccaaaccaaaatttacattaggaaaacgaaaagctacacccaccacctcatctacctctacaactgctaaacgcaaaaaacgtaagctgtaa'
    }
    args = {
        'out_dir': './puma_out_test',
        'data_dir': 'data_dir',
        'input_format': 'genbank',
        'min_prot_len': 25,
        'e_value': 1e-05,
        'for_user_dir': './puma_out_test/for_user',
        'program_files_dir': './puma_out_test/program_files'
    }
    blast_out = os.path.join('data_dir', 'tests', 'verify_E6',
                             'blast_result_E6.tab')
    assert os.path.isfile(blast_out)
    res = puma.blast_verify_e6(virus, args)
    assert os.path.isfile(res)


# --------------------------------------------------
def test_parse_blast_results_verify_e6():
    """Docstring"""

    warnings.simplefilter('ignore', BiopythonWarning)
    virus = {
        'E7': [
            1312, 1608,
            'atgcatggagatacacctacattgcatgaatatatgttagatttgcaaccagagacaactgatctctactgttatgagcaattaaatgacagctcagaggaggaggatgaaatagatggtccagctggacaagcagaaccggacagagcccattacaatattgtaaccttttgttgcaagtgtgactctacgcttcggttgtgcgtacaaagcacacacgtagacattcgtactttggaagacctgttaatgggcacactaggaattgtgtgccccatctgttctcagaaaccataa',
            Seq('MHGDTPTLHEYMLDLQPETTDLYCYEQLNDSSEEEDEIDGPAGQAEPDRAHYNI...KP*')
        ],
        'E1': [
            1615, 3564,
            'atggctgatcctgcaggtaccaatggggaagagggtacgggatgtaatggatggttttatgtagaggctgtagtggaaaaaaaaacaggggatgctatatcagatgacgagaacgaaaatgacagtgatacaggtgaagatttggtagattttatagtaaatgataatgattatttaacacaggcagaaacagagacagcacatgcgttgtttactgcacaggaagcaaaacaacatagagatgcagtacaggttctaaaacgaaagtatttgggtagtccacttagtgatattagtggatgtgtagacaataatattagtcctagattaaaagctatatgtatagaaaaacaaagtagagctgcaaaaaggagattatttgaaagcgaagacagcgggtatggcaatactgaagtggaaactcagcagatgttacaggtagaagggcgccatgagactgaaacaccatgtagtcagtatagtggtggaagtgggggtggttgcagtcagtacagtagtggaagtgggggagagggtgttagtgaaagacacactatatgccaaacaccacttacaaatattttaaatgtactaaaaactagtaatgcaaaggcagcaatgttagcaaaatttaaagagttatacggggtgagtttttcagaattagtaagaccatttaaaagtaataaatcaacgtgttgcgattggtgtattgctgcatttggacttacacccagtatagctgacagtataaaaacactattacaacaatattgtttatatttacacattcaaagtttagcatgttcatggggaatggttgtgttactattagtaagatataaatgtggaaaaaatagagaaacaattgaaaaattgctgtctaaactattatgtgtgtctccaatgtgtatgatgatagagcctccaaaattgcgtagtacagcagcagcattatattggtataaaacaggtatatcaaatattagtgaagtgtatggagacacgccagaatggatacaaagacaaacagtattacaacatagttttaatgattgtacatttgaattatcacagatggtacaatgggcctacgataatgacatagtagacgatagtgaaattgcatataaatatgcacaattggcagacactaatagtaatgcaagtgcctttctaaaaagtaattcacaggcaaaaattgtaaaggattgtgcaacaatgtgtagacattataaacgagcagaaaaaaaacaaatgagtatgagtcaatggataaaatatagatgtgatagggtagatgatggaggtgattggaagcaaattgttatgtttttaaggtatcaaggtgtagagtttatgtcatttttaactgcattaaaaagatttttgcaaggcatacctaaaaaaaattgcatattactatatggtgcagctaacacaggtaaatcattatttggtatgagtttaatgaaatttctgcaagggtctgtaatatgttttgtaaattctaaaagccatttttggttacaaccattagcagatgccaaaataggtatgttagatgatgctacagtgccctgttggaactacatagatgacaatttaagaaatgcattggatggaaatttagtttctatggatgtaaagcatagaccattggtacaactaaaatgccctccattattaattacatctaacattaatgctggtacagattctaggtggccttatttacataatagattggtggtgtttacatttcctaatgagtttccatttgacgaaaacggaaatccagtgtatgagcttaatgataagaactggaaatcctttttctcaaggacgtggtccagattaagtttgcacgaggacgaggacaaggaaaacgatggagactctttgccaacgtttaaatgtgtgtcaggacaaaatactaacacattatga',
            Seq('MADPAGTNGEEGTGCNGWFYVEAVVEKKTGDAISDDENENDSDTGEDLVDFIVN...TL*')
        ],
        'L2': [
            4987, 6408,
            'atgcgacacaaacgttctgcaaaacgcacaaaacgtgcatcggctacccaactttataaaacatgcaaacaggcaggtacatgtccacctgacattatacctaaggttgaaggcaaaactattgctgatcaaatattacaatatggaagtatgggtgtattttttggtgggttaggaattggaacagggtcgggtacaggcggacgcactgggtatattccattgggaacaaggcctcccacagctacagatacacttgctcctgtaagaccccctttaacagtagatcctgtgggcccttctgatccttctatagtttctttagtggaagaaactagttttattgatgctggtgcaccaacatctgtaccttccattcccccagatgtatcaggatttagtattactacttcaactgataccacacctgctatattagatattaataatactgttactactgttactacacataataatcccactttcactgacccatctgtattgcagcctccaacacctgcagaaactggagggcattttacactttcatcatccactattagtacacataattatgaagaaattcctatggatacatttattgttagcacaaaccctaacacagtaactagtagcacacccataccagggtctcgcccagtggcacgcctaggattatatagtcgcacaacacaacaggttaaagttgtagaccctgcttttgtaaccactcccactaaacttattacatatgataatcctgcatatgaaggtatagatgtggataatacattatatttttctagtaatgataatagtattaatatagctccagatcctgactttttggatatagttgctttacataggccagcattaacctctaggcgtactggcattaggtacagtagaattggtaataaacaaacactacgtactcgtagtggaaaatctataggtgctaaggtacattattattatgatttaagtactattgatcctgcagaagaaatagaattacaaactataacaccttctacatatactaccacttcacatgcagcctcacctacttctattaataatggattatatgatatttatgcagatgactttattacagatacttctacaaccccggtaccatctgtaccctctacatctttatcaggttatattcctgcaaatacaacaattccttttggtggtgcatacaatattcctttagtatcaggtcctgatatacccattaatataactgaccaagctccttcattaattcctatagttccagggtctccacaatatacaattattgctgatgcaggtgacttttatttacatcctagttattacatgttacgaaaacgacgtaaacgtttaccatattttttttcagatgtctctttggctgcctag',
            Seq('MRHKRSAKRTKRASATQLYKTCKQAGTCPPDIIPKVEGKTIADQILQYGSMGVF...AA*')
        ],
        'E6': [
            833, 1309,
            'atgcaccaaaagagaactgcaatgtttcaggacccacaggagcgacccagaaagttaccacagttatgcacagagctgcaaacaactatacatgatataatattagaatgtgtgtactgcaagcaacagttactgcgacgtgaggtatatgactttgcttttcgggatttatgcatagtatatagagatgggaatccatatgctgtatgtgataaatgtttaaagttttattctaaaattagtgagtatagacattattgttatagtttgtatggaacaacattagaacagcaatacaacaaaccgttgtgtgatttgttaattaggtgtattaactgtcaaaagccactgtgtcctgaagaaaagcaaagacatctggacaaaaagcaaagattccataatataaggggtcggtggaccggtcgatgtatgtcttgttgcagatcatcaagaacacgtagagaaacccagctgtaa',
            Seq('MHQKRTAMFQDPQERPRKLPQLCTELQTTIHDIILECVYCKQQLLRREVYDFAF...QL*')
        ],
        'E2': [
            3506, 4603,
            'atggagactctttgccaacgtttaaatgtgtgtcaggacaaaatactaacacattatgaaaatgatagtacagacctacgtgaccatatagactattggaaacacatgcgcctagaatgtgctatttattacaaggccagagaaatgggatttaaacatattaaccaccaggtggtgccaacactggctgtatcaaagaataaagcattacaagcaattgaactgcaactaacgttagaaacaatatataactcacaatatagtaatgaaaagtggacattacaagacgttagccttgaagtgtatttaactgcaccaacaggatgtataaaaaaacatggatatacagtggaagtgcagtttgatggagacatatgcaatacaatgcattatacaaactggacacatatatatatttgtgaagaagcatcagtaactgtggtagagggtcaagttgactattatggtttatattatgttcatgaaggaatacgaacatattttgtgcagtttaaagatgatgcagaaaaatatagtaaaaataaagtatgggaagttcatgcgggtggtcaggtaatattatgtcctacatctgtgtttagcagcaacgaagtatcctctcctgaaattattaggcagcacttggccaaccaccccgccgcgacccataccaaagccgtcgccttgggcaccgaagaaacacagacgactatccagcgaccaagatcagagccagacaccggaaacccctgccacaccactaagttgttgcacagagactcagtggacagtgctccaatcctcactgcatttaacagctcacacaaaggacggattaactgtaatagtaacactacacccatagtacatttaaaaggtgatgctaatactttaaaatgtttaagatatagatttaaaaagcattgtacattgtatactgcagtgtcgtctacatggcattggacaggacataatgtaaaacataaaagtgcaattgttacacttacatatgatagtgaatggcaacgtgaccaatttttgtctcaagttaaaataccaaaaactattacagtgtctactggatttatgtctatatga',
            Seq('METLCQRLNVCQDKILTHYENDSTDLRDHIDYWKHMRLECAIYYKAREMGFKHI...SI*')
        ],
        'L1': [
            6389, 7906,
            'atgtctctttggctgcctagtgaggccactgtctacttgcctcctgtcccagtatctaaggttgtaagcacggatgaatatgttgcacgcacaaacatatattatcatgcaggaacatccagactacttgcagttggacatccctattttcctattaaaaaacctaacaataacaaaatattagttcctaaagtatcaggattacaatacagggtatttagaatacatttacctgaccccaataagtttggttttcctgacacctcattttataatccagatacacagcggctggtttgggcctgtgtaggtgttgaggtaggtcgtggtcagccattaggtgtgggcattagtggccatcctttattaaataaattggatgacacagaaaatgctagtgcttatgcagcaaatgcaggtgtggataatagagaatgtatatctatggattacaaacaaacacaattgtgtttaattggttgcaaaccacctataggggaacactggggcaaaggatccccatgtaccaatgttgcagtaaatccaggtgattgtccaccattagagttaataaacacagttattcaggatggtgatatggttgatactggctttggtgctatggactttactacattacaggctaacaaaagtgaagttccactggatatttgtacatctatttgcaaatatccagattatattaaaatggtgtcagaaccatatggcgacagcttatttttttatttacgaagggaacaaatgtttgttagacatttatttaatagggctggtactgttggtgaaaatgtaccagacgatttatacattaaaggctctgggtctactgcaaatttagccagttcaaattattttcctacacctagtggttctatggttacctctgatgcccaaatattcaataaaccttattggttacaacgagcacagggccacaataatggcatttgttggggtaaccaactatttgttactgttgttgatactacacgcagtacaaatatgtcattatgtgctgccatatctacttcagaaactacatataaaaatactaactttaaggagtacctacgacatggggaggaatatgatttacagtttatttttcaactgtgcaaaataaccttaactgcagacgttatgacatacatacattctatgaattccactattttggaggactggaattttggtctacaacctcccccaggaggcacactagaagatacttataggtttgtaacatcccaggcaattgcttgtcaaaaacatacacctccagcacctaaagaagatccccttaaaaaatacactttttgggaagtaaatttaaaggaaaagttttctgcagacctagatcagtttcctttaggacgcaaatttttactacaagcaggattgaaggccaaaccaaaatttacattaggaaaacgaaaagctacacccaccacctcatctacctctacaactgctaaacgcaaaaaacgtaagctgtaa',
            Seq('MSLWLPSEATVYLPPVPVSKVVSTDEYVARTNIYYHAGTSRLLAVGHPYFPIKK...KL*')
        ],
        'name':
        'Human papillomavirus 16 (HPV16)',
        'accession':
        'HPV16REF',
        'genome':
        'gtattgtatgtatgttgaattagtgttgtttgttgtgtatatgtttgtatgtgcttgtatgtgcttgtaaatattaagttgtatgtgtgtttgtatgtatggtataataaacacgtgtgtatgtgtttttaaatgcttgtgtaactattgtgtcatgcaacataaataaacttattgtttcaacacctactaattgtgttgtggttattcattgtatataaactatatttgctacatcctgtttttgttttatatatactatattttgtagcgccagcggccattttgtagcttcaaccgaattcggttgcatgctttttggcacaaaatgtgtttttttaaatagttctatgtcagcaactatggtttaaacttgtacgtttcctgcttgccatgcgtgccaaatccctgttttcctgacctgcactgcttgccaaccattccattgttttttacactgcactatgtgcaactactgaatcactatgtacattgtgtcatataaaataaatcactatgcgccaacgccttacataccgctgttaggcacatatttttggcttgttttaactaacctaattgcatatttggcataaggtttaaacttctaaggccaactaaatgtcaccctagttcatacatgaactgtgtaaaggttagtcatacattgttcatttgtaaaactgcacatgggtgtgtgcaaaccgttttgggttacacatttacaagcaacttatataataatactaaactacaataattcatgtataaaactaagggcgtaaccgaaatcggttgaaccgaaaccggttagtataaaagcagacattttatgcaccaaaagagaactgcaatgtttcaggacccacaggagcgacccagaaagttaccacagttatgcacagagctgcaaacaactatacatgatataatattagaatgtgtgtactgcaagcaacagttactgcgacgtgaggtatatgactttgcttttcgggatttatgcatagtatatagagatgggaatccatatgctgtatgtgataaatgtttaaagttttattctaaaattagtgagtatagacattattgttatagtttgtatggaacaacattagaacagcaatacaacaaaccgttgtgtgatttgttaattaggtgtattaactgtcaaaagccactgtgtcctgaagaaaagcaaagacatctggacaaaaagcaaagattccataatataaggggtcggtggaccggtcgatgtatgtcttgttgcagatcatcaagaacacgtagagaaacccagctgtaatcatgcatggagatacacctacattgcatgaatatatgttagatttgcaaccagagacaactgatctctactgttatgagcaattaaatgacagctcagaggaggaggatgaaatagatggtccagctggacaagcagaaccggacagagcccattacaatattgtaaccttttgttgcaagtgtgactctacgcttcggttgtgcgtacaaagcacacacgtagacattcgtactttggaagacctgttaatgggcacactaggaattgtgtgccccatctgttctcagaaaccataatctaccatggctgatcctgcaggtaccaatggggaagagggtacgggatgtaatggatggttttatgtagaggctgtagtggaaaaaaaaacaggggatgctatatcagatgacgagaacgaaaatgacagtgatacaggtgaagatttggtagattttatagtaaatgataatgattatttaacacaggcagaaacagagacagcacatgcgttgtttactgcacaggaagcaaaacaacatagagatgcagtacaggttctaaaacgaaagtatttgggtagtccacttagtgatattagtggatgtgtagacaataatattagtcctagattaaaagctatatgtatagaaaaacaaagtagagctgcaaaaaggagattatttgaaagcgaagacagcgggtatggcaatactgaagtggaaactcagcagatgttacaggtagaagggcgccatgagactgaaacaccatgtagtcagtatagtggtggaagtgggggtggttgcagtcagtacagtagtggaagtgggggagagggtgttagtgaaagacacactatatgccaaacaccacttacaaatattttaaatgtactaaaaactagtaatgcaaaggcagcaatgttagcaaaatttaaagagttatacggggtgagtttttcagaattagtaagaccatttaaaagtaataaatcaacgtgttgcgattggtgtattgctgcatttggacttacacccagtatagctgacagtataaaaacactattacaacaatattgtttatatttacacattcaaagtttagcatgttcatggggaatggttgtgttactattagtaagatataaatgtggaaaaaatagagaaacaattgaaaaattgctgtctaaactattatgtgtgtctccaatgtgtatgatgatagagcctccaaaattgcgtagtacagcagcagcattatattggtataaaacaggtatatcaaatattagtgaagtgtatggagacacgccagaatggatacaaagacaaacagtattacaacatagttttaatgattgtacatttgaattatcacagatggtacaatgggcctacgataatgacatagtagacgatagtgaaattgcatataaatatgcacaattggcagacactaatagtaatgcaagtgcctttctaaaaagtaattcacaggcaaaaattgtaaaggattgtgcaacaatgtgtagacattataaacgagcagaaaaaaaacaaatgagtatgagtcaatggataaaatatagatgtgatagggtagatgatggaggtgattggaagcaaattgttatgtttttaaggtatcaaggtgtagagtttatgtcatttttaactgcattaaaaagatttttgcaaggcatacctaaaaaaaattgcatattactatatggtgcagctaacacaggtaaatcattatttggtatgagtttaatgaaatttctgcaagggtctgtaatatgttttgtaaattctaaaagccatttttggttacaaccattagcagatgccaaaataggtatgttagatgatgctacagtgccctgttggaactacatagatgacaatttaagaaatgcattggatggaaatttagtttctatggatgtaaagcatagaccattggtacaactaaaatgccctccattattaattacatctaacattaatgctggtacagattctaggtggccttatttacataatagattggtggtgtttacatttcctaatgagtttccatttgacgaaaacggaaatccagtgtatgagcttaatgataagaactggaaatcctttttctcaaggacgtggtccagattaagtttgcacgaggacgaggacaaggaaaacgatggagactctttgccaacgtttaaatgtgtgtcaggacaaaatactaacacattatgaaaatgatagtacagacctacgtgaccatatagactattggaaacacatgcgcctagaatgtgctatttattacaaggccagagaaatgggatttaaacatattaaccaccaggtggtgccaacactggctgtatcaaagaataaagcattacaagcaattgaactgcaactaacgttagaaacaatatataactcacaatatagtaatgaaaagtggacattacaagacgttagccttgaagtgtatttaactgcaccaacaggatgtataaaaaaacatggatatacagtggaagtgcagtttgatggagacatatgcaatacaatgcattatacaaactggacacatatatatatttgtgaagaagcatcagtaactgtggtagagggtcaagttgactattatggtttatattatgttcatgaaggaatacgaacatattttgtgcagtttaaagatgatgcagaaaaatatagtaaaaataaagtatgggaagttcatgcgggtggtcaggtaatattatgtcctacatctgtgtttagcagcaacgaagtatcctctcctgaaattattaggcagcacttggccaaccaccccgccgcgacccataccaaagccgtcgccttgggcaccgaagaaacacagacgactatccagcgaccaagatcagagccagacaccggaaacccctgccacaccactaagttgttgcacagagactcagtggacagtgctccaatcctcactgcatttaacagctcacacaaaggacggattaactgtaatagtaacactacacccatagtacatttaaaaggtgatgctaatactttaaaatgtttaagatatagatttaaaaagcattgtacattgtatactgcagtgtcgtctacatggcattggacaggacataatgtaaaacataaaagtgcaattgttacacttacatatgatagtgaatggcaacgtgaccaatttttgtctcaagttaaaataccaaaaactattacagtgtctactggatttatgtctatatgacaaatcttgatactgcatccacaacattactggcgtgctttttgctttgcttttgtgtgcttttgtgtgtctgcctattaatacgtccgctgcttttgtctgtgtctacatacacatcattaataatattggtattactattgtggataacagcagcctctgcgtttaggtgttttattgtatatattatatttgtttatataccattatttttaatacatacacatgcacgctttttaattacataatgtatatgtacataatgtaattgttacatataattgttgtataccataacttactattttttcttttttattttcatatataatttttttttttgtttgtttgtttgttttttaataaactgttattacttaacaatgcgacacaaacgttctgcaaaacgcacaaaacgtgcatcggctacccaactttataaaacatgcaaacaggcaggtacatgtccacctgacattatacctaaggttgaaggcaaaactattgctgatcaaatattacaatatggaagtatgggtgtattttttggtgggttaggaattggaacagggtcgggtacaggcggacgcactgggtatattccattgggaacaaggcctcccacagctacagatacacttgctcctgtaagaccccctttaacagtagatcctgtgggcccttctgatccttctatagtttctttagtggaagaaactagttttattgatgctggtgcaccaacatctgtaccttccattcccccagatgtatcaggatttagtattactacttcaactgataccacacctgctatattagatattaataatactgttactactgttactacacataataatcccactttcactgacccatctgtattgcagcctccaacacctgcagaaactggagggcattttacactttcatcatccactattagtacacataattatgaagaaattcctatggatacatttattgttagcacaaaccctaacacagtaactagtagcacacccataccagggtctcgcccagtggcacgcctaggattatatagtcgcacaacacaacaggttaaagttgtagaccctgcttttgtaaccactcccactaaacttattacatatgataatcctgcatatgaaggtatagatgtggataatacattatatttttctagtaatgataatagtattaatatagctccagatcctgactttttggatatagttgctttacataggccagcattaacctctaggcgtactggcattaggtacagtagaattggtaataaacaaacactacgtactcgtagtggaaaatctataggtgctaaggtacattattattatgatttaagtactattgatcctgcagaagaaatagaattacaaactataacaccttctacatatactaccacttcacatgcagcctcacctacttctattaataatggattatatgatatttatgcagatgactttattacagatacttctacaaccccggtaccatctgtaccctctacatctttatcaggttatattcctgcaaatacaacaattccttttggtggtgcatacaatattcctttagtatcaggtcctgatatacccattaatataactgaccaagctccttcattaattcctatagttccagggtctccacaatatacaattattgctgatgcaggtgacttttatttacatcctagttattacatgttacgaaaacgacgtaaacgtttaccatattttttttcagatgtctctttggctgcctagtgaggccactgtctacttgcctcctgtcccagtatctaaggttgtaagcacggatgaatatgttgcacgcacaaacatatattatcatgcaggaacatccagactacttgcagttggacatccctattttcctattaaaaaacctaacaataacaaaatattagttcctaaagtatcaggattacaatacagggtatttagaatacatttacctgaccccaataagtttggttttcctgacacctcattttataatccagatacacagcggctggtttgggcctgtgtaggtgttgaggtaggtcgtggtcagccattaggtgtgggcattagtggccatcctttattaaataaattggatgacacagaaaatgctagtgcttatgcagcaaatgcaggtgtggataatagagaatgtatatctatggattacaaacaaacacaattgtgtttaattggttgcaaaccacctataggggaacactggggcaaaggatccccatgtaccaatgttgcagtaaatccaggtgattgtccaccattagagttaataaacacagttattcaggatggtgatatggttgatactggctttggtgctatggactttactacattacaggctaacaaaagtgaagttccactggatatttgtacatctatttgcaaatatccagattatattaaaatggtgtcagaaccatatggcgacagcttatttttttatttacgaagggaacaaatgtttgttagacatttatttaatagggctggtactgttggtgaaaatgtaccagacgatttatacattaaaggctctgggtctactgcaaatttagccagttcaaattattttcctacacctagtggttctatggttacctctgatgcccaaatattcaataaaccttattggttacaacgagcacagggccacaataatggcatttgttggggtaaccaactatttgttactgttgttgatactacacgcagtacaaatatgtcattatgtgctgccatatctacttcagaaactacatataaaaatactaactttaaggagtacctacgacatggggaggaatatgatttacagtttatttttcaactgtgcaaaataaccttaactgcagacgttatgacatacatacattctatgaattccactattttggaggactggaattttggtctacaacctcccccaggaggcacactagaagatacttataggtttgtaacatcccaggcaattgcttgtcaaaaacatacacctccagcacctaaagaagatccccttaaaaaatacactttttgggaagtaaatttaaaggaaaagttttctgcagacctagatcagtttcctttaggacgcaaatttttactacaagcaggattgaaggccaaaccaaaatttacattaggaaaacgaaaagctacacccaccacctcatctacctctacaactgctaaacgcaaaaaacgtaagctgtaa'
    }
    args = {
        'out_dir': './puma_out_test',
        'data_dir': 'data_dir',
        'input_format': 'genbank',
        'min_prot_len': 25,
        'e_value': 1e-05,
        'for_user_dir': './puma_out_test/for_user',
        'program_files_dir': './puma_out_test/program_files'
    }
    found = {
        'D21208':
        'MFQDTDEKPRNLHELCEALETTVHEISLPCVQCKKTLDRNEVYDFLFTDLKIVYRCGNPYGVCKQCLRLLSKVSEYRYFNYSVYGNTLEEIVHKPLNEITIRCITCQRPLCPQEKQRHVDRKKRFHNISNRWTGRCSVCWRPQRTQTQV*',
        'D90400':
        'MFQDAEEKPRTLHDLCQALETSVHEIELKCVECKKTLQRSEVYDFVFADLRIVYRDGNPFAVCKVCLRLLSKISEYRHYNYSLYGDTLEQTLKKCLNEILIRCIICQRPLCPQEKKRHVDLNKRFHNISGRWTGRCAVCWRPRRRQTQV*',
        'DQ080079':
        'MALFHNPEERPYKLPDLCRTLDTTLHDVTIDCVYCRRQLQRTEVYEFAFSDLCVVYRDGVPFAACQSCIKFYAKIRELRYYSESVYATTLETITNTKLYNLLIRCMSCLKPLCPAEKLRHLTTKRRLHKIAGNFTGQCRHCWTSKREDRRRIRQETQV*',
        'J04353':
        'MFKNPAERPRKLHELSSALEIPYDELRLNCVYCKGQLTETEVLDFAFTDLTIVYRDDTPHGVCTKCLRFYSKVSEFRWYRYSVYGTTLEKLTNKGICDLLIRCITCQRPLCPEEKQRHLDKKKRFHNIGGRWTGRCIACWRRPRTETQV*',
        'K02718':
        'MFQDPQERPRKLPQLCTELQTTIHDIILECVYCKQQLLRREVYDFAFRDLCIVYRDGNPYAVCDKCLKFYSKISEYRHYCYSLYGTTLEQQYNKPLCDLLIRCINCQKPLCPEEKQRHLDKKQRFHNIRGRWTGRCMSCCRSSRTRRETQL*',
        'M12732':
        'MFQDTEEKPRTLHDLCQALETTIHNIELQCVECKKPLQRSEVYDFAFADLTVVYREGNPFGICKLCLRFLSKISEYRHYNYSVYGNTLEQTVKKPLNEILIRCIICQRPLCPQEKKRHVDLNKRFHNISGRWAGRCAACWRSRRRETAL*',
        'NC_001587': '',
        'X74477': '',
        'X74481': '',
        'X94165': ''
    }
    res = puma.parse_blast_results_verify_e6(virus, args)
    assert isinstance(res, dict)
    #assert res == found


# --------------------------------------------------
def test_align_verify_e6():
    """Docstring"""

    aligned = puma.align_verify_e6()
    expected = '\n'.join([
        'SingleLetterAlphabet() alignment with 11 rows and 587 columns',
        '--------------------------------------------...--- KU350625',
        '--------------------------------------------...--- NC_015268',
        '--------------------------------------------...TGA NC_011051',
        '--------------------------------------------...TGA KX954132',
        '--------------------------------------------...--- JX174437',
        '--------------------------------------------...--- X02346',
        '--------------------------------------------...--- unknown',
        '--------------------------------------------...--- M20219',
        'ATGATCACACCATCACCGTTTTTTCAAGCGGGAAAAAAAAAGAC...--- JQ798171',
        '------------------------------ATGGGAATCTCTGG...TGA KP276343',
        '--------------------------------------------...--- U83595',
    ])

    assert isinstance(aligned, Bio.Align.MultipleSeqAlignment)
    assert str(aligned) == expected


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
