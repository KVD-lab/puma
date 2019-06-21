"""
puma library

authors: Josh Pace, Ken Youens-Clark, Cordell Freeman, Koenraad Van Doorslaer
University of Arizona, KVD Lab & Hurwitz Lab
PuMA 0.4 New L1, E5 variant and verify E6 additions 6/19/2019
"""

from distutils.spawn import find_executable
from Bio import SeqIO, GenBank, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import MafftCommandline
import os, glob, re, csv, time, operator, argparse, sys
import warnings
from Bio import BiopythonWarning
import itertools
import shutil
import logging
from dna_features_viewer import GraphicFeature, GraphicRecord
from subprocess import getstatusoutput
import numpy as np
import pandas as pd
from pprint import pprint as pp
from dire import die


class NoBLASTHits(Exception):
    pass


# ----------------------------------------------------------------------------------------
def make_l1_end(l1Result, original_genome, original_genome_length):
    """
    This function makes the L1 stop postion the last nucleotide of the genome and L1
    stop +1 the start of the genome

    :param l1Result: protein information on L1
    :param original_genome: nonlinearized genome
    :param original_genome_length: nonlinearized genome length
    :return: linearized genome based on L1
    """
    start = l1Result[0]
    stop = l1Result[1]

    #When gene wraps around
    if start < original_genome_length and stop > original_genome_length:
        sequence = original_genome[original_genome_length + 1:stop]
        new_genome = original_genome + sequence  #Still has extra at the beginning
        new_start = len(sequence)
        new_genome = new_genome[new_start + 1:]
    elif start > original_genome_length and stop > original_genome_length:  #No wrap around
        stop = stop - original_genome_length
        sequence = original_genome[stop + 1:]
        new_genome = sequence + original_genome  # Still has extra at the beginning
        new_stop = stop + len(sequence)
        new_genome = new_genome[:new_stop + 1]
    elif start < original_genome_length and stop < original_genome_length:  #No wrap around
        sequence = original_genome[stop:]
        new_stop = stop + len(sequence)
        new_genome = sequence + original_genome  # Still has extra at the beginning
        new_genome = new_genome[:new_stop]
    elif stop == original_genome_length:
        new_genome = original_genome

    return new_genome


# ----------------------------------------------------------------------------------------
def linearize_genome(original_genome, args):
    """
    This function checks for possible reverse complement issue and uses make_l1_end() to
    linearize genome based on L1

    :param original_genome: user inputted PV genome
    :param args: command line arguments, used for evalue, data_dir, out_dir, min_prot_len
    :return: linearized based on L1 and verified genome
    """

    orignal_genome_len = len(original_genome)
    extended_genome = original_genome + original_genome[:2000]

    if str(extended_genome).upper().count("N") > 1:
        logging.warning("Number of n nucleotides found:{}".format(number_of_n))
        raise Exception("Genome has too many n nucleotides to be annotated")

    proteins = identify_main_proteins(extended_genome, args)
    l1_result_extended = proteins['L1']
    # Determining if the reverse complement of the genome is needed and if it is a virus
    #  at all
    if l1_result_extended[0] == 'RC?':
        reverse_complement_genome = Seq(extended_genome).complement()
        l1_result_extended = identify_main_proteins(reverse_complement_genome,
                                                    args)
        if l1_result_extended[0] != 'RC?':
            extended_genome = reverse_complement_genome
        elif l1_result_extended[0] == 'RC?':
            raise Exception("PuMA cannot find an L1 protein")
    altered_genome = make_l1_end(l1_result_extended, original_genome,
                                 orignal_genome_len)

    altered_genome_len = len(altered_genome)
    if altered_genome_len != orignal_genome_len:
        raise Exception(
            "Genome Lengths do not match, there was a problem linearizing "
            "genome")

    return altered_genome


# ----------------------------------------------------------------------------------------
def trans_orf(seq, min_protein_length):
    """
    This functions finds all the open reading frames and translates them

    :param seq: PV genome
    :param min_protein_length: default is 25 or can be changed by user, ORFs shorter
    then value will be discarded
    :return: a dictionary where the key is the ORF and the value is the start of the
    ORF in the genome
    """

    # 1 == mammal or summat
    trans_table = 1

    ORFs = {}
    has_m = re.compile('M')
    for frame in range(3):
        trans = str(seq[frame:].translate(trans_table))
        trans_len = len(trans)
        aa_start = 0
        aa_end = 0
        while aa_start < trans_len:
            aa_end = trans.find("*", aa_start)
            if aa_end == -1:
                aa_end = trans_len
            if aa_end - aa_start >= min_protein_length:
                start = frame + aa_start * 3
                aa_seq = trans[aa_start:aa_end]
                if has_m.search(aa_seq):
                    ORFs[aa_seq] = start
            aa_start = aa_end + 1

    return ORFs


# ----------------------------------------------------------------------------------------
def run_blastp(query, subject, outfile, evalue=1e-4):
    """
    This function runs blastp

    :param query: file that contains query sequences
    :param subject: file that contains the sequences used for blasting
    :param outfile: the path to the file that will contain blast results if the blastp
    is successful
    :param evalue: threshold value for blast results, default is 1e-4 or user defined
    :return: returns the amount of blast hits which says if the blastp was successful
    """
    cmd = blastp(query=query,
                 subject=subject,
                 evalue=evalue,
                 outfmt=6,
                 out=outfile)

    #print(cmd)
    stdout, stderr = cmd()

    #if stderr:
    #    logging.warn("STDERR = ", stderr)

    if not os.path.isfile(outfile):
        raise Exception('No BLAST output')

    if os.path.getsize(outfile) == 0:
        return 0

    return len(open(outfile).read().splitlines())


# ----------------------------------------------------------------------------------------
def blast_main_orfs(genome, args):
    """
    This functions writes the ORFs found from trans_orf() to a fasta file and calls
    run_blastp

    :param genome: linearized based on L1 genome
    :param args: command line arguments for data_dir, min_prot_len etc.
    :return: a list of two file paths, the first element being the blast output and the
    second being the ORFs found
    """

    file_results = []
    orfs = trans_orf(genome, args['min_prot_len'])

    if not orfs:
        raise Exception('No ORFs, must stop.')

    orfs_fa = os.path.join(args['out_dir'], 'orfs.fa')
    orfs_fh = open(orfs_fa, 'wt')

    for orf in orfs:
        orfs_fh.write('\n'.join(['>' + str(orfs[orf]), orf, '']))
    orfs_fh.close()

    blast_sub = os.path.join(args['data_dir'], 'main_blast.fa')
    blast_out = os.path.join(args['out_dir'], 'blast_results_main.tab')

    if os.path.isfile(blast_out):
        os.remove(blast_out)
    try:
        num_hits = run_blastp(orfs_fa, blast_sub, blast_out)
        #print('Got {} hits!'.format(num_hits))
    except Exception as e:
        die('EXCEPTION: {}'.format(e))

    #file_results.append(blast_out)
    #file_results.append(orfs_fa)
    #return file_results

    return blast_out, orfs_fa


# ----------------------------------------------------------------------------------------
def identify_main_proteins(genome, args):
    """
    This function uses the blast output and ORF files from blast_main_orfs() to
    identify proteins except L1 and E5 variants in the ORFs

    :param genome: linearized based on L1 genome
    :param args: command line arguments for data_dir etc.
    :return: returns a dictionary (found_proteins) where the key is the name of the
    protein and the value is a list of start and stop positions, nucleotide seq and
    translated seq
    """
    protein_start = {}
    protein_seq = {}
    found_proteins = {}

    blast_out, orfs_fa = blast_main_orfs(genome, args)

    if not os.path.isfile(blast_out) or not os.path.getsize(blast_out):
        found_proteins['L1'] = ['RC?']  # Possible RC?
        return found_proteins

    with open(blast_out) as tab_file:
        for line in csv.reader(tab_file, delimiter="\t"):
            protein_start[line[1]] = int(line[0])

    with open(orfs_fa) as fasta_file:
        for line in fasta_file:
            for num in protein_start:
                try:
                    start = int(line[1:])
                    if start == protein_start[num]:
                        seq = next(fasta_file)
                        protein_seq[num] = seq[:-1]
                except:
                    pass

    for seq in protein_seq:
        for start in protein_start:
            if seq == start:
                if seq == 'L1':
                    found_proteins.update(
                        identify_l1(genome, protein_seq[seq],
                                    protein_start[start]))
                else:
                    try:
                        M = re.search('M', protein_seq[seq])
                        real_start = protein_start[start] + M.start(
                        ) + M.start() + M.start()
                        end = protein_start[start] + (
                            (len(protein_seq[seq]) + 1) * 3)
                        sequence = str(
                            genome[int(real_start):int(end)]).lower()
                        translated = Seq(sequence).translate()
                        found_proteins[seq] = [
                            int(real_start) + 1,
                            int(end), sequence, translated
                        ]
                    except AttributeError:
                        pass

    return found_proteins


# ----------------------------------------------------------------------------------------
def identify_l1(genome, sequence, start):
    """
    This function identifies the L1 from the blast results and ORFs

    :param genome: linearized based on L1 genome
    :param sequence: ORF identfied as having a potential L1
    :param start: ORF start position in genome
    :return: dictionary with the key being L1 and the value is a list of start and stop
    positions, nucleotide seq and translated seq
    """
    found_proteins = {}
    M = re.search('M', sequence)
    real_start = start + M.start() + M.start() + M.start()
    end = start + ((len(sequence) + 1) * 3)
    found_proteins['L1'] = []
    L1_pre = genome[(start + 3 * M.start()):int(end)]
    splice = '(C|T)(C|T)(A|C|G|T)(C|T)AG(A)TG'
    spliced = re.search(splice, str(L1_pre))
    if spliced:
        start_L1 = int(spliced.start()) + 6
        if start_L1 % 3 == 0:
            if start_L1 > 600:
                L1_post = L1_pre
                found_proteins['L1'] = [
                    int(start_L1),
                    int(end),
                    str(L1_post).lower(),
                    Seq(str(L1_post)).translate()
                ]
            else:
                L1_post = L1_pre[start_L1:]
                found_proteins['L1'] = [
                    int(real_start) + 1 + int(start_L1),
                    int(end),
                    str(L1_post).lower(),
                    Seq(str(L1_post)).translate()
                ]
        else:
            L1_post = L1_pre
            found_proteins['L1'] = [
                int(real_start) + 1,
                int(end),
                str(L1_post).lower(),
                Seq(str(L1_post)).translate()
            ]
    else:
        L1_post = L1_pre
        found_proteins['L1'] = [
            int(real_start) + 1,
            int(end),
            str(L1_post).lower(),
            Seq(str(L1_post)).translate()
        ]
    return found_proteins


# ----------------------------------------------------------------------------------------
def blast_e5_variants(virus, args):
    """
    This function uses a sequence 200 nucleotides at the end of E2 to 200 nucleotides
    into L2 find ORFs to be blasted to potentially locate E5 variants

    :param virus: Dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: a list of two file paths, the first element being the blast output and the
    second being the ORFs found
    """
    file_results = []
    data_dir = args['data_dir']
    out_dir = args['out_dir']
    e5_sequence = virus['genome'][virus['E2'][1] - 200:virus['L2'][0] + 200]
    orfs = trans_orf(Seq(e5_sequence), 5)
    orfs_fa = os.path.join(out_dir, 'orfs_E5.fa')
    orfs_fh = open(orfs_fa, 'wt')
    for orf in orfs:
        orfs_fh.write('\n'.join(['>' + str(orfs[orf]), orf, '']))
    orfs_fh.close()
    blast_sub = os.path.join(data_dir, 'blast_E5.fa')
    blast_out = os.path.join(out_dir, 'blast_results_E5.tab')
    if os.path.isfile(blast_out):
        os.remove(blast_out)
    try:
        num_hits = run_blastp(orfs_fa, blast_sub, blast_out, evalue=2)
        #print('Got {} hits!'.format(num_hits))
    except Exception as e:
        die('EXCEPTION: {}'.format(e))

    file_results.append(blast_out)
    file_results.append(orfs_fa)

    return file_results


# ----------------------------------------------------------------------------------------
def identify_e5_variants(virus, args):
    """
    This function searches the linearized genome for possible E5_ALPHA, BETA, DELTA,
    EPSILON, GAMA, ZETA

    :param virus: Dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: dictionary, either empty of with found E5 variants
    """
    protein_seq = {}
    found_proteins = {}
    file_results = blast_e5_variants(virus, args)
    blast_out = file_results[0]
    orfs_fa = file_results[1]
    hdrs = ('qseqid sseqid pident length mismatch gapopen qstart '
            'qend sstart send evalue bitscore').split()
    df = pd.read_csv(blast_out, sep='\t', names=hdrs)
    evalue = 1e-1
    wanted = df[df.apply(lambda x: '_' in x['sseqid'] and x['evalue'] < evalue,
                         axis=1)]
    if len(wanted) < 1:
        return found_proteins
    protein_start = dict(zip(wanted['sseqid'], wanted['qseqid']))
    orfs_dict = {}
    with open(orfs_fa) as fasta_file:
        for line in fasta_file:
            if line[0] == '>':
                start = int(line[1:])
                orf_seq = next(fasta_file)
                orf_seq = str(orf_seq).replace("\n", "")
                orfs_dict[start] = orf_seq
    for start in protein_start:
        if protein_start[start] in orfs_dict:
            protein_seq[protein_start[start]] = orfs_dict[protein_start[start]]
    for seq in protein_seq:
        for start in protein_start:
            if seq == protein_start[start]:
                try:
                    M = re.search('M', protein_seq[seq])
                    real_start = protein_start[start] + M.start() + M.start(
                    ) + M.start()
                    end = protein_start[start] + (
                        (len(protein_seq[seq]) + 1) * 3)
                    sequence = str(
                        e5_sequence[int(real_start):int(end)]).lower()
                    translated = Seq(sequence).translate()
                    genome_start = re.search(sequence, virus['genome']).start()
                    genome_end = re.search(sequence, virus['genome']).end()
                    found_proteins[start] = [
                        int(genome_start) + 1,
                        int(genome_end), sequence, translated
                    ]
                except AttributeError:
                    pass
    for e5 in list(found_proteins):
        if len(found_proteins[e5][3]) < 40:
            del found_proteins[e5]
    return found_proteins


# ----------------------------------------------------------------------------------------
def blast_verify_e6(virus, args):
    """
    This function blasts the found E6 agaisnt all E6s in PaVE

    :param virus: Dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: file path to blast output
    """
    E6_trans = str(virus['E6'][3])
    out_dir = args['out_dir']
    data_dir = args['data_dir']
    verify_E6_dir = os.path.join(out_dir, 'verify_E6')
    ID = virus['accession']

    if not os.path.isdir(verify_E6_dir):
        os.makedirs(verify_E6_dir)
    blast_subject = os.path.join(data_dir, 'blast_E6_updated.fa')
    blast_out = os.path.join(verify_E6_dir, 'blast_result_E6.tab')
    if os.path.isfile(blast_out):
        os.remove(blast_out)
    query_file = os.path.join(verify_E6_dir, 'query.fa')
    with open(query_file, 'a') as query:
        query.write('>{}\n'.format(ID))
        query.write(E6_trans)
    try:
        num_hits = run_blastp(query_file, blast_subject, blast_out)
        # print('Got {} hits!'.format(num_hits))
    except Exception as e:
        die('EXCEPTION: {}'.format(e))
    return blast_out


# ----------------------------------------------------------------------------------------
def parse_blast_results_verify_e6(virus, args):
    """
    This functions parses the blast output from blast_verify_E6() for
    results that have an evalue less then or equal to 1e-43

    :param virus: Dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: dictionary, keys are the accession number and the values are the E6
    sequence of identified sequences that fall at or below the evalue of 1e-43
    """
    blast_out = blast_verify_e6(virus, args)
    blast_options = []
    last_resort_blast_options = []
    number_over_eval = 0
    count = 0
    e_values = {}
    data_dir = args['data_dir']

    with open(blast_out) as blast_file:
        blast_result = csv.reader(blast_file, delimiter='\t')
        for row in blast_result:
            count = count + 1
            e_values[row[1]] = float(row[-2])
            if count > 10:
                break
    for genome in e_values:
        last_resort_blast_options.append(genome)
        if (e_values[genome] > float(1e-43)):
            number_over_eval = number_over_eval + 1
        else:
            blast_options.append(genome)
    if number_over_eval > 0:
        print("Blast results for verifying E6 fall below "
              "the set confidence level. "
              "The results are still being reported. "
              "Number found below the confidence level is:{}.".format(
                  number_over_eval))
    if len(blast_options) == 0:
        blast_options = last_resort_blast_options[0:10]
    else:
        blast_options = blast_options[0:10]
    logging.debug("E6 blast options:{}".format(blast_options))
    known_E6 = {}  # Stores a dictionary of the 10 closest blast results

    csv_database = os.path.join(data_dir, 'all_pave_new_updated.csv')

    with open(csv_database, 'r') as csvfile:
        read = csv.DictReader(csvfile,
                              fieldnames=[
                                  'accession', 'gene', 'positions', 'seq',
                                  'translated seq'
                              ])
        for row in read:
            if row['accession'] in blast_options and row['gene'] == 'E6':
                known_E6[row['accession']] = str(row['translated seq'])
    logging.debug(known_E6)
    return known_E6


# ----------------------------------------------------------------------------------------
def align_verify_e6(virus, args):
    """
    This function uses MUSCLE to align the found E6 and the E6s from the blast results

    :param virus: Dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: Biopython alignment object
    """

    E6_trans = str(virus['E6'][3])
    out_dir = args['out_dir']
    known_E6 = parse_blast_results_verify_e6(virus, args)
    verify_E6_dir = os.path.join(out_dir, 'verify_E6')
    unaligned = os.path.join(verify_E6_dir, 'unaligned.fa')
    aligned = os.path.join(verify_E6_dir, 'aligned.fa')

    if os.path.isfile(unaligned):
        os.remove(unaligned)
    if os.path.isfile(aligned):
        os.remove(aligned)
    logging.debug("E6 before:{}".format(E6_trans))
    with open(unaligned, 'a') as sequence_file:
        sequence_file.write(">{}\n".format('unknown'))
        sequence_file.write("{}\n".format(E6_trans))

    for key in known_E6:
        with open(unaligned, 'a') as sequence_file:
            sequence_file.write(">{}\n".format(key))
            sequence_file.write("{}\n".format(known_E6[key]))

    if find_executable('muscle'):
        cline = MuscleCommandline(input=unaligned, out=aligned, verbose=False)
        stdout, stderr = cline()
    else:
        raise Exception('muscle not installed')

    alignment = AlignIO.read(aligned, 'fasta')
    logging.debug(alignment)

    return alignment


# ----------------------------------------------------------------------------------------
def verify_e6(virus, args):
    """
    This function use the alignment object from align_verify_e6() to identify if the
    found E6 is potentially too long at the beginning of the sequence

    :param virus: Dictionary that has all found proteins so far and linearized based on L1
    genome
    :return: dictionary, key is E6 and the value is a list of start and stop
    positions, nucleotide seq and translated seq of the updated E6
    """

    verified_E6 = {}
    dashes = {}
    conserved = {}
    alignment = align_verify_e6(virus, args)
    alignment_length = alignment.get_alignment_length()
    E6_seq = str(virus['E6'][2])
    E6_whole = virus['E6']
    ID = virus['accession']

    seq_by_id = dict([(rec.id, str(rec.seq)) for rec in alignment])
    for i, rec in enumerate(alignment):
        seq = str(rec.seq)
        match = re.search(r'^([-]+)', seq)
        num_dashes = 0
        if match:
            dash = match.group(1)
            num_dashes = len(dash)

        dashes[rec.id] = num_dashes

    for i in range(0, alignment_length):
        col = alignment[:, i]
        conserved[i] = col.count('M')

    max_num_met = max(conserved.values())
    seqs_at_max = list(filter(lambda t: t[1] == max_num_met,
                              conserved.items()))
    num_of_dashes = seqs_at_max[0][0]
    prefix = seq_by_id['unknown'][0:num_of_dashes]

    if len(prefix) * '-' == prefix:
        actual_start = 0
    elif prefix.count('-') == 0:
        #print('E6 = ', seq_by_id['unknown'][num_of_dashes:])
        actual_start = len(prefix) * 3
    else:
        no_dashes = len(prefix.replace('-', ''))
        actual_start = no_dashes * 3

    test_seq = E6_seq[actual_start:]
    trans = Seq(str(test_seq)).translate()
    logging.debug("test trans:{}".format(trans))
    if str(Seq(str(test_seq)).translate())[0] != 'M':
        test_trans = str(Seq(str(test_seq[3:])).translate())
        if test_trans[0] == 'M':
            actual_start = actual_start + 3
        else:
            actual_start = 0
    if len(trans) <= 110:
        actual_start = 0
    if len(trans) >= 184:
        logging.debug("SEQ IS TOO LONG")
    new_seq = E6_seq[actual_start:]

    verified_E6[ID] = [
        E6_whole[0] + actual_start, E6_whole[1],
        str(new_seq).lower(),
        Seq(str(new_seq)).translate()
    ]

    return verified_E6


# ----------------------------------------------------------------------------------------
def find_urr(virus):
    """
    This function locates the URR in the linearized genome based on the postions of the
    proteins found and the definition of the URR

    :param virus: Dictionary that has all found proteins so far and linearized based on L1
    genome
    :return: dictionary, where the key is URR and and the value list of start and stop
    positions and nucleotide seq
    """
    start_stop = []
    urr = {}
    blasted = {}
    blasted.update(virus)
    altered_genome = blasted['genome']
    del blasted['name']
    del blasted['accession']
    del blasted['genome']

    altered_genome_len = len(altered_genome)
    for protein in blasted:
        if protein == 'L1':
            start_stop.append((blasted[protein][1]))
            urr_start = blasted[protein][1]
        else:
            start_stop.append(blasted[protein][0])

    start_stop = sorted(start_stop)
    #print("start_stop:{}".format(start_stop))

    for numbers in start_stop:
        if numbers == urr_start:
            if numbers == start_stop[-1]:
                urr_stop = start_stop[0]
            else:
                position = start_stop.index(numbers)
                urr_stop = start_stop[position + 1]
    urr_start = int(urr_start)
    urr_stop = int(urr_stop) - 1

    if urr_start == altered_genome_len:
        urr_start = 1
    if urr_stop == 0:
        urr_stop = altered_genome_len
    if urr_stop > urr_start:
        urr_found = str(altered_genome[urr_start - 1:urr_stop]).lower()
        urr['URR'] = [int(urr_start), int(urr_stop), urr_found]
    else:
        urr_found = str(altered_genome[urr_start - 1:] +
                        altered_genome[:urr_stop]).lower()
        urr['URR'] = [
            int(urr_start),
            int(genomelen), 1,
            int(urr_stop), urr_found
        ]

    return urr


# ----------------------------------------------------------------------------------------
def find_E1BS(genome, URR, URRstart, ID, data_dir, out_dir):
    """
    Finds the E1 binding site in the genome using the URR

    :param genome:
    :param URR:
    :param URRstart:
    :param ID:
    :param data_dir:
    :param out_dir:
    :return:
    """

    genomeLength = len(genome)
    E1BS = {}  # Storing E1BS
    startListURR = []
    startURR = 0

    tmp = os.path.join(out_dir, "puma_urr.fa")
    with open(tmp, "w") as tempfile:
        tempfile.write('>URR for {}\n'.format(ID))
        tempfile.write(str(URR))

    fimo_exe = find_executable('fimo')
    if not fimo_exe:
        print('Please install "fimo" into your $PATH')
        return

    fimo_dir = os.path.join(out_dir, 'E1BS')
    if not os.path.isdir(fimo_dir):
        os.makedirs(fimo_dir)

    background = os.path.join(data_dir, 'background_model_E1BS.txt')
    motif = os.path.join(data_dir, 'E1BS_motif.txt')

    fimo_cmd = '{} --oc {} --norc --verbosity 1 --thresh 1.0E-1 --bgfile {} {} {}'
    fimo_out = os.path.join(fimo_dir, 'fimo.tsv')

    if os.path.isfile(fimo_out):
        os.remove(fimo_out)

    cline = (fimo_cmd.format(fimo_exe, fimo_dir, background, motif, tmp))

    rv, out = getstatusoutput(str(cline))
    if rv != 0:
        raise Exception('Failed to run fimo for E1BS')

    if not os.path.isfile(fimo_out):
        logging.warning('Failed to create fimo out "{}"'.format(fimo_out))
        return

    with open(fimo_out, "rU") as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            start = row['start']
            if not start is None:
                startListURR.append(start)

    startURR = int(startListURR[0])
    genomestart = (startURR + URRstart)

    if genomestart > genomeLength:
        genomestart = genomestart - genomeLength

    genomestop = genomestart + 19

    if genomestop > genomeLength:
        genomestop = genomestop - genomeLength
        sequence = str(genome[int(genomestart) - 2:] +
                       genome[:genomestop]).lower()
        E1BS['E1BS'] = [
            int(genomestart) - 1,
            int(genomeLength), 1,
            int(genomestop), sequence
        ]

    else:
        if genomestart == 1:
            sequence = str(genome[-1]).lower() + str(
                genome[int(genomestart - 1):int(genomestop)]).lower()
            E1BS['E1BS'] = [int(genomestart), int(genomestop), sequence]

        else:
            sequence = str(genome[int(genomestart -
                                      2):int(genomestop)]).lower()
            E1BS['E1BS'] = [int(genomestart) - 1, int(genomestop), sequence]
    return E1BS


# ----------------------------------------------------------------------------------------
def find_E2BS(genome, URR, URRstart, ID, data_dir, out_dir):
    """
    This function finds the E2BS in a genome using the URR
    :param genome:
    :param URR:
    :param URRstart:
    :param ID:
    :param data_dir:
    :param out_dir:
    :return:
    """
    genomeLength = len(genome)
    startListURR = []
    startListGenome = []
    E2BS = {}

    tmp = os.path.join(out_dir, "puma_urr.fa")
    with open(tmp, "w") as tempfile:
        tempfile.write('>URR for {}\n'.format(ID))
        tempfile.write(str(URR))

    fimo_exe = find_executable('fimo')
    if not fimo_exe:
        raise Exception('Cannot find fimo in your $PATH')

    fimo_dir = os.path.join(out_dir, 'E2BS')
    if not os.path.isdir(fimo_dir):
        os.makedirs(fimo_dir)

    motif = os.path.join(data_dir, 'E2BS_motif.txt')

    fimo_out = os.path.join(fimo_dir, 'fimo.tsv')

    if os.path.isfile(fimo_out):
        os.remove(fimo_out)

    fimo_cmd = '{} --oc {} --norc --verbosity 1 --thresh 1.0E-4 {} {}'
    cline = (fimo_cmd.format(fimo_exe, fimo_dir, motif, tmp))

    # os.system(str(cline))

    rv, out = getstatusoutput(str(cline))
    if rv != 0:
        raise Exception('Failed to run fimo for E2BS: {}'.format(out))

    if not os.path.isfile(fimo_out):
        logging.critical('Failed to create fimo out "{}"'.format(fimo_out))
        return

    with open(fimo_out, "rU") as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        try:
            for row in reader:
                start = row['start']
                if not start is None:
                    startListURR.append(row['start'])

            startListURR = list(map(int, set(startListURR)))

            startListGenome = list(map(int, startListGenome))

            for i in range(0, len(startListURR), 1):
                genomestart = startListURR[i]
                genomestart = (genomestart + URRstart)
                if genomestart > genomeLength:
                    genomestart = genomestart - genomeLength
                    startListGenome.append(genomestart - 1)
                else:
                    startListGenome.append(genomestart - 1)

            E2BS['E2BS'] = startListGenome
            #print("E2BS:{}".format(E2BS))
            return E2BS

        except KeyError:

            E2BS['E2BS'] = ["No E2BS found"]
            return E2BS


# ----------------------------------------------------------------------------------------
def find_E4(E2, genome, position):
    """
    Finds the E4, used for stop site for E1^E4

    :param E2:
    :param genome:
    :param position:
    :return:
    """
    E4 = {}
    E4_orfs = []
    trans_E2 = Seq(E2[1:len(E2)]).translate()
    E4_orfs = str(trans_E2).split("*")
    E4_orfs.sort(key=len)
    E4protein_long = E4_orfs[position]
    if 'M' in E4protein_long:
        M = re.search('M', E4protein_long)
        if M.start() > 41:
            E4protein = E4protein_long
        else:
            E4protein_tmp = E4protein_long.split('M', 1)[1]
            E4protein = 'M' + str(E4protein_tmp)
    else:
        E4protein = E4protein_long

    E4_start = re.search(str(E4protein), str(trans_E2)).start()
    E4_end = re.search(str(E4protein), str(trans_E2)).end()
    E4_nt = str(E2[(E4_start * 3) + 1:((E4_end + 1) * 3) + 1])
    E4_nt_start = re.search(E4_nt, str(genome)).start()
    E4_nt_end = E4_nt_start + len(E4_nt)
    sequence = str(genome[int(E4_nt_start):int(E4_nt_end)]).lower()
    translated = Seq(sequence).translate()
    E4['E4'] = [int(E4_nt_start) + 1, int(E4_nt_end) + 1]
    return E4


# ----------------------------------------------------------------------------------------
def find_blastResults_E1E8(E1_whole, ID, out_dir, data_dir):
    """
    Finds the blast results for E1 and E8

    :param E1_whole:
    :param ID:
    :param out_dir:
    :param data_dir:
    :return:
    """
    E1_trans = str(E1_whole[3])

    blastE1E8_dir = os.path.join(out_dir, 'blastE1E8')
    if not os.path.isdir(blastE1E8_dir):
        os.makedirs(blastE1E8_dir)

    blast_subject = os.path.join(data_dir, 'E1E8_blast.fa')
    blast_out = os.path.join(blastE1E8_dir, 'blast_result.tab')

    if os.path.isfile(blast_out):
        os.remove(blast_out)

    query_file = os.path.join(blastE1E8_dir, 'query.fa')

    with open(query_file, 'a') as query:
        query.write('>{}\n'.format(ID))
        query.write(E1_trans)

    try:
        num_hits = run_blastp(query_file, blast_subject, blast_out)
        #print('Got {} hits!'.format(num_hits))
        blastResult = True
    except Exception:
        blastResult = False
        return blastResult

    return blastResult


# ----------------------------------------------------------------------------------------
def find_E1E4(E1_whole, E2_whole, ID, genome, start_E4_nt, blastE1E8_dir,
              blastResult, data_dir):
    """
   Finds E1^E4

   :param E1_whole:
   :param E2_whole:
   :param ID:
   :param genome:
   :param start_E4_nt:
   :param blastE1E8_dir:
   :param blastResult:
   :param data_dir:
   :return:
   """
    E1_E4 = {}
    genome = str(genome).lower()
    E2_seq = str(E2_whole[2])
    E1_seq = str(E1_whole[2])

    if start_E4_nt == False:
        E1_E4['E1^E4'] = False
        return E1_E4

    if blastResult:
        blast_out = os.path.join(blastE1E8_dir, 'blast_result.tab')
    else:
        E1_E4['E1^E4'] = [0, 0, 0, 0, 'No Blast Output', '']
        return E1_E4

    blast_options = []
    with open(blast_out) as blast_file:
        blast_result = csv.reader(blast_file, delimiter='\t')
        for row in blast_result:
            blast_options.append(row[1])
    blast_options = blast_options[0:1]
    splice_sites = []
    E1_stop = []
    for options in blast_options:
        query = options
        known_E1 = {}
        csv_database = os.path.join(data_dir, 'all_pave_new_updated.csv')

        with open(csv_database, 'r') as csvfile:
            read = csv.DictReader(csvfile,
                                  ('accession', 'gene', 'positions', 'seq'))
            for row in read:
                if row['accession'] == query and row['gene'] == 'E1^E4':
                    E1_positions = row['positions']
                if row["accession"] == query and row["gene"] == 'E1':
                    known_E1[query] = str(row['seq']).lower()
                    known_E1_stop = str(row["positions"])
                if row['accession'] == query and row['gene'] == 'CG':
                    known_CG = str(row['seq']).lower()

        try:
            E1_stop_genome = E1_positions.split('+')[0]
            E1_stop_genome = E1_stop_genome.split('(')[1]
            # print(E8_stop_genome)

            E1_stop_genome = int(str(E1_stop_genome).split('..')[1])

            known_E1_stop = int(known_E1_stop.split('..')[0])

            E1_stop_known = (E1_stop_genome - known_E1_stop)

            splice_sites.append(E1_stop_known)

        except UnboundLocalError:
            print('Query does not have E1^E4 or E8^E2')
            stopE1_nt = False
            return stopE1_nt

        unaligned = os.path.join(blastE1E8_dir, 'unaligned.fa')
        aligned = os.path.join(blastE1E8_dir, 'aligned.fa')

        if os.path.isfile(unaligned):
            os.remove(unaligned)

        if os.path.isfile(aligned):
            os.remove(aligned)

        for key in known_E1:
            with open(unaligned, 'a') as sequence_file:
                sequence_file.write(">{}\n".format(ID))
                sequence_file.write("{}\n".format(E1_seq))
                sequence_file.write(">{}\n".format(key))
                sequence_file.write("{}\n".format(known_E1[key]))

        cline = MuscleCommandline(input=unaligned, out=aligned, verbose=False)

        stdout, stderr = cline()

        align_seq = []
        for aln in AlignIO.read(aligned, 'fasta'):
            align_seq.append(aln.seq)

        unknown_seq = str(align_seq[0]).lower()
        known_seq = str(align_seq[1]).lower()

        j = 0
        aligned_E1_stop = 0

        for position in known_seq:
            aligned_E1_stop = aligned_E1_stop + 1
            # print(aligned_E8_stop)
            if position.lower() in ['a', 'c', 't', 'g']:
                j = j + 1
                if j == E1_stop_known:
                    break

        E1_stop.append(aligned_E1_stop)

        #print("E1_Stop Options:{}".format(E1_stop))

    aligned_E1_stop = E1_stop[-1]
    search_seq = unknown_seq[aligned_E1_stop:aligned_E1_stop + 50].replace(
        '-', '')

    stopE1_nt = (re.search(search_seq, str(genome).lower()).start()) + 1
    startE1_nt = E1_whole[0]

    whole_E4 = find_E4(E2_seq, genome, -1)
    #For when it may not always be the longest nonexistent E4, checks to make sure that
    #  the splice acceptor site is within the found E4
    # logging.debug("E4 Start:{}".format(whole_E4['E4'][0]))
    # logging.debug("E4 Stop:{}".format(whole_E4['E4'][1]))
    # logging.debug("Splice Site Staart:{}".format(start_E4_nt))
    position = -1
    while ((start_E4_nt < whole_E4['E4'][0])
           and (start_E4_nt < whole_E4['E4'][1])):
        whole_E4 = find_E4(E2_seq, genome, position)
        position = position - 1

    stopE4_nt = whole_E4['E4'][1] - 1

    E1_E4_seq = str(genome[startE1_nt - 1:stopE1_nt] +
                    genome[start_E4_nt:stopE4_nt])
    E1_E4_trans = Seq(E1_E4_seq).translate()

    #Fixing "*" issue
    E1_part = Seq(genome[startE1_nt - 1:stopE1_nt]).translate()

    E1_E4['E1^E4'] = [
        startE1_nt, stopE1_nt, start_E4_nt + 1, stopE4_nt, E1_E4_seq,
        E1_E4_trans
    ]
    if "*" in E1_E4_trans[:-1]:
        logging.debug("E1 part:{}".format(E1_part))
        logging.debug("E1_E4 Translated:{}".format(E1_E4_trans))
        #print("E1:{}".format(E1_whole))
    return E1_E4


# ----------------------------------------------------------------------------------------
def find_E8E2(E1_whole, E2_whole, ID, genome, startE2_nt, blastE1E8_dir,
              blastResult, data_dir):
    """
    Finds E8^E2

    :param E1_whole:
    :param E2_whole:
    :param ID:
    :param genome:
    :param startE2_nt:
    :param blastE1E8_dir:
    :param blastResult:
    :param data_dir:
    :return:
    """

    E8_E2 = {}
    E1_seq = str(E1_whole[2])
    E1_trans = str(E1_whole[3])
    genome = str(genome).lower()
    startE8List = []
    startMotif = 'MKL'

    if startE2_nt == False:

        E8_E2['E8^E2'] = False
        return E8_E2

    if blastResult:
        blast_out = os.path.join(blastE1E8_dir, 'blast_result.tab')
    else:
        E8_E2['E8^E2'] = [0, 0, 0, 0, 'No Blast Output', '']
        return E8_E2

    blast_options = []
    with open(blast_out) as blast_file:
        blast_result = csv.reader(blast_file, delimiter='\t')
        for row in blast_result:
            blast_options.append(row[1])
    blast_options = blast_options[0:1]
    splice_sites = []
    E8_stop = []
    for options in blast_options:
        query = options
        known_E8 = {}
        csv_database = os.path.join(data_dir, 'all_pave_new_updated.csv')

        with open(csv_database, 'r') as csvfile:
            read = csv.DictReader(csvfile,
                                  ('accession', 'gene', 'positions', 'seq'))
            for row in read:
                if row['accession'] == query and row['gene'] == 'E8^E2':
                    E8_positions = row['positions']
                if row["accession"] == query and row["gene"] == 'E1':
                    known_E8[query] = str(row['seq']).lower()
                    known_E8_stop = str(row["positions"])
                if row['accession'] == query and row['gene'] == 'CG':
                    known_CG = str(row['seq']).lower()

        try:
            E8_stop_genome = E8_positions.split('+')[0]
            E8_stop_genome = E8_stop_genome.split('(')[1]
            E8_stop_genome = int(str(E8_stop_genome).split('..')[1])
            known_E8_stop = int(known_E8_stop.split('..')[0])
            E8_stop_known = (E8_stop_genome - known_E8_stop)
            splice_sites.append(E8_stop_known)
        except UnboundLocalError:
            E8_E2['E8^E2'] = False
            return E8_E2

        unaligned = os.path.join(blastE1E8_dir, 'unaligned.fa')
        aligned = os.path.join(blastE1E8_dir, 'aligned.fa')
        if os.path.isfile(unaligned):
            os.remove(unaligned)

        if os.path.isfile(aligned):
            os.remove(aligned)

        for key in known_E8:
            with open(unaligned, 'a') as sequence_file:
                sequence_file.write(">{}\n".format(ID))
                sequence_file.write("{}\n".format(E1_seq))
                sequence_file.write(">{}\n".format(key))
                sequence_file.write("{}\n".format(known_E8[key]))

        cline = MuscleCommandline(input=unaligned, out=aligned, verbose=False)

        stdout, stderr = cline()

        align_seq = []
        for aln in AlignIO.read(aligned, 'fasta'):
            align_seq.append(aln.seq)

        unknown_seq = str(align_seq[0]).lower()
        known_seq = str(align_seq[1]).lower()

        j = 0
        aligned_E8_stop = 0

        for position in known_seq:
            aligned_E8_stop = aligned_E8_stop + 1
            #print(aligned_E8_stop)
            if position.lower() in ['a', 'c', 't', 'g']:
                j = j + 1
                if j == E8_stop_known:
                    break

        E8_stop.append(aligned_E8_stop)

    aligned_E8_stop = E8_stop[-1]
    #logging.debug("Aligned_E8_Stops:{}".format(aligned_E8_stop))
    search_seq = unknown_seq[aligned_E8_stop:aligned_E8_stop + 50].replace(
        '-', '')

    stopE8_nt = (re.search(search_seq, str(genome).lower()).start()) + 1

    stopE8 = (stopE8_nt - E1_whole[0]) - 1
    tempStart = stopE8 - 70

    search_seq = E1_seq[tempStart:stopE8]

    for match in re.finditer('atg', search_seq):
        startE8List.append(match.start() + tempStart)
    for i in range(0, len(startE8List)):
        try:
            if (stopE8 - startE8List[i]) < 15:
                del startE8List[i]
            elif (stopE8 - startE8List[i]) > 65:
                del startE8List[i]
        except IndexError:
            break
    #logging.debug("StartlistE8:{}".format(startE8List))
    for i in range(0, len(startE8List)):
        try:
            testStart = startE8List[i] + E1_whole[0]
            check_seq = Seq(genome[testStart - 1:stopE8_nt]).translate()
            if "*" in check_seq:
                del startE8List[i]
            elif check_seq.startswith(startMotif):
                startE8_nt = startE8List[i] + E1_whole[0]
        except IndexError:
            break
    #logging.debug("StartlistE8:{}".format(startE8List))
    if len(startE8List) == 0:
        E8_E2['E8^E2'] = ["No E8 found"]
        return E8_E2

    if check_seq.startswith(startMotif):
        pass
    else:
        try:
            startE8_nt = startE8List[0] + E1_whole[0]
        except IndexError:
            E8_E2['E8^E2'] = ["No E8 found"]
            return E8_E2

    stopE8_nt = (stopE8 + E1_whole[0]) + 1

    if startE2_nt != False:
        stopE2_nt = E2_whole[1]

    E8_E2_seq = Seq(genome[startE8_nt - 1:stopE8_nt] +
                    genome[startE2_nt:stopE2_nt])
    E8_E2_trans = E8_E2_seq.translate()
    E8_part = Seq(genome[startE8_nt - 1:stopE8_nt]).translate()

    E8_E2['E8^E2'] = [
        startE8_nt, stopE8_nt, startE2_nt + 1, stopE2_nt, E8_E2_seq,
        E8_E2_trans
    ]

    if "*" in E8_E2_trans[:-1]:
        logging.debug(E8_E2)
        logging.debug(E2_whole)
        logging.debug("E8 part:{}".format(E8_part))
        logging.debug("E8_E2 Translated:{}".format(E8_E2_trans))
    return E8_E2


# ----------------------------------------------------------------------------------------
def find_splice_acceptor(E2_whole, ID, genome, data_dir, out_dir):
    """
    Finds splice acceptor posistion for E1^E4 and E8^E2

    :param E2_whole:
    :param ID:
    :param genome:
    :param data_dir:
    :param out_dir:
    :return:
    """

    E2_trans = str(E2_whole[3])
    E2_seq = E2_whole[2]

    splice_acceptor_dir = os.path.join(out_dir, 'splice_acceptor')
    if not os.path.isdir(splice_acceptor_dir):
        os.makedirs(splice_acceptor_dir)

    blast_subject = os.path.join(data_dir, 'splice_acceptor_blast.fa')
    blast_out = os.path.join(splice_acceptor_dir,
                             'blast_result_splice_acceptor.tab')

    if os.path.isfile(blast_out):
        os.remove(blast_out)

    query_file = os.path.join(splice_acceptor_dir, 'query.fa')

    with open(query_file, 'a') as query:
        query.write('>{}\n'.format(ID))
        query.write(E2_trans)

    try:
        num_hits = run_blastp(query_file, blast_subject, blast_out)
        #print('Got {} hits!'.format(num_hits))
    except Exception as e:
        startE2_nt = "No blast results for unkown E2"
        return startE2_nt
        #die('EXCEPTION: {}'.format(e))

    blast_options = []
    with open(blast_out) as blast_file:
        blast_result = csv.reader(blast_file, delimiter='\t')
        for row in blast_result:
            blast_options.append(row[1])
    blast_options = blast_options[0:1]
    splice_sites = []
    aligned_starts = []
    for options in blast_options:
        query = options
        known_E2 = {}
        csv_database = os.path.join(data_dir, 'all_pave_new_updated.csv')

        with open(csv_database, 'r') as csvfile:
            read = csv.DictReader(csvfile,
                                  ('accession', 'gene', 'positions', 'seq'))
            for row in read:
                if row['accession'] == query and row['gene'] == 'E8^E2':
                    splice_acceptor_positions = row['positions']
                if row["accession"] == query and row["gene"] == 'E2':
                    known_E2[query] = str(row['seq']).lower()
                    known_E2_start = str(row["positions"])
                if row['accession'] == query and row['gene'] == 'CG':
                    known_CG = str(row['seq']).lower()

        try:
            splice_start_genome = splice_acceptor_positions.split('+')[1]

            splice_start_genome = int(str(splice_start_genome).split('..')[0])

            known_E2_start = int(known_E2_start.split('..')[0])

            splice_start_known = (splice_start_genome - known_E2_start)

            splice_sites.append(splice_start_known)

        except UnboundLocalError:
            print('Query does not have E1^E4 or E8^E2')
            startE2_nt = False
            return startE2_nt

        unaligned = os.path.join(splice_acceptor_dir, 'unaligned.fa')
        aligned = os.path.join(splice_acceptor_dir, 'aligned.fa')

        if os.path.isfile(unaligned):
            os.remove(unaligned)

        if os.path.isfile(aligned):
            os.remove(aligned)

        for key in known_E2:
            with open(unaligned, 'a') as sequence_file:
                sequence_file.write(">{}\n".format(ID))
                sequence_file.write("{}\n".format(E2_seq))
                sequence_file.write(">{}\n".format(key))
                sequence_file.write("{}\n".format(known_E2[key]))

        cline = MuscleCommandline(input=unaligned, out=aligned, verbose=False)

        stdout, stderr = cline()

        align_seq = []
        for aln in AlignIO.read(aligned, 'fasta'):
            align_seq.append(aln.seq)

        unknown_seq = str(align_seq[0]).lower()
        known_seq = str(align_seq[1]).lower()

        j = 0
        aligned_splice_start = 0

        for position in known_seq:
            aligned_splice_start = aligned_splice_start + 1
            if position.lower() in ['a', 'c', 't', 'g']:
                j = j + 1
                if j == splice_start_known:
                    break

        aligned_starts.append(aligned_splice_start)
        #print("Aligned Splice Starts:{}".format(aligned_starts))
    aligned_splice_start = aligned_starts[-1]

    search_seq = unknown_seq[aligned_splice_start:aligned_splice_start +
                             50].replace('-', '')

    startE2_nt = re.search(search_seq, str(genome).lower()).start()
    #makelogging.debug("Splice Site Start:{}".format(startE2_nt))
    return startE2_nt


# ----------------------------------------------------------------------------------------
def to_gff3(dict, genomelen, out_dir):
    """
    Outputs gff3 format

    :param dict:
    :param genomelen:
    :param out_dir:
    :return:
    """
    del dict['URR']

    name = dict['accession']
    dict['name'] = name
    gff3_out = os.path.join(out_dir, '{}.gff3'.format(dict['name']))

    with open(gff3_out, 'a') as out_file:
        out_file.write("##gff-version 3\n")
        out_file.write("##sequence-region {} 1 {}\n".format(
            dict['name'], genomelen))

    for protein in dict:
        if protein == 'name' or protein == 'accession':
            pass
        elif "^" in protein:
            frameNumber = dict[protein][0] % 3
            if frameNumber == 1:
                frame = 1
            elif frameNumber == 2:
                frame = 2
            elif frameNumber == 0:
                frame = 3
            with open(gff3_out, 'a') as out_file:
                out_file.write(
                    "{}\tPuMA\tsplice_site\t{}\t{}\t{}\t{}\t.\t+\t{}\tID={};"
                    "Note=[{}-{} + {}-{"
                    "}]\n".format(dict['name'], dict[protein][0],
                                  dict[protein][1], dict[protein][2],
                                  dict[protein][3], frame, protein,
                                  dict[protein][0], dict[protein][1],
                                  dict[protein][2], dict[protein][3]))
        else:
            frameNumber = dict[protein][0] % 3
            if frameNumber == 1:
                frame = 1
            elif frameNumber == 2:
                frame = 2
            elif frameNumber == 0:
                frame = 3
            with open(gff3_out, 'a') as out_file:
                out_file.write("{}\tPuMA\tCDS\t{}\t{}\t.\t+\t{}\tID={};"
                               "Note=[{}-{}]\n".format(dict['name'],
                                                       dict[protein][0],
                                                       dict[protein][1], frame,
                                                       protein,
                                                       dict[protein][0],
                                                       dict[protein][1]))

    return


# ----------------------------------------------------------------------------------------
def to_pdf(annotations, out_dir):
    """
    Genome visualization to pdf

    :param annotations:
    :param out_dir:
    :return:
    """

    sequence_length = len(annotations['genome'])
    try:
        del annotations['genome']
        del annotations['E1BS']
        del annotations['E2BS']
    except KeyError:
        pass

    features = []

    pdf_out = os.path.join(out_dir, '{}.pdf'.format(annotations['accession']))

    for value in annotations:
        if value == 'name' or value == 'accession':
            pass
        elif 'E8^' in value:
            features.append(
                GraphicFeature(start=annotations[value][0],
                               end=annotations[value][1],
                               strand=+1,
                               color="#800080",
                               label='E8'))
            features.append(
                GraphicFeature(start=annotations[value][2],
                               end=annotations[value][3],
                               strand=+1,
                               color="#800080",
                               label='E2'))
        elif 'E1^' in value:
            features.append(
                GraphicFeature(start=annotations[value][0],
                               end=annotations[value][1],
                               strand=+1,
                               color="#40e0d0",
                               label='E1'))
            features.append(
                GraphicFeature(start=annotations[value][2],
                               end=annotations[value][3],
                               strand=+1,
                               color="#40e0d0",
                               label='E4'))
        elif 'URR' in value:
            try:
                features.append(
                    GraphicFeature(start=annotations[value][0],
                                   end=annotations[value][1],
                                   strand=+1,
                                   color="#ffff00",
                                   label='URR'))
                features.append(
                    GraphicFeature(start=annotations[value][2],
                                   end=annotations[value][3],
                                   strand=+1,
                                   color="#ffff00",
                                   label='URR'))

            except IndexError:
                features.append(
                    GraphicFeature(start=annotations[value][0],
                                   end=annotations[value][1],
                                   strand=+1,
                                   color="#ffff00",
                                   label='URR'))
        elif 'L' in value:
            features.append(
                GraphicFeature(start=annotations[value][0],
                               end=annotations[value][1],
                               strand=+1,
                               color="#0000ff",
                               label=value))
        elif 'E' in value:
            features.append(
                GraphicFeature(start=annotations[value][0],
                               end=annotations[value][1],
                               strand=+1,
                               color="#ffa500",
                               label=value))

    record = GraphicRecord(sequence_length=sequence_length, features=features)
    ax, _ = record.plot(figure_width=10)
    ax.figure.savefig(pdf_out)
    return


# ----------------------------------------------------------------------------------------


def export_to_csv(annotations, out_dir):
    """
        Output to csv file

        :param annotations:
        :param out_dir:
        :return:
        """

    csv_out = os.path.join(out_dir, 'puma_results.csv')
    try:
        if annotations['E2BS'] == ['No E2BS found']:
            del annotations['E2BS']
    except KeyError:
        pass

    with open(csv_out, 'a') as out:
        out_file = csv.writer(out)

        for value in annotations:
            if value == 'genome':
                out_file.writerows([[
                    annotations['accession'], 'CG', "",
                    str(annotations[value]).lower()
                ]])
            elif value == 'URR':
                try:
                    out_file.writerows([[
                        annotations['accession'], value,
                        'join(' + str(annotations[value][0]) + '..' +
                        str(annotations[value][1]) + '+' +
                        str(annotations[value][2]) + ".." +
                        str(annotations[value][3]) + ')', annotations[value][4]
                    ]])
                except IndexError:
                    out_file.writerows([[
                        annotations['accession'], value,
                        str(annotations[value][0]) + '..' +
                        str(annotations[value][1]), annotations[value][2]
                    ]])
            elif value == 'E1BS':
                try:
                    out_file.writerows([[
                        annotations['accession'], value,
                        'join(' + str(annotations[value][0]) + '..' +
                        str(annotations[value][1]) + '+' +
                        str(annotations[value][2]) + ".." +
                        str(annotations[value][3]) + ')', annotations[value][4]
                    ]])
                except IndexError:
                    out_file.writerows([[
                        annotations['accession'], value,
                        str(annotations[value][0]) + '..' +
                        str(annotations[value][1]), annotations[value][2]
                    ]])
            elif value == 'E2BS':
                for i in range(0, len(annotations[value]), 1):
                    out_file.writerows([[
                        annotations['accession'], value,
                        str(annotations[value][i]) + '..' +
                        str(annotations[value][i] + 11),
                        str(annotations['genome'][annotations[value][i] -
                                                  1:annotations[value][i] +
                                                  11]).lower()
                    ]])
            elif value == 'E1^E4':
                out_file.writerows([[
                    annotations['accession'], value,
                    'join(' + str(annotations[value][0]) + '..' +
                    str(annotations[value][1]) + '+' +
                    str(annotations[value][2]) + ".." +
                    str(annotations[value][3]) + ')', annotations[value][4],
                    annotations[value][5]
                ]])
            elif value == 'E8^E2':
                out_file.writerows([[
                    annotations['accession'], value,
                    'join(' + str(annotations[value][0]) + '..' +
                    str(annotations[value][1]) + '+' +
                    str(annotations[value][2]) + ".." +
                    str(annotations[value][3]) + ')', annotations[value][4],
                    annotations[value][5]
                ]])
            elif value == 'name':
                pass
            elif value == 'accession':
                pass
            else:
                out_file.writerows([[
                    annotations['accession'], value,
                    str(annotations[value][0]) + '..' +
                    str(annotations[value][1]), annotations[value][2],
                    annotations[value][3]
                ]])

    return


# ----------------------------------------------------------------------------------------
def to_results(dict):
    """
    :param dict:
    :return:
    """
    try:
        del dict['genome']
        del dict['accession']
        del dict['E1BS']
        del dict['E2BS']
    except KeyError:
        pass
    all = dict['name']
    #short_name = re.search('\(([^)]+)', all).group(1)
    #dict['name'] = short_name

    results_dir = os.path.join('puma_results')

    results = os.path.join(results_dir, 'puma_results_6_19_19.fa')

    for protein in dict:
        if protein == 'name':
            pass
        elif protein == 'accession':
            pass
        elif protein == 'URR':
            try:
                if type(dict[protein][3]) == int:
                    with open(results, 'a') as out_file:
                        out_file.write(">{}, {}\n".format(
                            dict['name'], protein))
                        out_file.write("{}\n".format(dict[protein][4]))
                else:
                    with open(results, 'a') as out_file:
                        out_file.write(">{}, {}\n".format(
                            dict['name'], protein))
                        out_file.write("{}\n".format(dict[protein][2]))
            except IndexError:
                with open(results, 'a') as out_file:
                    out_file.write(">{}, {}\n".format(dict['name'], protein))
                    out_file.write("{}\n".format(dict[protein][2]))

        elif '^' in protein:
            with open(results, 'a') as out_file:
                out_file.write(">{}, {} gene\n".format(dict['name'], protein))
                out_file.write("{}\n".format(dict[protein][4]))
        else:
            with open(results, 'a') as out_file:
                out_file.write(">{}, {} gene\n".format(dict['name'], protein))
                out_file.write("{}\n".format(dict[protein][2]))

    return


# ----------------------------------------------------------------------------------------
def print_genome_info(args, virus):

    print(
        "\nThis is the gene information for {} after making L1 end of Genome:".
        format(virus['accession']))
    # if args['sites'] == ['ALL']:
    #     sites = {}
    #     sites.update(virus)
    #     del sites['name']
    #     del sites['genome']
    #     del sites['accession']
    for name in sites:
        if name == 'E2BS':
            print("\n{} E2 binding sites found:".format(len(virus['E2BS'])))
            for i in range(0, len(virus['E2BS'])):
                print('\n{} start and stop position:\n{},{}\n'.format(
                    name, virus[name][i], virus[name][i] + 11))
                print('{} sequence:\n{}\n'.format(
                    name,
                    str(virus['genome'][virus['E2BS'][i] - 1:virus['E2BS'][i] +
                                        11]).lower()))
        elif name == 'E1BS':
            if type(virus[name][2]) == int:
                print('\n{} start and stop position:\n{},{},{},{}\n'.format(
                    name, virus[name][0], virus[name][1], virus[name][2],
                    virus[name][3]))
                print('{} seqeunce:\n{}\n'.format(name, virus[name][4]))
            else:
                print('\n{} start and stop position:\n{},{}\n'.format(
                    name, virus[name][0], virus[name][1]))
                print('{} sequence:\n{}\n'.format(name, virus[name][2]))
        else:
            try:
                if type(virus[name][3]) == int:
                    print(
                        '\n{} start and stop position:\n{},{},{},{}\n'.format(
                            name, virus[name][0], virus[name][1],
                            virus[name][2], virus[name][3]))
                    print('{} seqeunce:\n{}\n'.format(name, virus[name][4]))
                    if name != 'URR':
                        print('{} translated sequnce:\n{}\n'.format(
                            name, virus[name][5][:-1]))
                else:
                    print('\n{} start and stop position:\n{},{}\n'.format(
                        name, virus[name][0], virus[name][1]))
                    print('{} sequence:\n{}\n'.format(name, virus[name][2]))
                    if name != 'URR':
                        print('{} translated seqeunce:\n{}\n'.format(
                            name, virus[name][3][:-1]))
            except IndexError:
                print('\n{} start and stop position:\n{},{}\n'.format(
                    name, virus[name][0], virus[name][1]))
                print('{} sequence:\n{}\n'.format(name, virus[name][2]))
                if name != 'URR':
                    print('{} translated seqeunce:\n{}\n'.format(
                        name, virus[name][3][:-1]))

    return


# ----------------------------------------------------------------------------------------
def validate_args(args):
    """
    Validate arguments

    Returns a `dict` of arguments
    """

    #sites = args['sites'] if 'sites' in args else ['ALL']
    #sites = args.get('sites') or []
    input_file = args.get('input')
    out_dir = args.get('outdir')
    data_dir = args.get('data_dir')
    input_format = args.get('format').lower()
    min_prot_len = args.get('min_prot_len')
    e_value = args.get('evalue')
    create_gff3 = args.get('gff3')
    create_csv = args.get('csv')
    blast_e1_e8_dir = os.path.join(out_dir, 'blastE1E8')
    if not os.path.isdir(blast_e1_e8_dir):
        os.makedirs(blast_e1_e8_dir)

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    return {
        'input_file': input_file,
        'out_dir': out_dir,
        'data_dir': data_dir,
        'input_format': input_format,
        'min_prot_len': min_prot_len,
        'e_value': e_value,
        'create_gff3': create_gff3,
        'create_csv': create_csv,
        'blast_e1_e8_dir': blast_e1_e8_dir
    }


# --------------------------------------------------
def run(args):
    """main"""
    virus = {}
    blasted = {}
    logging.info('run = {}\n'.format(args))
    warnings.simplefilter('ignore', BiopythonWarning)
    args = validate_args(args)
    for seq_record in SeqIO.parse(args['input_file'], args['input_format']):
        original_genome = seq_record.seq
        name = seq_record.description.split(",")[0]
        ID = seq_record.name
    #full_name = name.split("|")[1]
    #ID = name.split("|")[0]
    virus['name'] = name  # full_name
    virus['accession'] = ID
    print(ID)
    altered_genome = linearize_genome(original_genome, args)
    altered_genome_len = len(altered_genome)
    virus['genome'] = str(altered_genome).lower()
    blasted.update(identify_main_proteins(altered_genome, args))
    virus.update(blasted)
    #print(blasted)
    virus.update(identify_e5_variants(virus, args))

    if 'E6' in virus.keys():
        verified_E6 = verify_e6(virus, args)
        del virus['E6']
        virus['E6'] = verified_E6[ID]
        # print("E6 from genome:{}".format(virus['genome'][virus['E6'][0]-1:virus['E6'][1]]))
        # print("E6 seq:{}".format(virus['E6'][2]))
    virus.update(find_urr(virus))

    E2BS = find_E2BS(altered_genome, virus['URR'][-1], virus['URR'][0], ID,
                     args['data_dir'], args['out_dir'])
    if E2BS['E2BS'] != ["No E2BS found"]:
        virus.update(E2BS)

    E1BS = find_E1BS(altered_genome, virus['URR'][-1], virus['URR'][0], ID,
                     args['data_dir'], args['out_dir'])
    if E1BS:
        virus.update(E1BS)
    try:
        start_splice_site = find_splice_acceptor(virus['E2'], ID,
                                                 altered_genome,
                                                 args['data_dir'],
                                                 args['out_dir'])
    except KeyError:
        pass
    blast_result = find_blastResults_E1E8(virus['E1'], ID, args['out_dir'],
                                          args['data_dir'])
    try:
        E1_E4 = find_E1E4(virus['E1'], virus['E2'], ID, altered_genome,
                          start_splice_site, args['blast_e1_e8_dir'],
                          blast_result, args['data_dir'])
        E8_E2 = find_E8E2(virus['E1'], virus['E2'], ID, altered_genome,
                          start_splice_site, args['blast_e1_e8_dir'],
                          blast_result, args['data_dir'])
    except KeyError:
        pass
    try:
        if E8_E2['E8^E2'] == False:
            pass
        elif E8_E2['E8^E2'] == ["No E8 found"]:
            pass
        else:
            virus.update(E8_E2)
        if E1_E4['E1^E4'] == False:
            pass
        else:
            virus.update(E1_E4)
    except UnboundLocalError:
        pass
    print_genome_info(args, virus)

    if args['create_csv']:
        export_to_csv(virus, args['out_dir'])
    else:
        pass
    to_pdf(virus, args['out_dir'])
    if args['create_gff3']:
        to_gff3(virus, altered_genome_len, args['out_dir'])
    else:
        pass
    #to_results(virus)
    return 1
