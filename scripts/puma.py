"""
puma library

authors: Josh Pace, Ken Youens-Clark, Cordell Freeman, Koenraad Van Doorslaer
University of Arizona, KVD Lab & Hurwitz Lab

Refer to https://puma-docs.readthedocs.io/en/latest/?badge=latest for
documentation

github: https://github.com/KVD-lab/puma

Contact Dr. Koenraad Van Doorslaer at vandoorslaer@arizona.edu

PuMA 1.2 release 7/24/2020
"""


import os
import glob
import re
import csv
import time
import operator
import argparse
import sys
import random
import warnings
import itertools
import shutil
import logging
import numpy as np
import pandas as pd
from collections import Counter
from Bio import SeqIO
from Bio import GenBank
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import CompoundLocation
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio.Align.Applications import MuscleCommandline
from Bio import BiopythonWarning
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from subprocess import getstatusoutput
from distutils.spawn import find_executable
from pprint import pprint as pp

# --------------------------------------------------
def make_l1_end(l1Result, original_genome):
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
    new_genome = ''
    original_genome_length = len(original_genome)
    #When gene wraps around
    if start < original_genome_length and stop > original_genome_length:
        sequence = original_genome[original_genome_length + 1:stop]
        new_genome = original_genome + sequence  #Still has extra at the beginning
        new_start = len(sequence)
        new_genome = new_genome[new_start + 1:]
    elif start > original_genome_length and stop > original_genome_length: #No wrap around
        stop = stop - original_genome_length
        sequence = original_genome[stop + 1:]
        new_genome = sequence + original_genome  # Still has extra at the beginning
        new_stop = stop + len(sequence)
        new_genome = new_genome[:new_stop + 1]
    elif start < original_genome_length and stop < original_genome_length: #No wrap around
        sequence = original_genome[stop:]
        new_stop = stop + len(sequence)
        new_genome = sequence + original_genome  # Still has extra at the beginning
        new_genome = new_genome[:new_stop]
    # else?
    elif stop == original_genome_length:
        new_genome = original_genome
    return new_genome
# --------------------------------------------------
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
    other_nts = len(re.findall('[^ATCG]', str(original_genome), re.IGNORECASE))
    if other_nts > 0:
        logging.warning("Number of non ACTG nucleotides found:{}".format(other_nts))
    if other_nts >= 5:
        logging.critical("Genome has too many n nucleotides to be annotated")
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
    altered_genome = make_l1_end(l1_result_extended, original_genome)
    altered_genome_len = len(altered_genome)
    if altered_genome_len != orignal_genome_len:
        raise Exception(
            "Genome Lengths do not match, there was a problem linearizing "
            "genome")

    return altered_genome
# --------------------------------------------------
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
# --------------------------------------------------
def run_blastp(query, subject, outfile, evalue=1e-5):
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
    stdout, stderr = cmd()

    if not os.path.isfile(outfile):
        raise Exception('No BLAST output')

    if os.path.getsize(outfile) == 0:
        return 0

    return len(open(outfile).read().splitlines())
# --------------------------------------------------
def blast_main_orfs(genome, args):
    """
    This functions writes the ORFs found from trans_orf() to a fasta file and calls
    run_blastp

    :param genome: linearized based on L1 genome
    :param args: command line arguments for data_dir, min_prot_len etc.
    :return: a list of two file paths, the first element being the blast output and the
    second being the ORFs found
    """
    orfs = trans_orf(genome, args['min_prot_len'])
    if not orfs:
        raise Exception('No ORFs, must stop.')
    out_dir = args['program_files_dir']
    main_blast_dir = os.path.join(out_dir, 'main_blast')
    if not os.path.isdir(main_blast_dir):
        os.makedirs(main_blast_dir)
    orfs_fa = os.path.join(main_blast_dir, 'orfs.fa')
    orfs_fh = open(orfs_fa, 'wt')

    for orf in orfs:
        orfs_fh.write('\n'.join(['>' + str(orfs[orf]), orf, '']))
    orfs_fh.close()

    blast_sub = os.path.join(args['data_dir'], 'main_blast.fa')
    blast_out = os.path.join(main_blast_dir, 'blast_results_main.tab')

    if os.path.isfile(blast_out):
        os.remove(blast_out)
    num_hits = run_blastp(orfs_fa, blast_sub, blast_out)
    return blast_out, orfs_fa
# --------------------------------------------------
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
                try:
                    M = re.search('M', protein_seq[seq])
                    real_start = protein_start[start] + M.start(
                    ) + M.start() + M.start()
                    end = protein_start[start] + (
                        (len(protein_seq[seq]) + 1) * 3)
                    sequence = str(
                        genome[int(real_start):int(end)]).lower()
                    translated = Seq(sequence).translate()
                    if len(translated) < 25:
                        logging.info("{} is of length {} and therefore is too short to be "
                                  "considered viable".format(seq, len(translated)))
                    else:
                        found_proteins[seq] = [
                            int(real_start) + 1,
                            int(end), sequence, translated
                        ]
                except AttributeError:
                    pass

    return found_proteins
# --------------------------------------------------
def blast_verify_gene(gene, accession, name, args):

    """
    This function blasts the current gene being verified against corresponding
    gene in blast_subject_all.fa
    :param gene: amino acid sequence of unknown gene
    :param accession: accession number of unknown gene
    :param name: gene being ran (L1, E1 etc.)
    :param args: command line arguments for data_dir etc
    :return: file path to blast output
    """

    out_dir = args['program_files_dir']
    data_dir = args['data_dir']
    verify_gene_dir = os.path.join(out_dir, 'verify_{}'.format(name))
    if not os.path.isdir(verify_gene_dir):
        os.makedirs(verify_gene_dir)
    blast_out = os.path.join(verify_gene_dir, 'blast_result_{}.tab'.format(name))
    if os.path.isfile(blast_out):
        os.remove(blast_out)
    query_file = os.path.join(verify_gene_dir, 'query.fa')
    with open(query_file, 'a') as query:
        query.write('>{}\n'.format(accession))
        query.write(gene)
    blast_subject = os.path.join(data_dir, 'blast_subject_all.fa')
    num_hits = run_blastp(query_file, blast_subject, blast_out)
    return blast_out
# --------------------------------------------------
def parse_blast_results_verify_gene(gene, accession, name, args):

    """
    This functions parses the blast output from blast_verify_gene() for the top
    10 best evalue hits

    :param gene: amino acid sequence of unknown gene
    :param accession: accession number of unknown gene
    :param name: gene being ran (L1, E1 etc.)
    :param args: command line arguments for data_dir etc
    :return: dictionary, keys are the accession number and the values are the
    amino acid sequence of blast result sequences
    """
    blast_out = blast_verify_gene(gene, accession, name, args)
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
        if (e_values[genome] > float(1e-34)):
            number_over_eval = number_over_eval + 1
        else:
            blast_options.append(genome)

    for genome in e_values:
            blast_options.append(genome)
    if number_over_eval > 0:
        logging.info("Blast results for verifying {} fall below "
                     "the set confidence level. "
                     "Number found below the confidence level is:{}.".format(
            name,number_over_eval))
    if len(blast_options) == 0:
        blast_options = last_resort_blast_options[0:10]
    else:
        blast_options = blast_options[0:10]

    known_gene = {}  # Stores a dictionary of the 10 closest blast results
    csv_database = os.path.join(data_dir, 'PaVE.csv')
    with open(csv_database, 'r') as csvfile:
        read = csv.DictReader(csvfile,
                              fieldnames=[
                                  'accession', 'gene', 'positions', 'seq',
                                  'translated seq'])
        for row in read:
            if row['accession'] in blast_options and row['gene'] == name:
                known_gene[row['accession']] = str(Seq(row['seq']).translate())[:-1]

    return known_gene
# --------------------------------------------------
def align_verify_gene(gene, accession, name, args):
    """
        This function uses MUSCLE to align the unkown gene and the genes from
        the blast results

        :param gene: amino acid sequence of unknown gene
        :param accession: accession number of unknown gene
        :param name: gene being ran (L1, E1 etc.)
        :param args: command line arguments for data_dir etc
        :return: Biopython alignment object
        """

    out_dir = args['program_files_dir']
    known_gene = parse_blast_results_verify_gene(gene, accession, name, args)
    verify_gene_dir = os.path.join(out_dir, 'verify_{}'.format(name))
    unaligned = os.path.join(verify_gene_dir, 'unaligned.fa')
    aligned = os.path.join(verify_gene_dir, 'aligned.fa')

    if os.path.isfile(unaligned):
        os.remove(unaligned)
    if os.path.isfile(aligned):
        os.remove(aligned)
    with open(unaligned, 'a') as sequence_file:
        sequence_file.write(">{}\n".format('unknown'))
        sequence_file.write("{}\n".format(gene))

    for key in known_gene:
        with open(unaligned, 'a') as sequence_file:
            sequence_file.write(">{}\n".format(key))
            sequence_file.write("{}\n".format(known_gene[key]))

    if find_executable('muscle'):
        cline = MuscleCommandline(input=unaligned, out=aligned, verbose=False)
        stdout, stderr = cline()
    else:
        raise Exception('muscle not installed')

    alignment = AlignIO.read(aligned, 'fasta')
    return alignment
# --------------------------------------------------
def verify_gene(virus, name, args):
    """
    This function use the alignment object from align_verify_gene() to identify
    if the unknown gene is potentially too long at the beginning of the sequence.
    This is done by looking at the alignment and seeing is a M is a match across
    all the the sequences at a certain location.

    What needs to be worked on is at what point should we not take the new start
    position if it's to far downstream. This is where potentially a seperate
    function with known parameters for each gene could be implemented.

    :param virus: dictionary that has all found proteins so far and linearized
    based on L1 genome
    :param name: gene being ran (L1, E1 etc.)
    :param args: command line arguments for data_dir etc
    :return: dictionary, key is unknown gene being currently ran  and the value
    is a list of start and stop positions, nucleotide seq and translated seq of
    the updated gene


    This function needs to be verified more to ensure that it works for all
    genes being verified. Could possibly need a check function for each gene.
    For example still use all the above related 'verifiy_' functions (that blast
    ,align etc.) for all but then have a specfic verify_e1() etc. for the final
    check
    """
    verified_gene = {}
    dashes = {}
    conserved = {}
    gene_trans = str(virus[name][3][:-1])
    gene_seq = str(virus[name][2])
    gene_all = virus[name]
    verified_gene[name] = [
            gene_all[0], gene_all[1],
            str(gene_seq).lower(),
            Seq(str(gene_seq)).translate()]
    alignment = align_verify_gene(gene_trans, virus['accession'], name, args)
    lengths = []
    max_num_met = ''
    for a in alignment:
        if a.id != 'unknown':
            lengths.append(len(str(a.seq).replace("-","")))
            
        else:
            unknown_len=len(str(a.seq).replace("-",""))
            unknown_seq=a.seq
    if len(lengths) != 0:
        average_len = (sum(lengths)/len(lengths))
        stdev_len =  np.std(lengths)
        if stdev_len < 2:
            stdev_len = 2
    else:
        average_len = 0
        stdev_len = 2
        
    alignment_length = alignment.get_alignment_length()

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
    x=0
    for c in conserved.values():
        if c > 2:
            if str(unknown_seq[x]) == 'M':
                max_num_met=c
                break
        x=x+1
    if  max_num_met != '':
        seqs_at_max = list(filter(lambda t: t[1] == max_num_met,
                                  conserved.items()))
        num_of_dashes = seqs_at_max[0][0]
        prefix = seq_by_id['unknown'][0:num_of_dashes]
        if len(prefix) * '-' == prefix:
            actual_start = 0
        elif prefix.count('-') == 0:
            actual_start = len(prefix) * 3
        else:
            no_dashes = len(prefix.replace('-', ''))
            actual_start = no_dashes * 3

        corrected_len = unknown_len - (actual_start/3)
        low_len=average_len-(3*stdev_len)
        hi_len=average_len+(3*stdev_len)
        if actual_start > 12:
            if low_len > corrected_len:
                actual_start = 0

        test_seq = gene_seq[actual_start:]
        trans = Seq(str(test_seq)).translate()
        if str(Seq(str(test_seq)).translate())[0] != 'M':
            test_trans = str(Seq(str(test_seq[3:])).translate())
            if test_trans[0] == 'M':
                actual_start = actual_start + 3
            else:
                actual_start = 0

        if len(trans) <= low_len:
            actual_start = 0
        new_seq = gene_seq[actual_start:]
    
        verified_gene[name] = [
            gene_all[0] + actual_start, gene_all[1],
            str(new_seq).lower(),
            Seq(str(new_seq)).translate()]

    return verified_gene
# --------------------------------------------------
def blast_e5_variants(virus, args):
    """
    This function uses a sequence 200 nucleotides at the end of E2 to 200 nucleotides
    into L2 find ORFs to be blasted to potentially locate E5 variants

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: a list of two file paths, the first element being the blast output and the
    second being the ORFs found
    """
    data_dir = args['data_dir']
    out_dir = args['program_files_dir']
    e5_variants_dir = os.path.join(out_dir, 'E5_variants')
    if not os.path.isdir(e5_variants_dir):
        os.makedirs(e5_variants_dir)
    e5_sequence = virus['genome'][virus['E2'][1] - 200:virus['L2'][0] + 200]
    orfs = trans_orf(Seq(e5_sequence), 5)
    orfs_fa = os.path.join(e5_variants_dir, 'orfs_E5.fa')
    orfs_fh = open(orfs_fa, 'wt')
    for orf in orfs:
        orfs_fh.write('\n'.join(['>' + str(orfs[orf]), orf, '']))
    orfs_fh.close()
    blast_sub = os.path.join(data_dir, 'blast_E5.fa')
    blast_out = os.path.join(e5_variants_dir, 'blast_results_E5.tab')
    if os.path.isfile(blast_out):
        os.remove(blast_out)

    num_hits = run_blastp(orfs_fa, blast_sub, blast_out, evalue=2)
    return blast_out, orfs_fa, e5_sequence
# --------------------------------------------------
def identify_e5_variants(virus, args):
    """
    This function searches the linearized genome for possible E5_ALPHA, BETA, DELTA,
    EPSILON, GAMA, ZETA

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: dictionary, either empty of with found E5 variants
    """
    protein_seq = {}
    found_proteins = {}
    blast_out , orfs_fa, e5_sequence = blast_e5_variants(virus,args)
    hdrs = ('qseqid sseqid pident length mismatch gapopen qstart '
            'qend sstart send evalue bitscore').split()
    df = pd.read_csv(blast_out, sep='\t', names=hdrs)
    evalue = 1e-1
    wanted = df[df.apply(lambda x: '_' in x['sseqid'] and x['evalue'] < evalue,
                         axis=1)]
    c = Counter(wanted['qseqid'])
    for qseqid, count in c.items():
        if count > 1:
            best = sorted(wanted[wanted['qseqid'] == qseqid]['evalue'])[0]
            indexes = wanted.index[(wanted['qseqid'] == qseqid)
                                   & (wanted['evalue'] > best)].tolist()
            wanted = wanted.drop(indexes)
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
# --------------------------------------------------
def find_urr(virus):
    """
    This function locates the URR in the linearized genome based on the postions of the
    proteins found and the definition of the URR

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :return: dictionary, where the key is URR and and the value is list of start and stop
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
            int(urr_stop), urr_found]
    return urr
# --------------------------------------------------
def fimo_e1bs(virus, args):
    """
    This function executes fimo on the URR for the E1 Binding Sites

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: TSV fimo output file
    """
    URR = virus['URR'][-1]
    ID = virus['accession']
    out_dir = args['program_files_dir']
    data_dir = args['data_dir']

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

    background = os.path.join(data_dir, 'background_model_E1BS_new.txt')
    motif = os.path.join(data_dir, 'E1BS_motif_new.txt')

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
    return fimo_out
# --------------------------------------------------
def find_e1bs(virus, args):
    """
    This function uses the fimo output from fimo_e1bs() along with the URR to find the
    E1 Binding Sites

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: dictionary, where the key is E1BS and and the value is list of start and
    stop
    positions and nucleotide seq
    """
    genome = virus['genome']
    genome_length = len(genome)
    E1BS = {}  # Storing E1BS
    start_list_urr = []
    fimo_out = fimo_e1bs(virus, args)

    with open(fimo_out, "r") as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            start = row['start']
            if not start is None:
                start_list_urr.append(start)
    start_urr = int(start_list_urr[0])
    genome_start = (start_urr + virus['URR'][0])
    if genome_start > genome_length:
        genome_start = genome_start - genome_length
    genome_stop = genome_start + 19
    if genome_stop > genome_length:
        genome_stop = genome_stop - genome_length
        sequence = str(genome[int(genome_start) - 1:] +
                       genome[:genome_stop]).lower()
        E1BS['E1BS'] = [
            int(genome_start),
            int(genome_length), 1,
            int(genome_stop), sequence]
    else:
        if genome_start == 1:
            sequence = str(genome[-1]).lower() + str(
                genome[int(genome_start - 1):int(genome_stop)]).lower()
            E1BS['E1BS'] = [int(genome_start), int(genome_stop), sequence]
        else:
            sequence = str(genome[int(genome_start -
                                      1):int(genome_stop)]).lower()
            E1BS['E1BS'] = [int(genome_start), int(genome_stop), sequence]
    return E1BS
# --------------------------------------------------
def fimo_e2bs(virus, args):
    """
    This function executes fimo on the URR for the E2 Binding Sites

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: TSV fimo output file
    """
    URR = virus['URR'][-1]
    ID = virus['accession']
    data_dir = args['data_dir']
    out_dir = args['program_files_dir']
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

    motif = os.path.join(data_dir, 'E2BS_motif_new.txt')

    fimo_out = os.path.join(fimo_dir, 'fimo.tsv')

    if os.path.isfile(fimo_out):
        os.remove(fimo_out)

    fimo_cmd = '{} --oc {} --norc --verbosity 1 --thresh 1.0E-4 {} {}'
    cline = (fimo_cmd.format(fimo_exe, fimo_dir, motif, tmp))

    rv, out = getstatusoutput(str(cline))
    if rv != 0:
        raise Exception('Failed to run fimo for E2BS: {}'.format(out))

    if not os.path.isfile(fimo_out):
        logging.critical('Failed to create fimo out "{}"'.format(fimo_out))
        return
    return fimo_out
# --------------------------------------------------
def find_e2bs(virus, args):
    """
    This function uses the fimo output from fimo_e2bs() along with the URR to find the
    E2 Binding Sites

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: dictionary, where the key is E2BS and and the value is list of start
    positions
    """

    genome_length = len(virus['genome'])
    start_list_urr = []
    start_list_genome = []
    E2BS = {}
    fimo_out = fimo_e2bs(virus, args)

    with open(fimo_out, "r") as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        try:
            for row in reader:
                start = row['start']
                if not start is None:
                    start_list_urr.append(row['start'])

            start_list_urr = list(map(int, set(start_list_urr)))

            start_list_genome = list(map(int, start_list_genome))

            for i in range(0, len(start_list_urr), 1):
                genome_start = start_list_urr[i]
                genome_start = (genome_start + virus['URR'][0])
                if genome_start > genome_length:
                    genome_start = genome_start - genome_length
                    start_list_genome.append(genome_start - 1)
                else:
                    start_list_genome.append(genome_start - 1)

            E2BS['E2BS'] = start_list_genome
        except KeyError:
            E2BS['E2BS'] = []
    return E2BS
# --------------------------------------------------
def blast_splice_acceptor(virus, args):
    """
    This function blasts the query E2 agaisnt all E2s in PaVE to start the process of
    finding the splice acceptor postion for E1^E4 and E8^E2

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: file path to blast output
    """
    E2_trans = str(virus['E2'][-1])
    ID = virus['accession']
    out_dir = args['program_files_dir']
    data_dir = args['data_dir']

    splice_acceptor_dir = os.path.join(out_dir, 'splice_acceptor')
    if not os.path.isdir(splice_acceptor_dir):
        os.makedirs(splice_acceptor_dir)

    blast_subject = os.path.join(data_dir, 'splice_acceptor_blast_new.fa')
    blast_out = os.path.join(splice_acceptor_dir,
                             'blast_result_splice_acceptor.tab')

    if os.path.isfile(blast_out):
        os.remove(blast_out)

    query_file = os.path.join(splice_acceptor_dir, 'query.fa')

    with open(query_file, 'a') as query:
        query.write('>{}\n'.format(ID))
        query.write(E2_trans)
    num_hits = run_blastp(query_file, blast_subject, blast_out)
    return blast_out
# --------------------------------------------------
def locate_known_splice_acceptor(virus, args):
    """
    This function uses the blast results from blast_splice_acceptor() to find the
    splice acceptor information in the E2 that was the closest blast match

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: If the splice acceptor information is found then the known splice start
    and a dictionary containing the known E2 protein information is returned. If the
    splice acceptor information is not found, which means there is a high chance the
    query virus does not have a E1^E4 or E8^E2 then False is returned.
    """
    blast_out = blast_splice_acceptor(virus, args)
    data_dir = args['data_dir']
    query = ''

    with open(blast_out) as blast_result:
        row = blast_result.readline()
        query=row.split("\t")[1]

    splice_start_known = -1
    known_E2 = {}
    csv_database = os.path.join(data_dir, 'PaVE.csv')

    with open(csv_database, 'r') as csvfile:
        read = csv.DictReader(
            csvfile, fieldnames=['accession', 'gene', 'positions', 'seq'])
        for row in read:
            if row['accession'] == query and row['gene'] == 'E8^E2':
                splice_acceptor_positions = row['positions']
            if row['accession'] == query and row["gene"] == 'E2':
                known_E2[query] = str(row['seq']).lower()
                known_E2_start = str(row["positions"])
    try:
        splice_start_genome = splice_acceptor_positions.split('+')[1]

        splice_start_genome = int(str(splice_start_genome).split('..')[0])

        known_E2_start = int(known_E2_start.split('..')[0])

        splice_start_known = (splice_start_genome - known_E2_start)
        
    except UnboundLocalError:
        print('Query does not have E1^E4 or E8^E2')
        logging.info("Strong Possibility that there is not an E1^E4 or E8^E2 because "
                     "closest blast match does not have either.")
    return splice_start_known, known_E2
# --------------------------------------------------
def align_splice_acceptor(virus, args):
    """
    This function aligns the query and the closest E2 using MUSCLE

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: Biopython alignment object
    """
    E2_seq = virus['E2'][2]
    out_dir = args['program_files_dir']
    ID = virus['accession']
    splice_acceptor_dir = os.path.join(out_dir, 'splice_acceptor')

    unaligned = os.path.join(splice_acceptor_dir, 'unaligned.fa')
    aligned = os.path.join(splice_acceptor_dir, 'aligned.fa')

    if os.path.isfile(unaligned):
        os.remove(unaligned)

    if os.path.isfile(aligned):
        os.remove(aligned)

    known_E2 = virus['known_E2']

    for key in known_E2:
        with open(unaligned, 'a') as sequence_file:
            sequence_file.write(">{}\n".format(ID))
            sequence_file.write("{}\n".format(E2_seq))
            sequence_file.write(">{}\n".format(key))
            sequence_file.write("{}\n".format(known_E2[key]))

    cline = MuscleCommandline(input=unaligned, out=aligned, verbose=False)

    stdout, stderr = cline()
    return aligned
# --------------------------------------------------
def find_splice_acceptor(virus, args):
    """
    This function uses the output from align_splice_acceptor(),
    locate_known_splice_acceptor() and blast_splice_acceptor() to locate the splice
    acceptor start position for E1^E4 and E8^E2

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: the splice acceptor genome start position for E1^E4 and E8^E2
    """
    splice_start_known, known_E2 = locate_known_splice_acceptor(virus, args)
    virus['known_E2'] = known_E2
    
    if splice_start_known < 0:
        return -1
    aligned = align_splice_acceptor(virus, args)

    genome = str(virus['genome'])
    aligned_starts = []
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
    search_seq = unknown_seq[aligned_splice_start:aligned_splice_start +
                                                  50].replace('-', '')
    startE2_nt = re.search(search_seq, str(genome).lower()).start()

    return startE2_nt
# --------------------------------------------------
def blast_spliced_e1_e8(virus,args):
    """
    This function blasts the query E1 against all E1s in PaVE to start the
    process of finding the splice donor position for E1^E4 and E8^E2

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: file path to blast output
    """
    E1_trans = str(virus['E1'][-1])
    ID = virus['accession']
    out_dir = args['program_files_dir']
    data_dir = args['data_dir']
    blastE1E8_dir = os.path.join(out_dir, 'blast_E1E8')
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
    num_hits = run_blastp(query_file, blast_subject, blast_out)
    return blast_out
# --------------------------------------------------
def locate_known_e1_splice_donor(virus, args):
    """
    This function uses the blast results from blast_spliced_e1_e8() to find the
    splice donor information in the E1 that was the closest blast match

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: The known E1 splice donor position and a dictionary with the known E1
    protein information
    """

    blast_out = blast_spliced_e1_e8(virus, args)
    data_dir = args['data_dir']
    blast_options = []
    with open(blast_out) as blast_file:
        blast_result = csv.reader(blast_file, delimiter='\t')
        for row in blast_result:
            blast_options.append(row[1])
    blast_options = blast_options[0:1]
    for options in blast_options:
        query = options
        known_E1 = {}
        csv_database = os.path.join(data_dir, 'PaVE.csv')
        with open(csv_database, 'r') as csvfile:
            read = csv.DictReader(csvfile,
                                  ('accession', 'gene', 'positions', 'seq'))
            for row in read:
                if row['accession'] == query and row['gene'] == 'E1^E4':
                    E1_positions = row['positions']
                if row["accession"] == query and row["gene"] == 'E1':
                    known_E1[query] = str(row['seq']).lower()
                    known_E1_stop = str(row["positions"])
        E1_stop_genome = E1_positions.split('+')[0]
        E1_stop_genome = E1_stop_genome.split('(')[1]
        E1_stop_genome = int(str(E1_stop_genome).split('..')[1])
        known_E1_stop = int(known_E1_stop.split('..')[0])
        E1_stop_known = (E1_stop_genome - known_E1_stop)
    return E1_stop_known, known_E1
# --------------------------------------------------
def align_splice_donor_e1(virus, args):
    """
    This function aligns the query E1 and the closest E1 using MUSCLE

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: Biopython alignment object
    """
    #_, known_E1 = locate_known_e1_splice_donor(virus,args)
    out_dir = args['program_files_dir']
    ID = virus['accession']
    E1_seq = str(virus['E1'][2])
    blastE1E8_dir = os.path.join(out_dir, 'blast_E1E8')
    unaligned = os.path.join(blastE1E8_dir, 'unaligned.fa')
    aligned = os.path.join(blastE1E8_dir, 'aligned.fa')

    if os.path.isfile(unaligned):
        os.remove(unaligned)

    if os.path.isfile(aligned):
        os.remove(aligned)
    known_E1 = virus['known_E1']
    for key in known_E1:
        with open(unaligned, 'a') as sequence_file:
            sequence_file.write(">{}\n".format(ID))
            sequence_file.write("{}\n".format(E1_seq))
            sequence_file.write(">{}\n".format(key))
            sequence_file.write("{}\n".format(known_E1[key]))

    cline = MuscleCommandline(input=unaligned, out=aligned, verbose=False)

    stdout, stderr = cline()

    return aligned
# --------------------------------------------------
def find_e4(E2, genome, position):
    """
    This function finds the E4, used for stop site for E1^E4, CURRENTLY NOT IN USE

    :param E2: The E2 nucleotide sequence
    :param genome: linearized genome based on L1
    :param position: number used to choose different E2 ORFS (based on length)
    :return: Dictionary, key being E4 and the value is a list of the start and stop
    position of the found E4
    """
    E4 = {}
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
# --------------------------------------------------
def find_e1_e4(virus,start_E4_nt, args):
    """
    This function uses the output from blast_spliced_e1_e8(),
    locate_known_e1_splice_donor(), align_splice_donor_e1(), find_e4() to locate the
    splice donor site for E1^E4

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param start_E4_nt: splice acceptor position
    :param args: command line arguments for data_dir etc
    :return: dictionary, key being E1^E4 and the value is a list of all the protein
    information
    """
    E1_E4 = {}
    E1_stop_known, known_E1 = locate_known_e1_splice_donor(virus, args)
    virus['known_E1'] = known_E1
    aligned = align_splice_donor_e1(virus, args)
    genome = virus['genome']
    E1_stop = []
    align_seq = []
    j = 0
    aligned_E1_stop = 0

    for aln in AlignIO.read(aligned, 'fasta'):
        align_seq.append(aln.seq)
    unknown_seq = str(align_seq[0]).lower()
    known_seq = str(align_seq[1]).lower()

    for position in known_seq:
        aligned_E1_stop = aligned_E1_stop + 1
        if position.lower() in ['a', 'c', 't', 'g']:
            j = j + 1
            if j == E1_stop_known:
                break
    E1_stop.append(aligned_E1_stop)
    aligned_E1_stop = E1_stop[-1]
    search_seq = unknown_seq[aligned_E1_stop:aligned_E1_stop + 50].replace(
        '-', '')
    stopE1_nt = (re.search(search_seq, str(genome).lower()).start()) + 1
    startE1_nt = virus['E1'][0]
    e4_part = Seq(genome[start_E4_nt +2:]).translate().split("*")[0]
    stopE4_nt = (start_E4_nt + ((len(e4_part) + 1) * 3)) + 2
    E1_E4_seq = str(genome[startE1_nt - 1:stopE1_nt] +
                    genome[start_E4_nt:stopE4_nt])
    E1_E4_trans = Seq(E1_E4_seq).translate()
    E1_E4['E1^E4'] = [startE1_nt, stopE1_nt, start_E4_nt + 1, stopE4_nt, E1_E4_seq,
                      E1_E4_trans]
    return E1_E4
# --------------------------------------------------
def locate_known_e8_splice_donor(virus, args):
    """
    This function uses the blast results to find the E8 stops positions of the closest
    match

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: known E8 stop position, known E1 sequence
    """
    blast_out = blast_spliced_e1_e8(virus, args)
    data_dir = args['data_dir']
    blast_options = []
    with open(blast_out) as blast_file:
        blast_result = csv.reader(blast_file, delimiter='\t')
        for row in blast_result:
            blast_options.append(row[1])
    blast_options = blast_options[0:1]
    E8_stop_known = -1
    for options in blast_options:
        query = options
        known_E8 = {}
        csv_database = os.path.join(data_dir, 'PaVE.csv')
        with open(csv_database, 'r') as csvfile:
            read = csv.DictReader(csvfile,
                                  ('accession', 'gene', 'positions', 'seq'))
            for row in read:
                if row['accession'] == query and row['gene'] == 'E8^E2':
                    E8_positions = row['positions']
                if row["accession"] == query and row["gene"] == 'E1':
                    known_E8[query] = str(row['seq']).lower()
                    known_E8_stop = str(row["positions"])
            try:
                E8_stop_genome = E8_positions.split('+')[0]
                E8_stop_genome = E8_stop_genome.split('(')[1]
                E8_stop_genome = int(str(E8_stop_genome).split('..')[1])
                known_E8_stop = int(known_E8_stop.split('..')[0])
                E8_stop_known = (E8_stop_genome - known_E8_stop)
            except UnboundLocalError:
                logging.info("There was not a blast match for E1 in the process of "
                             "finding E1^E4 and E8^E2.")
    return E8_stop_known, known_E8
# --------------------------------------------------
def align_splice_donor_e8(virus, args):
    """
    This function aligns the unknown E1 and closest E1 blast match

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: MUSCLE alignment of the unknown E1 and closest E1 blast match
    """
    #_, known_E8 = locate_known_e8_splice_donor(virus, args)
    out_dir = args['program_files_dir']
    ID = virus['accession']
    blastE1E8_dir = os.path.join(out_dir, 'blast_E1E8')
    E1_seq = str(virus['E1'][2])
    unaligned = os.path.join(blastE1E8_dir, 'unaligned.fa')
    aligned = os.path.join(blastE1E8_dir, 'aligned.fa')
    if os.path.isfile(unaligned):
        os.remove(unaligned)
    if os.path.isfile(aligned):
        os.remove(aligned)
    known_E8 = virus['known_E8']
    for key in known_E8:
        with open(unaligned, 'a') as sequence_file:
            sequence_file.write(">{}\n".format(ID))
            sequence_file.write("{}\n".format(E1_seq))
            sequence_file.write(">{}\n".format(key))
            sequence_file.write("{}\n".format(known_E8[key]))

    cline = MuscleCommandline(input=unaligned, out=aligned, verbose=False)
    stdout, stderr = cline()
    return aligned
# --------------------------------------------------
def find_e8(virus, args, e8_start):
    """
    This function finds the E8 portion of E8^E2 with in the +1 frame of E1
    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :param e8_start: corresponds to the index of the E8 start postion list
    :return: start and stop positions of the found E8 portion of E8^E2
    """
    E8_stop_known, known_E8 = locate_known_e8_splice_donor(virus, args)
    virus['known_E8'] = known_E8
    if E8_stop_known < 0:
        return -1, -1
    aligned = align_splice_donor_e8(virus, args)
    align_seq = []
    startE8List = []
    startMotif = 'MKL'
    genome = virus['genome']
    E1_whole = virus['E1']
    E1_seq = str(E1_whole[2])
    E8_stop = []
    j = 0
    aligned_E8_stop = 0
    for aln in AlignIO.read(aligned, 'fasta'):
        align_seq.append(aln.seq)
    unknown_seq = str(align_seq[0]).lower()
    known_seq = str(align_seq[1]).lower()
    for position in known_seq:
        aligned_E8_stop = aligned_E8_stop + 1
        if position.lower() in ['a', 'c', 't', 'g']:
            j = j + 1
            if j == E8_stop_known:
                break
    E8_stop.append(aligned_E8_stop)
    aligned_E8_stop = E8_stop[-1]
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
    if len(startE8List) > 0:
        for i in range(0, len(startE8List)):
            try:
                testStart = startE8List[i] + E1_whole[0]
                check_seq = Seq(genome[testStart - 1:stopE8_nt]).translate()
                if "*" in check_seq:
                    del startE8List[i]
                elif check_seq.startswith(startMotif):
                    startE8_nt = startE8List[i] + E1_whole[0]
                    break
            except IndexError:
                break
        if "*" in check_seq:
            startE8_nt = startE8List[0] + E1_whole[0]
        elif check_seq.startswith(startMotif):
            pass
        else:
            startE8_nt = startE8List[e8_start] + E1_whole[0]
    else:
        startE8_nt = E1_whole[0]
    stopE8_nt = (stopE8 + E1_whole[0]) + 1
    return startE8_nt, stopE8_nt
# --------------------------------------------------
def find_e8_e2(virus, startE2_nt, args):
    """
    This function uses all the return values from the above E8^E2 functions to create a
    dictionary with the E8^E2 protein info

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param startE2_nt: nucleotide start position of the splice acceptor site
    :param args: command line arguments for data_dir etc
    :return: dictionary containing the start and stop positions, sequence and
    translated sequence of the E8^E2
    """

    E8_E2 = {}
    genome = virus['genome']
    stopE2_nt = virus['E2'][1]
    try:
        startE8_nt, stopE8_nt = find_e8(virus,args, 0)
        if startE8_nt and stopE8_nt < 0:
            E8_E2['E8^E2'] = []
            return E8_E2
        e2_of_e8_e2 = Seq(genome[startE2_nt +1:stopE2_nt]).translate()
        if e2_of_e8_e2[-20:] == virus['E2'][3][-20:]:#Checking Splice Acceptor
            E8_E2_seq = Seq(genome[startE8_nt - 1:stopE8_nt] + genome[startE2_nt:stopE2_nt])
            E8_E2_trans = E8_E2_seq.translate()
            if E8_E2_trans[-20:] == virus['E2'][3][-20:]:#Checking all of E8^E2
                E8_E2['E8^E2'] = [startE8_nt, stopE8_nt, startE2_nt + 1, stopE2_nt, E8_E2_seq,
                                  E8_E2_trans]
            else:#Case where the E8 start postion may be the other option
                try:
                    startE8_nt, stopE8_nt = find_e8(virus, args, 1)
                    E8_E2_seq = Seq(genome[startE8_nt - 1:stopE8_nt] +
                                    genome[startE2_nt:stopE2_nt])
                    E8_E2_trans = E8_E2_seq.translate()
                    if E8_E2_trans[-20:] == virus['E2'][3][-20:]:
                        E8_E2['E8^E2'] = [startE8_nt, stopE8_nt, startE2_nt + 1, stopE2_nt,
                                          E8_E2_seq,E8_E2_trans]
                    else:
                        logging.info("Problem finding E8 of E8^E2.The results of E1^E4 and E8^E2 "
                                     "are not displayed")
                        E8_E2['E8^E2'] = []
                except:
                    E8_E2['E8^E2'] = []
        else:
            logging.info("Last 20 aa of E8^E2 and E2 do not match. "
                         "The results of E1^E4 and E8^E2 are not displayed")
            E8_E2['E8^E2'] = []
    except:
        E8_E2['E8^E2'] = []
    return E8_E2
# --------------------------------------------------
def to_gff3(virus, for_user_dir):
    """
    This function creates a General Feature Format file

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param for_user_dir: path to the for_user dictionary that is in the puma_out
    directory
    :return: nothing
    """
    dict = {}
    dict.update(virus)
    del dict['URR']
    name = dict['accession']
    dict['name'] = name
    gff3_out = os.path.join(for_user_dir, '{}.gff3'.format(dict['name']))
    genomelen = len(str(virus['genome']))
    with open(gff3_out, 'a') as out_file:
        out_file.write("##gff-version 3\n")
        out_file.write("##sequence-region {} 1 {}\n".format(
            dict['name'], genomelen))

    for protein in dict:
        if protein == 'name' or protein == 'accession' or protein == 'genome' or \
                protein == 'E1BS' or protein == 'E2BS':
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
#---------------------------------------------------
def to_graphic(virus, for_user_dir):
    """
    This functions creates a graph of the found annotations within the genome

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param for_user_dir: path to the for_user dictionary that is in the puma_out
    directory
    :return: nothing
    """
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    virus_copy = {}
    virus_copy.update(virus)
    genome_length = len(virus_copy['genome'])
    pdf_out = os.path.join(for_user_dir,'{}.pdf'.format(virus_copy['accession']))
    color_choices = ['#924900',
                     '#490092',
                     '#24ff24',
                     '#ffb6db']


    fig, ax = plt.subplots()
    with PdfPages(pdf_out) as pdf:
        for gene in virus_copy:
            if gene == 'name' or gene == 'accession' or gene == 'genome':
                pass
            else:

                if gene == 'E1BS':
                    start = virus_copy[gene][0]
                    length = virus_copy[gene][1] - virus_copy[gene][0]
                    ax.broken_barh([(start, length)], (27, 6), facecolors='#000000')
                elif gene == 'E2BS':
                    for binding_site in virus_copy[gene]:
                        ax.broken_barh([(binding_site,12)], (12, 6),
                                       facecolors='#004949')
                elif gene == 'URR':
                    start = virus_copy[gene][0]
                    length = virus_copy[gene][1] - virus_copy[gene][0]
                    ax.broken_barh([(start, length)], (42, 6), facecolors='#009292')
                elif "E1^E4" in gene:
                    start_1 = virus_copy[gene][0]
                    length_1 = virus_copy[gene][1] - virus_copy[gene][0]
                    start_2 = virus_copy[gene][2]
                    length_2 = virus_copy[gene][3] - virus_copy[gene][2]
                    ax.broken_barh([(start_1, length_1),
                                    (start_2, length_2)],
                                    (57, 6),
                                    facecolors='#ff6db6',
                                    label=gene)
                elif "E8^E2" in gene:
                    start_1 = virus_copy[gene][0]
                    length_1 = virus_copy[gene][1] - virus_copy[gene][0]
                    start_2 = virus_copy[gene][2]
                    length_2 = virus_copy[gene][3] - virus_copy[gene][2]
                    ax.broken_barh([(start_1, length_1),
                                       (start_2, length_2)],
                                   (50, 6),
                                   facecolors='#ffff6d',
                                   label=gene)
                elif "E6" in gene:
                    start = virus_copy[gene][0]
                    length = virus_copy[gene][1] - virus_copy[gene][0]
                    frameNumber = virus_copy[gene][0] % 3
                    if frameNumber == 1:
                        frame = 1
                    elif frameNumber == 2:
                        frame = 2
                    elif frameNumber == 0:
                        frame = 3
                    if frame == 1:
                        ax.broken_barh([(start, length)],
                                       (72, 6),
                                       facecolors='#db6d00',
                                       label=gene)
                    elif frame == 2:
                        ax.broken_barh([(start, length)],
                                       (87, 6),
                                       facecolors='#db6d00',
                                       label=gene)
                    elif frame == 3:
                        ax.broken_barh([(start, length)],
                                       (102.5, 6),
                                       facecolors='#db6d00',
                                       label=gene)
                elif "E5" in gene:
                    start = virus_copy[gene][0]
                    length = virus_copy[gene][1] - virus_copy[gene][0]
                    frameNumber = virus_copy[gene][0] % 3
                    if frameNumber == 1:
                        frame = 1
                    elif frameNumber == 2:
                        frame = 2
                    elif frameNumber == 0:
                        frame = 3
                    if frame == 1:
                        ax.broken_barh([(start, length)],
                                       (72, 6),
                                       facecolors='#006ddb',
                                       label=gene)
                    elif frame == 2:
                        ax.broken_barh([(start, length)],
                                       (87, 6),
                                       facecolors='#006ddb',
                                       label=gene)
                    elif frame == 3:
                        ax.broken_barh([(start, length)],
                                       (102.5, 6),
                                       facecolors='#006ddb',
                                       label=gene)
                elif "E7" in gene:
                    start = virus_copy[gene][0]
                    length = virus_copy[gene][1] - virus_copy[gene][0]
                    frameNumber = virus_copy[gene][0] % 3
                    if frameNumber == 1:
                        frame = 1
                    elif frameNumber == 2:
                        frame = 2
                    elif frameNumber == 0:
                        frame = 3
                    if frame == 1:
                        ax.broken_barh([(start, length)],
                                       (72, 6),
                                       facecolors='#b66dff',
                                       label=gene)
                    elif frame == 2:
                        ax.broken_barh([(start, length)],
                                       (87, 6),
                                       facecolors='#b66dff',
                                       label=gene)
                    elif frame == 3:
                        ax.broken_barh([(start, length)],
                                       (102.5, 6),
                                       facecolors='#b66dff',
                                       label=gene)
                elif "E1" in gene:
                    start = virus_copy[gene][0]
                    length = virus_copy[gene][1] - virus_copy[gene][0]
                    frameNumber = virus_copy[gene][0] % 3
                    if frameNumber == 1:
                        frame = 1
                    elif frameNumber == 2:
                        frame = 2
                    elif frameNumber == 0:
                        frame = 3
                    if frame == 1:
                        ax.broken_barh([(start, length)],
                                       (72, 6),
                                       facecolors='#6db6ff',
                                       label=gene)
                    elif frame == 2:
                        ax.broken_barh([(start, length)],
                                       (87, 6),
                                       facecolors='#6db6ff',
                                       label=gene)
                    elif frame == 3:
                        ax.broken_barh([(start, length)],
                                       (102.5, 6),
                                       facecolors='#6db6ff',
                                       label=gene)
                elif "E2" in gene:
                    start = virus_copy[gene][0]
                    length = virus_copy[gene][1] - virus_copy[gene][0]
                    frameNumber = virus_copy[gene][0] % 3
                    if frameNumber == 1:
                        frame = 1
                    elif frameNumber == 2:
                        frame = 2
                    elif frameNumber == 0:
                        frame = 3
                    if frame == 1:
                        ax.broken_barh([(start, length)],
                                       (72, 6),
                                       facecolors='#b6dbff',
                                       label=gene)
                    elif frame == 2:
                        ax.broken_barh([(start, length)],
                                       (87, 6),
                                       facecolors='#b6dbff',
                                       label=gene)
                    elif frame == 3:
                        ax.broken_barh([(start, length)],
                                       (102.5, 6),
                                       facecolors='#b6dbff',
                                       label=gene)

                elif "L1" in gene:
                    start = virus_copy[gene][0]
                    length = virus_copy[gene][1] - virus_copy[gene][0]
                    frameNumber = virus_copy[gene][0] % 3
                    if frameNumber == 1:
                        frame = 1
                    elif frameNumber == 2:
                        frame = 2
                    elif frameNumber == 0:
                        frame = 3
                    if frame == 1:
                        ax.broken_barh([(start, length)],
                                       (72, 6),
                                       facecolors='#984ea3',
                                       label=gene)
                    elif frame == 2:
                        ax.broken_barh([(start, length)],
                                       (87, 6),
                                       facecolors='#984ea3',
                                       label=gene)
                    elif frame == 3:
                        ax.broken_barh([(start, length)],
                                       (102.5, 6),
                                       facecolors='#984ea3',
                                       label=gene)
                elif "L2" in gene:
                    start = virus_copy[gene][0]
                    length = virus_copy[gene][1] - virus_copy[gene][0]
                    frameNumber = virus_copy[gene][0] % 3
                    if frameNumber == 1:
                        frame = 1
                    elif frameNumber == 2:
                        frame = 2
                    elif frameNumber == 0:
                        frame = 3
                    if frame == 1:
                        ax.broken_barh([(start, length)],
                                       (72, 6),
                                       facecolors='#920000',
                                       label=gene)
                    elif frame == 2:
                        ax.broken_barh([(start, length)],
                                       (87, 6),
                                       facecolors='#920000',
                                       label=gene)
                    elif frame == 3:
                        ax.broken_barh([(start, length)],
                                       (102.5, 6),
                                       facecolors='#920000',
                                       label=gene)


                else:
                    color = random.choice(color_choices)
                    color_choices.remove(color)
                    start = virus_copy[gene][0]
                    length = virus_copy[gene][1] - virus_copy[gene][0]
                    frameNumber = virus_copy[gene][0] % 3
                    if frameNumber == 1:
                        frame = 1
                    elif frameNumber == 2:
                        frame = 2
                    elif frameNumber == 0:
                        frame = 3
                    if frame == 1:
                        ax.broken_barh([(start, length)],
                                       (72, 6),
                                       facecolors=color,
                                       label=gene)
                    elif frame == 2:
                        ax.broken_barh([(start, length)],
                                       (87, 6),
                                       facecolors=color,
                                       label=gene)
                    elif frame == 3:
                        ax.broken_barh([(start, length)],
                                       (102.5, 6),
                                       facecolors=color,
                                       label=gene)
        ax.set_ylim(0, 120)
        ax.set_xlim(0, genome_length)
        ax.set_xlabel('Genome Position')
        ax.set_yticks([15, 30, 45, 60, 75, 90, 105])
        ax.set_yticklabels(['E2BS', 'E1BS', 'URR', "SPLICED", 'ORF 1', 'ORF 2', 'ORF 3'])
        plt.legend(loc='lower right', ncol=3, shadow=True)
        plt.title(virus_copy['name'] + " Linearized Based on L1")
        pdf.savefig()
        plt.close()
        firstPage = plt.figure(figsize=(6, 5))
        firstPage.clf()
        text = []
        bs_list = []
        txt = ""
        genome_name = "Genome    1-{}".format(len(virus_copy['genome']))
        for gene in virus_copy:
            if gene == 'name' or gene == 'accession' or gene == 'genome':
                pass
            else:
                if "^" in gene:
                    start_1 = virus_copy[gene][0]
                    stop_1 = virus_copy[gene][1]
                    start_2 = virus_copy[gene][2]
                    stop_2 = virus_copy[gene][3]
                    text.append("{}     {}-{}, {}-{}".format(gene, start_1, stop_1,
                                                               start_2, stop_2))
                elif "BS" in gene:
                    if gene == 'E1BS':
                        bs_list.append("{}     {}-{}".format(gene,virus_copy[gene][0],
                                                             virus_copy[gene][1]))
                    elif gene == 'E2BS':
                        for bs in virus_copy[gene]:
                            bs_list.append("{}     {}-{}".format(gene, bs, bs+12))
                else:
                    text.append("{}     {}-{}".format(gene, virus_copy[gene][0],
                                                      virus_copy[gene][1]))
        text.sort(key= len)
        text.insert(0,genome_name)
        text.extend(bs_list)
        txt = "\n".join(text)
        firstPage.text(0, 1, txt, horizontalalignment='left',
                       verticalalignment='top',fontsize=12)
        pdf.savefig()
        plt.close()

    return
# --------------------------------------------------
def to_csv(virus, for_user_dir):
    """
    This function creates a comma seperated values file

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param for_user_dir: path to the for_user dictionary that is in the puma_out
    directory
    :return: nothing
    """
    virus_copy = {}
    virus_copy.update(virus)

    csv_out = os.path.join(for_user_dir, '{}.csv'.format(virus_copy['accession']))
    try:
        if virus_copy['E2BS'] == ['No E2BS found']:
            del virus_copy['E2BS']
    except KeyError:
        pass

    with open(csv_out, 'a') as out:
        out_file = csv.writer(out)

        for value in virus_copy:
            if value == 'genome':
                out_file.writerows([[
                    virus_copy['accession'], 'CG', "",
                    str(virus_copy[value]).lower()
                ]])
            elif value == 'URR':
                try:
                    out_file.writerows([[
                        virus_copy['accession'], value,
                        'join(' + str(virus_copy[value][0]) + '..' +
                        str(virus_copy[value][1]) + '+' +
                        str(virus_copy[value][2]) + ".." +
                        str(virus_copy[value][3]) + ')', virus_copy[value][4]
                    ]])
                except IndexError:
                    out_file.writerows([[
                        virus_copy['accession'], value,
                        str(virus_copy[value][0]) + '..' +
                        str(virus_copy[value][1]), virus_copy[value][2]
                    ]])
            elif value == 'E1BS':
                try:
                    out_file.writerows([[
                        virus_copy['accession'], value,
                        'join(' + str(virus_copy[value][0]) + '..' +
                        str(virus_copy[value][1]) + '+' +
                        str(virus_copy[value][2]) + ".." +
                        str(virus_copy[value][3]) + ')', virus_copy[value][4]
                    ]])
                except IndexError:
                    out_file.writerows([[
                        virus_copy['accession'], value,
                        str(virus_copy[value][0]) + '..' +
                        str(virus_copy[value][1]), virus_copy[value][2]
                    ]])
            elif value == 'E2BS':
                for i in range(0, len(virus_copy[value]), 1):
                    out_file.writerows([[
                        virus_copy['accession'], value,
                        str(virus_copy[value][i]) + '..' +
                        str(virus_copy[value][i] + 12),
                        str(virus_copy['genome'][virus_copy[value][i] -
                                                  1:virus_copy[value][i] +
                                                  11]).lower()
                    ]])
            elif value == 'E1^E4':
                out_file.writerows([[
                    virus_copy['accession'], value,
                    'join(' + str(virus_copy[value][0]) + '..' +
                    str(virus_copy[value][1]) + '+' +
                    str(virus_copy[value][2]) + ".." +
                    str(virus_copy[value][3]) + ')', virus_copy[value][4],
                    virus_copy[value][5]
                ]])
            elif value == 'E8^E2':
                out_file.writerows([[
                    virus_copy['accession'], value,
                    'join(' + str(virus_copy[value][0]) + '..' +
                    str(virus_copy[value][1]) + '+' +
                    str(virus_copy[value][2]) + ".." +
                    str(virus_copy[value][3]) + ')', virus_copy[value][4],
                    virus_copy[value][5]
                ]])
            elif value == 'name':
                pass
            elif value == 'accession':
                pass
            else:
                out_file.writerows([[
                    virus_copy['accession'], value,
                    str(virus_copy[value][0]) + '..' +
                    str(virus_copy[value][1]), virus_copy[value][2],
                    virus_copy[value][3]
                ]])
    return
# --------------------------------------------------
def to_genbank(virus, for_user_dir):
    """
    This functions utilizes biopython to create a GenBank formatted output file

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param for_user_dir: path to the for_user dictionary that is in the puma_out
    directory
    :return: nothing
    """
    ID = virus['accession']
    name = virus['name']
    if len(name) < 1:
        name = ID
        if len(name) > 15:
            name = name[0:15]
            name = name.replace(" ","_")
    else:
        if len(name) > 15:
            name = name[0:15]
            name = name.replace(" ", "_")


    sequence_string = virus['genome']
    sequence_object = Seq(sequence_string, IUPAC.ambiguous_dna)
    record = SeqRecord(sequence_object,
                       id=ID,  # random accession number
                       name=name,  # replace this with a sensible name
                       description='GenBank file for %s generated by PuMA' % (virus[
                           'name']))
    for gene in virus:
        if gene == 'name' or gene == 'accession' or gene == 'genome':
            pass
        elif "^" not in gene:
            if "E1BS" == gene:
                start = virus[gene][0] - 1
                end = virus[gene][1]
                notes = {"note": gene}
                feature = SeqFeature(FeatureLocation(start=start, end=end),
                                     type='misc_feaure', qualifiers=notes)
                record.features.append(feature)
            elif "E2BS" == gene:
                for i in range(0, len(virus['E2BS'])):
                    start = virus[gene][i] -1
                    end = virus[gene][i] + 12
                    notes = {"note": gene}
                    feature = SeqFeature(FeatureLocation(start=start, end=end),
                                         type='misc_feaure', qualifiers=notes)
                    record.features.append(feature)
            elif gene == 'URR':
                start = virus[gene][0] -1
                end = virus[gene][1]
                notes = {"note": gene}
                feature = SeqFeature(FeatureLocation(start=start, end=end),
                                     type='misc_feaure', qualifiers=notes)
                record.features.append(feature)
            else:
                start = virus[gene][0] -1
                end = virus[gene][1]
                notes = {"gene": gene, "protein_id": ID + "_" + gene,
                    "translation": sequence_object[start:end].translate()}
                feature = SeqFeature(FeatureLocation(start=start, end=end), type='CDS',
                                     qualifiers=notes)
                record.features.append(feature)
        else:
            start = virus[gene][0] -1
            s_d = virus[gene][1]
            s_a = virus[gene][2] -1
            end = virus[gene][3]
            f_1 = FeatureLocation(start, s_d)
            f_2 = FeatureLocation(s_a, end)
            join = CompoundLocation([f_1, f_2])
            spliced = sequence_object[start:s_d] + sequence_object[s_a:end]
            notes = {"gene": gene, "protein_id": ID + "_" + gene,
                "translation": spliced.translate()}
            feature = SeqFeature(join, type='CDS', qualifiers=notes)
            record.features.append(feature)
    gb_out = os.path.join(for_user_dir, '{}.gb'.format(virus['accession']))
    output_file = open(gb_out, 'a')
    SeqIO.write(record, output_file, 'genbank')
    return
# --------------------------------------------------
def to_sequin(virus, for_user_dir):
    """
    This function uses the created genbank file to create .fsa and .tbl files to
    help with the genbank submission process

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param for_user_dir: path to the for_user dictionary that is in the puma_out
    directory
    :return: nothing
    """
    genbank = os.path.join(for_user_dir, '{}.gb'.format(virus['accession']))
    if os.path.isfile(genbank):
        pass
    else:
        logging.info("Problem with genbank file. Cannot create .fsa and .tbl files for "
                     "genbank submission process")
        return
    genbank_sub_dir = os.path.join(for_user_dir, 'genbank_submission')
    if not os.path.isdir(genbank_sub_dir):
        os.makedirs(genbank_sub_dir)

    for r in SeqIO.parse(genbank,"genbank"):
       fasta_name = r.id + '.fsa'
       table_name = r.id + '.tbl'
    fasta_file = os.path.join(genbank_sub_dir, '{}'.format(fasta_name))
    table_file = os.path.join(genbank_sub_dir, '{}'.format(table_name))
    for r in SeqIO.parse(genbank,"genbank"):
      with open(fasta_file, 'w') as fasta, open(table_file, 'w') as table:
       fasta.write(">{} [organism= {}]\n{}".format(r.id,r.id,r.seq))
       table.write(">Feature {}\n".format(r.id))
       for (index, feature) in enumerate(r.features):
        if feature.type ==  'CDS':
            m=re.search('\[(\d*)\:(\d*)\]',str(feature.location))
            if m:
               start = int(m.groups()[0])
               stop= int(m.groups()[1])
            table.write("{}\t{}\tgene\n\t\t\tgene\t{}\n".format(start+1,
                        stop,feature.qualifiers['gene'][0]))
            table.write("{}\t{}\tCDS\n\t\t\tproduct\t{}\n\t\t\tgene\t{}\n\t\t"
                        "\tcodon_start\t{}\n".format(int(m.groups()[0])+1, m.groups()[1],
                        feature.qualifiers['gene'][0], feature.qualifiers['gene'][0], "1"))
        if feature.type ==  'mRNA':
               m=re.search('\[(\d*)\:(\d*)\].*\[(\d*)\:(\d*)\]',str(feature.location))
               if m:
                  start, SD, SA, stop = m.groups()[0],m.groups()[1],m.groups()[2],m.groups()[3]
                  table.write("{}\t{}\tmRNA\n{}\t{}\n\t\t\tproduct\t{}\n"
                              .format(int(start)+1, SD, int(SA)+1, stop,
                                      feature.qualifiers['gene'][0]))
        if "misc" in feature.type:
            m=re.search('\[(\d*)\:(\d*)\]',str(feature.location))
            if m:
               start = int(m.groups()[0])
               stop= int(m.groups()[1])
               table.write("{}\t{}\tmisc_feature\n\t\t\tnote\t{}\n"
                           .format(start+1,stop, feature.qualifiers['note'][0]))
    read_me = os.path.join(genbank_sub_dir, 'README.txt')
    with open(read_me, "w") as rd:
        rd.write("The .fsa and .tbl files in this folder are for "
                 "the genbank submission process.\n"
                 "Please refer to https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/ to "
                 "utilize the files created with PuMA along with user created files to "
                 "start the genbank submission process.\n"
                 "\n")
    return
# --------------------------------------------------
def print_genome_info(virus):
    """
    This function prints all the found annotations

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :return: nothing
    """
    print("\nThis is the gene information for {} after making L1 end of Genome:".
        format(virus['accession']))
    if "known_E2" in virus:
        del virus["known_E2"]
    if "known_E1":
        del virus['known_E1']
    if "known_E8":
        del virus['known_E8']
    for name in virus:
        if name == 'genome':
            pass
        elif name == 'accession':
            pass
        elif name == 'name':
            pass
        elif name == 'E2BS':
            print("\n{} E2 binding sites found:".format(len(virus['E2BS'])))
            for i in range(0, len(virus['E2BS'])):
                print('\n{} start and stop position:\n{},{}\n'.format(
                    name, virus[name][i], virus[name][i] + 12))
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
                        print('{} translated sequence:\n{}\n'.format(
                            name, virus[name][5][:-1]))
                else:
                    print('\n{} start and stop position:\n{},{}\n'.format(
                        name, virus[name][0], virus[name][1]))
                    print('{} sequence:\n{}\n'.format(name, virus[name][2]))
                    if name != 'URR':
                        print('{} translated sequence:\n{}\n'.format(
                            name, virus[name][3][:-1]))
            except IndexError:
                print('\n{} start and stop position:\n{},{}\n'.format(
                    name, virus[name][0], virus[name][1]))
                print('{} sequence:\n{}\n'.format(name, virus[name][2]))
                if name != 'URR':
                    print('{} translated seqeunce:\n{}\n'.format(
                        name, virus[name][3][:-1]))
    print("All annotations displayed above.\n")
    return
# --------------------------------------------------
def validate_args(args):
    """
    Validates command line arguments

    :param args: command line arguments
    :return: dictionary, keys are the command argument names, values are the command
    line values
    """

    input_file = args.get('input')
    out_dir = args.get('out_dir')
    data_dir = args.get('data_dir')
    input_format = args.get('format').lower()
    min_prot_len = args.get('min_prot_len')
    e_value = args.get('evalue')
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    return {
        'input_file': input_file,
        'out_dir': out_dir,
        'data_dir': data_dir,
        'input_format': input_format,
        'min_prot_len': min_prot_len,
        'e_value': e_value
    }
# --------------------------------------------------
def puma_output(virus, args):
    """
    This functions calls all the functions that are related to creating the output files

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param startE2_nt: nucleotide start position of the splice acceptor site
    :param args: command line arguments for data_dir etc
    :return: nothing
    """

    print_genome_info(virus)
    to_graphic(virus, args['for_user_dir'])
    to_csv(virus, args['for_user_dir'])
    to_gff3(virus, args['for_user_dir'])
    to_genbank(virus, args['for_user_dir'])
    to_sequin(virus, args['for_user_dir'])
    return
# --------------------------------------------------
def run(args):
    """main"""
    genomes_from_file = {}
    out_dir = args['out_dir']
    logging.info('run = {}\n'.format(args))
    warnings.simplefilter('ignore', BiopythonWarning)
    args = validate_args(args)
    for seq_record in SeqIO.parse(args['input_file'], args['input_format']):
        original_genome = seq_record.seq
        name = seq_record.description
        try:
            full_name = name.split("|")[1]
            ID = name.split("|")[0]
            genomes_from_file[ID] = [full_name, original_genome]
        except IndexError:
            raise Exception("Please format the input fasta file correctly.\n"
                            "The formatting is:\n"
                            ">Short Name|Full Name\n"
                            "Sequence (atcg)\n")

    logging.info(
        "Total of {} genome(s) in input file\n\n".format(
            len(genomes_from_file)))
    for name in genomes_from_file:
        virus = {}
        virus_dir = os.path.join(out_dir, name)
        if not os.path.isdir(virus_dir):
            os.makedirs(virus_dir)
        for_user_dir = os.path.join(virus_dir, 'for_user')
        if not os.path.isdir(for_user_dir):
            os.makedirs(for_user_dir)
        args['for_user_dir'] = for_user_dir
        program_dir = os.path.join(virus_dir, 'program_files')
        if not os.path.isdir(program_dir):
            os.makedirs(program_dir)
        args['program_files_dir'] = program_dir
        logging.info("\n\nBelow, find info from annotation of {}" .format(name))
        virus['accession'] = name
        virus['name'] = genomes_from_file[name][0]
        original_genome = genomes_from_file[name][1]
        altered_genome = linearize_genome(original_genome, args)
        virus['genome'] = str(altered_genome).lower()
        virus.update(identify_main_proteins(altered_genome, args))
        if 'E1' in virus.keys():
            verified_gene = verify_gene(virus,'E1',args)
            del virus['E1']
            virus['E1'] = verified_gene['E1']
        else:
            logging.warning("No E1 found")

        if 'E2' in virus.keys():
            verified_gene = verify_gene(virus,'E2',args)
            del virus['E2']
            virus['E2'] = verified_gene['E2']
        else:
            logging.warning("No E2 found")
        if 'E5' in virus.keys():
            verified_gene = verify_gene(virus, 'E5', args)
            del virus['E5']
            virus['E5'] = verified_gene['E5']
        else:
            logging.warning("No E5 found")

        if 'E5_ALPHA' in virus.keys():
            verified_gene = verify_gene(virus, 'E5_ALPHA', args)
            del virus['E5_ALPHA']
            virus['E5_ALPHA'] = verified_gene['E5_ALPHA']
        else:
            logging.warning("No E5_ALPHA found")

        if 'E5_BETA' in virus.keys():
            verified_gene = verify_gene(virus, 'E5_BETA', args)
            del virus['E5_BETA']
            virus['E5_BETA'] = verified_gene['E5_BETA']
        else:
            logging.warning("No E5_BETA found")

        if 'E5_GAMMA' in virus.keys():
            verified_gene = verify_gene(virus, 'E5_GAMMA', args)
            del virus['E5_GAMMA']
            virus['E5_GAMMA'] = verified_gene['E5_GAMMA']
        else:
            logging.warning("No E5_GAMMA found")

        if 'E5_DELTA' in virus.keys():
            verified_gene = verify_gene(virus, 'E5_DELTA', args)
            del virus['E5_DELTA']
            virus['E5_DELTA'] = verified_gene['E5_DELTA']
        else:
            logging.warning("No E5_DELTA found")

        if 'E5_EPSILON' in virus.keys():
            verified_gene = verify_gene(virus, 'E5_EPSILON', args)
            del virus['E5_EPSILON']
            virus['E5_EPSILON'] = verified_gene['E5_EPSILON']
        else:
            logging.warning("No E5_EPSILON found")

        if 'E5_ZETA' in virus.keys():
            verified_gene = verify_gene(virus, 'E5_ZETA', args)
            del virus['E5_ZETA']
            virus['E5_ZETA'] = verified_gene['E5_ZETA']
        else:
            logging.warning("No E5_ZETA found")
        if 'E6' in virus.keys():
            verified_gene = verify_gene(virus,'E6',args)
            del virus['E6']
            virus['E6'] = verified_gene['E6']
        else:
            logging.warning("No E6 found")
        if 'E7' in virus.keys():
            verified_gene = verify_gene(virus,'E7',args)
            del virus['E7']
            virus['E7'] = verified_gene['E7']
        else:
            logging.warning("No E7 found")
        if 'E10' in virus.keys():
            verified_gene = verify_gene(virus, 'E10', args)
            del virus['E10']
            virus['E10'] = verified_gene['E10']
        else:
            logging.warning("No E10 found")

        if 'L1' in virus.keys():
            verified_gene = verify_gene(virus,'L1',args)
            del virus['L1']
            virus['L1'] = verified_gene['L1']
        else:
            logging.warning("No L1 found")
        if 'L2' in virus.keys():
            verified_gene = verify_gene(virus,'L2',args)
            del virus['L2']
            virus['L2'] = verified_gene['L2']
        else:
            logging.warning("No L2 found")
        if 'E2' and 'L2' in virus.keys():
            virus.update(identify_e5_variants(virus, args))
        else:
            logging.info("E5 variants function not executed because PuMA did not "
                         "find E2, L2 or both.")
        virus.update(find_urr(virus))
        E1BS = find_e1bs(virus, args)
        if E1BS:
            virus.update(E1BS)
        else:
            logging.info("NO E1BS found")
        E2BS = find_e2bs(virus, args)
        if len(E2BS['E2BS']) > 0:
            virus.update(E2BS)
        else:
            logging.info("NO E2BS found")
        if 'E2' in virus.keys():
            start_splice_site = find_splice_acceptor(virus, args)
            if start_splice_site > 0:
                E8_E2 = find_e8_e2(virus, start_splice_site, args)
                if len(E8_E2['E8^E2']) > 0:
                    E1_E4 = find_e1_e4(virus, start_splice_site, args)
                    virus.update(E8_E2)
                    virus.update(E1_E4)
        puma_output(virus,args)
    return 1
# ------------------------------------------------------------------------------