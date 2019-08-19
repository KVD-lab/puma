"""
puma poly library

Authors: Josh Pace, Ken Youens-Clark, Koenraad Van Doorslaer
University of Arizona, KVD Lab & Hurwitz Lab
PuMA Poly 0.1 7/25/19
"""

import os
import glob
import re
import csv
import time
import operator
import argparse
import sys
import warnings
import random
import itertools
import shutil
import logging
import numpy as np
import pandas as pd
from Bio import SeqIO, GenBank, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import MafftCommandline
from Bio import BiopythonWarning
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from subprocess import getstatusoutput
from distutils.spawn import find_executable
from pprint import pprint as pp
from dire import die

# ----------------------------------------------------------------------------------------
def find_orfs_with_trans(seq, trans_table, min_protein_length):
    """
        This functions finds all the open reading frames and translates them

        :param seq: PV genome
        :param min_protein_length: default is 25 or can be changed by user, ORFs shorter
        then value will be discarded
        :return: a dictionary where the key is the ORF and the value is the start of the
        ORF in the genome
        """
    ORFs = {}
    seq_len = len(seq)
    has_m = re.compile('M')
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end-aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame+aa_start*3
                        end = min(seq_len,frame+aa_end*3+3)
                    else:
                        start = seq_len-frame-aa_end*3-3
                        end = seq_len-frame-aa_start*3
                    aa_seq = trans[aa_start:aa_end]
                    if has_m.search(aa_seq):
                        ORFs[aa_seq] = start
                aa_start = aa_end+1

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
    orfs = find_orfs_with_trans(genome, 1, args['min_prot_len'])
    if not orfs:
        raise Exception('No ORFs, must stop.')

    orfs_fa = os.path.join(args['program_files_dir'], 'orfs.fa')
    orfs_fh = open(orfs_fa, 'wt')

    for orf in orfs:
        orfs_fh.write('\n'.join(['>' + str(orfs[orf]), orf, '']))
    orfs_fh.close()

    blast_sub = os.path.join(args['data_dir'], 'main_blast.fa')
    blast_out = os.path.join(args['program_files_dir'], 'blast_results_main.tab')

    if os.path.isfile(blast_out):
        os.remove(blast_out)
    try:
        num_hits = run_blastp(orfs_fa, blast_sub, blast_out)
        #print('Got {} hits!'.format(num_hits))
    except Exception as e:
        die('EXCEPTION: {}'.format(e))


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
        print("no blast output for main proteins")
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
    strand = ""
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
                    strand = "+"
                    if "*" in translated[:-1]:
                        diff = abs(protein_start[start] - real_start)
                        real_start = real_start - diff
                        sequence = str(
                            genome[int(real_start):int(end)]).lower()
                        sequence = Seq(sequence).reverse_complement()
                        translated = sequence.translate()
                        strand = "-"
                    if len(translated) < 25:
                        logging.info("{} is of length {} and therefore is too short to be "
                                  "considered viable".format(seq, len(translated)))
                    else:
                        found_proteins[seq] = [
                            int(real_start) + 1,
                            int(end), sequence, translated, strand
                        ]
                except AttributeError:
                    pass

    return found_proteins

def blast_splice_acceptor(virus, args):
    """
    This function blasts the query Large T  against all Large T in csv to start the
    process of finding the splice acceptor postion for Large T

    :param virus: dictionary that has all found proteins
    :param args: command line arguments for data_dir etc
    :return: file path to blast output
    """
    large_t_orf = str(virus['Large_T_orf_trans'])
    ID = virus['accession']
    out_dir = args['program_files_dir']
    data_dir = args['data_dir']

    splice_acceptor_dir = os.path.join(out_dir, 'splice_acceptor')
    if not os.path.isdir(splice_acceptor_dir):
        os.makedirs(splice_acceptor_dir)

    blast_subject = os.path.join(data_dir, 'large_t_blast_half.fa')
    blast_out = os.path.join(splice_acceptor_dir,
                             'blast_result_splice_acceptor.tab')

    if os.path.isfile(blast_out):
        os.remove(blast_out)

    query_file = os.path.join(splice_acceptor_dir, 'query.fa')

    with open(query_file, 'a') as query:
        query.write('>{}\n'.format(ID))
        query.write(large_t_orf)

    try:
        num_hits = run_blastp(query_file, blast_subject, blast_out)
        # print('Got {} hits!'.format(num_hits))
    except Exception as e:
        die('EXCEPTION: {}'.format(e))

    return blast_out

# ----------------------------------------------------------------------------------------
def locate_known_splice_acceptor(virus, args):
    """
    This function uses the blast results from blast_splice_acceptor() to find the
    splice acceptor information in the Large T that was the closest blast match

    :param virus: dictionary that has all found proteins so far
    :param args: command line arguments for data_dir etc
    :return: If the splice acceptor information is found then the known splice start
    and a dictionary containing the known E2 protein information is returned. If the
    splice acceptor information is not found, which means there is a high chance the
    query virus does not have a E1^E4 or E8^E2 then False is returned.
    """
    blast_out = blast_splice_acceptor(virus, args)
    data_dir = args['data_dir']
    query = ''
    with open(blast_out) as blast_file:
        blast_result = csv.reader(blast_file, delimiter='\t')
        for row in blast_result:
            query = row[1]
            break
    splice_start_known = -1
    known_large_t2 = {}
    csv_database = os.path.join(data_dir, 'results_poly.csv')

    with open(csv_database, 'r') as csvfile:
        read = csv.DictReader(
            csvfile, fieldnames=['name', 'gene', 'positions', 'seq','strand'])
        for row in read:
            if query in row['name'] and row['gene'] == 'Large T':
                splice_acceptor_positions = row['positions']
            if query in row['name'] and row['gene'] == 'CG':
                genome = str(row['seq']).lower()
        try:

            splice_acceptor_positions = splice_acceptor_positions.split('+')[1]
            stop_genome = int(str(splice_acceptor_positions).split('..')[0])
            splice_acceptor = str(splice_acceptor_positions).split('..')[1]
            splice_acceptor = int(str(splice_acceptor).split(")")[0])
            known_large_t2[query] = Seq(genome[stop_genome -1:splice_acceptor + 2
                                        ]).reverse_complement()
            splice_start_known = (splice_acceptor - stop_genome )
        except UnboundLocalError:
            pass
    return splice_start_known, known_large_t2
# ----------------------------------------------------------------------------------------
def align_splice_acceptor(virus, args):
    """
    This function aligns the query and the closest E2 using MUSCLE

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: Biopython alignment object
    """
    _, known_large_t2 = locate_known_splice_acceptor(virus, args)
    large_t_orf_seq = virus['Large_T_orf_seq']
    out_dir = args['program_files_dir']
    ID = virus['accession']
    splice_acceptor_dir = os.path.join(out_dir, 'splice_acceptor')

    unaligned = os.path.join(splice_acceptor_dir, 'unaligned.fa')
    aligned = os.path.join(splice_acceptor_dir, 'aligned.fa')

    if os.path.isfile(unaligned):
        os.remove(unaligned)

    if os.path.isfile(aligned):
        os.remove(aligned)

    for key in known_large_t2:
        with open(unaligned, 'a') as sequence_file:
            sequence_file.write(">{}\n".format(ID))
            sequence_file.write("{}\n".format(large_t_orf_seq))
            sequence_file.write(">{}\n".format(key))
            sequence_file.write("{}\n".format(known_large_t2[key]))

    cline = MuscleCommandline(input=unaligned, out=aligned, verbose=False)

    stdout, stderr = cline()

    return aligned
# ----------------------------------------------------------------------------------------
def percent_identity(alignment_file):
    alignment = AlignIO.parse(alignment_file,'fasta')
    for align_record in alignment:
        y=0
        z=0
        for r in range(len(align_record[0])):
            if list(align_record[0])[r] == list(align_record[1])[r]:
                y=y+1
    return 100*(float(y)/float(r))
# ----------------------------------------------------------------------------------------
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
    splice_start_known, _ = locate_known_splice_acceptor(virus, args)
    if splice_start_known < 0:
        return -1

    aligned = align_splice_acceptor(virus, args)
    genome = str(virus['genome']).lower()
    rc_genome = str(Seq(genome).reverse_complement())
    large_t_stop = virus['Large_T'][0]
    align_seq = []
    for aln in AlignIO.read(aligned, 'fasta'):
        align_seq.append(aln.seq)

    unknown_seq = str(align_seq[0]).lower()
    known_seq = str(align_seq[1]).lower()

    num_dashes = 0
    for position in known_seq:
        if position == '-':
            num_dashes = num_dashes + 1
        else:
            position_index = known_seq.index(position)
            if known_seq[position_index:position_index+2].lower() == 'ag':
                break
    if virus['Large_T_orf_seq'][num_dashes:num_dashes+2].lower() == 'ag':
        search_seq = virus['Large_T_orf_seq'][num_dashes:]
        splice_acceptor_gene = re.search(str(search_seq),str(rc_genome)).start()
        splice_acceptor = (len(genome) - splice_acceptor_gene) -2
    else:
        raise Exception("There was a problem finding the splice acceptor for Large T")
    return splice_acceptor
# ----------------------------------------------------------------------------------------
def blast_for_known_spliced_donor(virus,args):
    """
    This function blasts the query small t against all E1s in the csv file start the
    process of finding the splice donor position for Large T

    :param virus: dictionary that has all found proteins so far
    :param args: command line arguments for data_dir etc
    :return: file path to blast output
    """
    small_t_trans = str(virus['small_t'][-2])
    ID = virus['accession']
    out_dir = args['program_files_dir']
    data_dir = args['data_dir']

    blast_small_t_dir = os.path.join(out_dir, 'blast_small_t')
    if not os.path.isdir(blast_small_t_dir):
        os.makedirs(blast_small_t_dir)

    blast_subject = os.path.join(data_dir, 'small_t_blast_half.fa')
    blast_out = os.path.join(blast_small_t_dir, 'blast_result.tab')

    if os.path.isfile(blast_out):
        os.remove(blast_out)

    query_file = os.path.join(blast_small_t_dir, 'query.fa')

    with open(query_file, 'a') as query:
        query.write('>{}\n'.format(ID))
        query.write(small_t_trans)

    try:
        num_hits = run_blastp(query_file, blast_subject, blast_out)
        #print('Got {} hits!'.format(num_hits))
    except Exception as e:
        die('EXCEPTION: {}'.format(e))
    return blast_out
# ----------------------------------------------------------------------------------------
def locate_known_splice_donor(virus, args):
    """
    This function uses the blast results from blast_spliced_e1_e8() to find the
    splice donor information in the E1 that was the closest blast match

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: The known E1 splice donor position and a dictionary with the known E1
    protein information
    """

    blast_out = blast_for_known_spliced_donor(virus, args)
    data_dir = args['data_dir']
    blast_options = []
    with open(blast_out) as blast_file:
        blast_result = csv.reader(blast_file, delimiter='\t')
        for row in blast_result:
            blast_options.append(row[1])
    blast_options = blast_options[0:1]
    splice_sites = []
    known_splice_donor = -1
    for options in blast_options:
        query = options
        known_small_t = {}
        csv_database = os.path.join(data_dir, 'results_poly.csv')

        with open(csv_database, 'r') as csvfile:
            read = csv.DictReader(csvfile,
                                  ('name', 'gene', 'positions', 'seq'))
            for row in read:
                if query in row['name'] and row['gene'] == 'small t':
                    small_t_positions = row['positions']
                    known_small_t[query] = str(row['seq']).lower()
                if query in row['name'] and row['gene'] == 'Large T':
                    large_t_positions = str(row['positions'])
        try:
            small_t_start_genome = int(str(small_t_positions).split('..')[1])
            splice_donor = large_t_positions.split("(")[1]
            splice_donor = int(splice_donor.split("..")[0])
            known_splice_donor = small_t_start_genome - splice_donor
        except UnboundLocalError:
            pass
    return known_splice_donor, known_small_t
# ----------------------------------------------------------------------------------------
def align_splice_donor(virus, args):
    """
    This function aligns the query E1 and the closest E1 using MUSCLE

    :param virus: dictionary that has all found proteins so far and linearized based on L1
    genome
    :param args: command line arguments for data_dir etc
    :return: Biopython alignment object
    """
    _, known_small_t = locate_known_splice_donor(virus,args)
    out_dir = args['program_files_dir']
    ID = virus['accession']
    small_t_seq = str(virus['small_t'][2])
    blast_small_t_dir = os.path.join(out_dir, 'blast_small_t')
    unaligned = os.path.join(blast_small_t_dir, 'unaligned.fa')
    aligned = os.path.join(blast_small_t_dir, 'aligned.fa')

    if os.path.isfile(unaligned):
        os.remove(unaligned)

    if os.path.isfile(aligned):
        os.remove(aligned)

    for key in known_small_t:
        with open(unaligned, 'a') as sequence_file:
            sequence_file.write(">{}\n".format(ID))
            sequence_file.write("{}\n".format(small_t_seq))
            sequence_file.write(">{}\n".format(key))
            sequence_file.write("{}\n".format(known_small_t[key]))

    cline = MuscleCommandline(input=unaligned, out=aligned, verbose=False)

    stdout, stderr = cline()

    return aligned
# ----------------------------------------------------------------------------------------
def locate_large_t(virus,splice_acceptor, args):
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
    large_t = {}
    aligned = align_splice_donor(virus, args)
    known_splice_donor, _ = locate_known_splice_donor(virus, args)
    if known_splice_donor < 1:
        known_splice_donor = -1
        large_t['Large_T'] = [known_splice_donor]
        return large_t
    genome = str(virus['genome'])
    align_seq = []
    j = 0
    aligned_splice_donor_stop = 0

    for aln in AlignIO.read(aligned, 'fasta'):
        align_seq.append(aln.seq)
    unknown_seq = str(align_seq[0]).lower()
    known_seq = str(align_seq[1]).lower()

    for position in known_seq:
        aligned_splice_donor_stop = aligned_splice_donor_stop + 1
        if position.lower() in ['a', 'c', 't', 'g']:
            j = j + 1
            if j == known_splice_donor:
                break
    search_seq = Seq(unknown_seq[aligned_splice_donor_stop:aligned_splice_donor_stop +
                                                           50].replace('-', '')).reverse_complement()
    splice_donor_nt = (re.search(str(search_seq), str(genome).lower()).end())
    start_large_t_1 = virus['small_t'][1]
    large_t_stop = virus['Large_T'][0]

    large_t_seq = str(Seq(str(genome[large_t_stop -1:splice_acceptor] +
                      genome[splice_donor_nt -1: start_large_t_1]).lower(
    )).reverse_complement())
    large_t_trans = Seq(large_t_seq).translate()
    large_t['Large_T'] = [splice_donor_nt, start_large_t_1, large_t_stop,
        splice_acceptor, large_t_seq, large_t_trans]
    return large_t
# ----------------------------------------------------------------------------------------
def print_genome_info(virus):
    """
    This function prints all the found annotations

    :param virus: dictionary that has found proteins VP1, VP2, small t and Large T
    :return: nothing
    """
    print("\nThis is the gene information for {}:".
        format(virus['accession']))
    for name in virus:
        if name == 'genome':
            pass
        elif name == 'accession':
            pass
        elif name == 'name':
            pass
        elif name == 'Large_T_orf_trans':
            pass
        elif name == 'Large_T_orf_seq':
            pass
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
    print("All annotations displayed above.\n")
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
        del virus_copy['Large_T_orf_trans']
        del virus_copy['Large_T_orf_seq']
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
            elif value == 'Large_T':
                out_file.writerows([[
                    virus_copy['accession'], value,
                    'join(' + str(virus_copy[value][0]) + '..' +
                    str(virus_copy[value][1]) + '+' +
                    str(virus_copy[value][2]) + ".." +
                    str(virus_copy[value][3]) + ')', virus_copy[value][4],
                    virus_copy[value][5], virus_copy[value][-1]
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
                    virus_copy[value][3], virus_copy[value][-1]
                ]])
    return
# ----------------------------------------------------------------------------------------
def validate_args(args):
    """
    Validates command line arguments

    :param args: command line arguments
    :return: dictionary, keys are the command argument names, values are the command
    line values
    """

    #sites = args['sites'] if 'sites' in args else ['ALL']
    #sites = args.get('sites') or []
    input_file = args.get('input')
    out_dir = args.get('out_dir')
    data_dir = args.get('data_dir')
    input_format = args.get('format').lower()
    min_prot_len = args.get('min_prot_len')
    e_value = args.get('evalue')
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    for_user_dir = os.path.join(out_dir, 'for_user')
    if not os.path.isdir(for_user_dir):
        os.makedirs(for_user_dir)
    program_dir = os.path.join(out_dir, 'program_files')
    if not os.path.isdir(program_dir):
        os.makedirs(program_dir)
    return {
        'input_file': input_file,
        'out_dir': out_dir,
        'data_dir': data_dir,
        'input_format': input_format,
        'min_prot_len': min_prot_len,
        'e_value': e_value,
        'for_user_dir': for_user_dir,
        'program_files_dir': program_dir
    }
# ----------------------------------------------------------------------------------------
def run(args):
    """main"""
    virus = {}
    logging.info('run = {}\n'.format(args))
    warnings.simplefilter('ignore', BiopythonWarning)
    args = validate_args(args)
    for seq_record in SeqIO.parse(args['input_file'], args['input_format']):
        original_genome = seq_record.seq
        name = seq_record.description
    virus['genome'] = original_genome
    full_name = name.split("|")[1]
    ID = name.split("|")[0]
    virus['name'] = full_name
    virus['accession'] = ID
    print("\n" + ID + "\n")
    found = identify_main_proteins(original_genome, args)
    virus.update(found)
    orfs = find_orfs_with_trans(original_genome, 1, args['min_prot_len'])

    if ('Large_T' in virus) and ('small_t' in virus):
        for orf in orfs:
            if str(found['Large_T'][3][:-1]) in orf:
                virus['Large_T_orf_trans'] = orf
                orf_len = len(orf)
                start = orfs[orf]
                virus['Large_T_orf_seq'] = str(virus['genome'][start : (start +3) +(
                        orf_len+1)*3].reverse_complement()).lower()

        splice_acceptor = find_splice_acceptor(virus,args)
        if splice_acceptor > 1:
            large_t = locate_large_t(virus,splice_acceptor,args)
            if len(large_t['Large_T']) > 1:
                strand = virus['Large_T'][-1]
                del virus['Large_T']
                virus['Large_T'] = large_t['Large_T']
                virus['Large_T'].append(strand)
            else:
                pass

    if 'small_t' not in virus:
        del virus['Large_T']
    print_genome_info(virus)
    to_csv(virus, args['for_user_dir'])
    return 1