#!/usr/bin/env python3
# Josh Pace, Ken Youens-Clark, Cordell Freeman, Koenraad Van Doorslaer
#University of Arizona, KVD Lab & Hurwitz Lab
# PuMA 0.3-beta 1/29/2019

from distutils.spawn import find_executable
from Bio import SeqIO, GenBank, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio.Blast.Applications import NcbiblastnCommandline as blastn
from Bio.Align.Applications import MuscleCommandline
import os, glob, re, csv, time, operator, argparse, sys
import warnings
from Bio import BiopythonWarning
import itertools
import shutil

# ----------------------------------------------------------------------------------------
#
#This functions finds all the open reading frames and translates them
#
def trans_orf(seq, trans_table, min_protein_length):
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

# ----------------------------------------------------------------------------------------
#
# This function uses blast to find the different proteins in
# the translated genome
#

def blast_proteins(genome,min_prot_len,evalue,data_dir, out_dir):
    protein_start = {}
    protein_seq = {}
    found_proteins = {}

    orfs = trans_orf(genome, 1, min_prot_len)

    if not orfs:
        print('No ORFs, must stop.')
        sys.exit(1)

    orfs_fa = os.path.join(out_dir, 'orfs.fa')
    orfs_fh = open(orfs_fa, 'wt')

    for orf in orfs:
        orfs_fh.write('\n'.join(['>' + str(orfs[orf]), orf, '']))
    orfs_fh.close()


    blast_sub = os.path.join(data_dir, 'main_blast.fa')
    blast_out = os.path.join(out_dir, 'blast_results.tab')

    if os.path.isfile(blast_out):
        os.remove(blast_out)


    cmd = blastp(query=orfs_fa,
                 subject=blast_sub,
                 evalue=evalue,
                 outfmt=6,
                 out=blast_out)

    stdout, stderr = cmd()
    if stdout:
        print("STDOUT = ", stdout)
    if stderr:
        print("STDERR = ", stderr)

    if not os.path.isfile(blast_out) or not os.path.getsize(blast_out):
        print('No BLAST output "{}" (must have failed)'.format(blast_out))
        sys.exit(1)


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
                    M = re.search('M', protein_seq[seq])
                    real_start = protein_start[start] + M.start() + M.start() + M.start()
                    end = protein_start[start] + ((len(protein_seq[seq])+1) * 3)
                    found_proteins['L1'] = []
                    L1_pre = genome[(protein_start[start] + 3 * M.start()):int(end)]
                    splice = '(C|T)(C|T)(A|C|G|T)(C|T)AG(A)TG'
                    spliced = re.search(splice, str(L1_pre))
                    if spliced:
                        start_L1 = int(spliced.start()) + 6
                        if start_L1 % 3 == 0:
                            if start_L1 > 600:
                                L1_post = L1_pre
                                found_proteins['L1'] = [int(start_L1),
                                    int(end),str(L1_post).lower(),Seq(str(
                                        L1_post)).translate()]
                            else:
                                L1_post = L1_pre[start_L1:]
                                found_proteins['L1'] = [int(real_start) + 1 +
                                       int(start_L1),int(end), str(L1_post).lower(),
                                    Seq(str(L1_post)).translate()]
                        else:
                            L1_post = L1_pre
                            found_proteins['L1'] = [int(real_start) + 1, int(end),
                                str(L1_post).lower(),
                            Seq(str(L1_post)).translate()]
                    else:
                        L1_post = L1_pre
                        found_proteins['L1'] = [int(real_start) + 1,
                            int(end), str(L1_post).lower(), Seq(str(L1_post)).translate()]
                else:
                    try:
                        M = re.search('M', protein_seq[seq])
                        real_start = protein_start[start] + M.start() + M.start() + M.start()
                        end = protein_start[start] + ((len(protein_seq[seq])+1) * 3)
                        sequence = str(genome[int(real_start):int(end)]).lower()
                        translated = Seq(sequence).translate()
                        found_proteins[seq] = [int(real_start) + 1, int(end), sequence,
                        translated]
                    except AttributeError:
                        pass




    return found_proteins

# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
#
#Finds the E1 binding site in the genome using the URR
#

def find_E1BS(genome, URR, URRstart, ID,data_dir, out_dir):
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
    cline = (fimo_cmd.format(fimo_exe, fimo_dir, background ,
                             motif, tmp))

    os.system(str(cline))

    fimo_out = os.path.join(fimo_dir, 'fimo.tsv')

    if not os.path.isfile(fimo_out):
        print('Failed to create fimo out "{}"'.format(fimo_out))
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
        sequence = str(genome[int(genomestart)-2:] + genome[:genomestop]).lower()
        E1BS['E1BS'] = [int(genomestart), int(genomeLength), 1, int(genomestop),sequence]

    else:
        if genomestart == 1:
            sequence = str(genome[-1]).lower() + str(genome[int(genomestart - 1):int(
                genomestop)]).lower()
            E1BS['E1BS'] = [int(genomestart), int(genomestop), sequence]

        else:
            sequence = str(genome[int(genomestart-2):int(genomestop)]).lower()
            E1BS['E1BS'] = [int(genomestart), int(genomestop),sequence]

    return E1BS
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
#
# This function finds the E2BS in a genome using the URR
#
def find_E2BS(genome, URR, URRstart, ID, data_dir, out_dir):
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
        print('Please install "fimo" into your $PATH')
        return

    fimo_dir = os.path.join(out_dir, 'E2BS')
    if not os.path.isdir(fimo_dir):
        os.makedirs(fimo_dir)

    motif = os.path.join(data_dir,'E2BS_motif.txt' )

    fimo_cmd = '{} --oc {} --norc --verbosity 1 --thresh 1.0E-4 {} {}'
    cline = (fimo_cmd.format(fimo_exe, fimo_dir,motif , tmp))

    os.system(str(cline))

    fimo_out = os.path.join(fimo_dir, 'fimo.tsv')

    if not os.path.isfile(fimo_out):
        print('Failed to create fimo out "{}"'.format(fimo_out))
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
                    startListGenome.append(genomestart-1)
                else:
                    startListGenome.append(genomestart-1)

            E2BS['E2BS'] = startListGenome

            return E2BS

        except KeyError:

            E2BS['E2BS'] = ["No E2BS found"]
            return E2BS



# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
#
#Finds the E4, used for stop site for E1^E4
#
def find_E4(E2, genome):
    E4 = {}

    trans_E2 = Seq(E2[1:len(E2)]).translate()
    E4protein_long = str(max(trans_E2.split("*"),key=len))
    if 'M' in E4protein_long:
        M = re.search('M', E4protein_long)
        if M.start() > 41:
            E4protein = E4protein_long
        else:
            E4protein_tmp = E4protein_long.split('M',1)[1]
            E4protein = 'M' + str(E4protein_tmp)
    else:
        E4protein = E4protein_long

    E4_start = re.search(str(E4protein),str(trans_E2)).start()
    E4_end = re.search(str(E4protein), str(trans_E2)).end()
    E4_nt = str(E2[(E4_start * 3) + 1:((E4_end + 1) * 3) + 1])
    E4_nt_start = re.search(E4_nt, str(genome)).start()
    E4_nt_end = E4_nt_start + len(E4_nt)
    sequence = str(genome[int(E4_nt_start):int(E4_nt_end)]).lower()
    translated = Seq(sequence).translate()
    E4['E4'] = [int(E4_nt_start) + 1, int(E4_nt_end) + 1]
    return E4
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
#
#Finds the blast results for E1 and E8
#
def find_blastResults_E1E8(E1_whole,ID, out_dir,data_dir):
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

    cmd = blastp(query=query_file,
                 subject=blast_subject,
                 evalue=1e-10,
                 outfmt=6,
                 out=blast_out)

    stdout, stderr = cmd()

    if not os.path.isfile(blast_out) or not os.path.getsize(blast_out):
        print('No BLAST output "{}" (must have failed)'.format(blast_out))
        blastResult = False
        return blastResult

    else:
        blastResult = True
    return blastResult
# ----------------------------------------------------------------------------------------
#
#Finds E1^E4
#
def find_E1E4(E1_whole,E2_whole,ID,genome,start_E4_nt,blastE1E8_dir,blastResult, data_dir):
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
        E1_E4['E1^E4'] = [0, 0, 0, 0, 'No Blast Output','']
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
        csv_database = os.path.join(data_dir, 'all_pave.csv')

        with open(csv_database, 'r') as csvfile:
            read = csv.DictReader(csvfile, ('accession', 'gene', 'positions', 'seq'))
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

    aligned_E1_stop = E1_stop[-1]
    search_seq = unknown_seq[aligned_E1_stop:aligned_E1_stop + 50].replace('-', '')

    stopE1_nt = (re.search(search_seq, str(genome).lower()).start()) + 1
    startE1_nt = E1_whole[0]

    whole_E4 = find_E4(E2_seq, genome)
    stopE4_nt = whole_E4['E4'][1] - 1

    E1_E4_seq = str(genome[startE1_nt-1:stopE1_nt]+ genome[
    start_E4_nt:stopE4_nt])
    E1_E4_trans = Seq(E1_E4_seq).translate()[:-1]
    E1_E4['E1^E4'] = [startE1_nt,stopE1_nt,start_E4_nt + 1,stopE4_nt,E1_E4_seq,
        E1_E4_trans]
    return E1_E4
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
#
#Finds E8^E2
#
def find_E8E2(E1_whole, E2_whole, ID, genome,startE2_nt,blastE1E8_dir, blastResult, data_dir):
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
        csv_database = os.path.join(data_dir, 'all_pave.csv')

        with open(csv_database, 'r') as csvfile:
            read = csv.DictReader(csvfile, ('accession', 'gene', 'positions', 'seq'))
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
            #print(E8_stop_genome)


            E8_stop_genome = int(str(E8_stop_genome).split('..')[1])

            known_E8_stop = int(known_E8_stop.split('..')[0])

            E8_stop_known = (E8_stop_genome - known_E8_stop)

            splice_sites.append(E8_stop_known)

        except UnboundLocalError:
            stopE8_nt = False
            return stopE8_nt

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

    search_seq = unknown_seq[aligned_E8_stop:aligned_E8_stop + 50].replace('-', '')

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

    if check_seq.startswith(startMotif):
        pass
    else:
        startE8_nt = startE8List[0] + E1_whole[0]

    stopE8_nt = (stopE8 + E1_whole[0]) + 1

    if startE2_nt != False:
        stopE2_nt = E2_whole[1]

    E8_E2_seq = Seq(genome[startE8_nt - 1:stopE8_nt] + genome[
    startE2_nt:stopE2_nt])
    E8_E2_trans = E8_E2_seq.translate()


    E8_E2['E8^E2'] = [startE8_nt, stopE8_nt, startE2_nt + 1, stopE2_nt, E8_E2_seq,
            E8_E2_trans]

    return E8_E2
# ----------------------------------------------------------------------------------------
#
# ----------------------------------------------------------------------------------------
#
#Finds splice acceptor posistion for E1^E4 and E8^E2
#
def find_splice_acceptor( E2_whole, ID,genome, data_dir, out_dir):

    E2_trans = str(E2_whole[3])
    E2_seq = E2_whole[2]

    splice_acceptor_dir = os.path.join(out_dir, 'splice_acceptor')
    if not os.path.isdir(splice_acceptor_dir):
        os.makedirs(splice_acceptor_dir)

    blast_subject = os.path.join(data_dir, 'splice_acceptor_blast.fa')
    blast_out = os.path.join(splice_acceptor_dir, 'blast_result.tab')

    if os.path.isfile(blast_out):
        os.remove(blast_out)

    query_file = os.path.join(splice_acceptor_dir, 'query.fa')

    with open(query_file, 'a') as query:
        query.write('>{}\n'.format(ID))
        query.write(E2_trans)

    cmd = blastp(query=query_file,
                 subject=blast_subject,
                 evalue=1e-10,
                 outfmt=6,
                 out=blast_out)

    stdout, stderr = cmd()

    if not os.path.isfile(blast_out) or not os.path.getsize(blast_out):
        print('No BLAST output "{}" (must have failed)'.format(blast_out))
        startE2_nt = "No blast results for unkown E2"
        return startE2_nt

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
        csv_database = os.path.join(data_dir, 'all_pave.csv')


        with open(csv_database, 'r') as csvfile:
            read = csv.DictReader(csvfile,('accession', 'gene', 'positions', 'seq'))
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
            if position.lower() in ['a','c','t','g']:
                j = j + 1
                if j == splice_start_known:
                    break

        aligned_starts.append(aligned_splice_start)

    aligned_splice_start = aligned_starts[-1]


    search_seq = unknown_seq[aligned_splice_start:aligned_splice_start + 50].replace('-','')

    startE2_nt = re.search(search_seq,str(genome).lower()).start()


    return startE2_nt

# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
#
#Outputs gff3 format
#
def to_gff3(dict, genomelen, out_dir):
    del dict['genome']
    del dict['E1BS']
    del dict['E2BS']
    del dict['URR']

    name = dict['accession']
    dict['name'] = name
    gff3_out = os.path.join(out_dir, '{}.gff3'.format(dict['name']))

    with open(gff3_out, 'a') as out_file:
        out_file.write("##gff-version 3\n")
        out_file.write("##sequence-region {} 1 {}\n".format(dict['name'], genomelen))

    for protein in dict:
        if protein == 'name' or protein == 'accession':
            pass
        elif "^" in protein:
            with open(gff3_out, 'a') as out_file:
                out_file.write(
            "{}\tPuMA\tCDS\t{}\t{}\t{}\t{}\t.\t+\t.\tID={};Note=[{}-{} + {}-{}]\n".format(
                 dict['name'], dict[protein][0], dict[protein][1],
                dict[protein][2],dict[protein][3],protein, dict[protein][0],
                dict[protein][1],dict[protein][2],dict[protein][3]))
        else:
            with open(gff3_out, 'a') as out_file:
                out_file.write(
            "{}\tPuMA\tsplice_site\t{}\t{}\t.\t+\t.\tID={};Note=[{}-{}]\n".format(
                 dict['name'], dict[protein][0], dict[protein][1],
                protein, dict[protein][ 0], dict[protein][1]))

    return
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
#
#Output to csv file
#
def export_to_csv(annotations,out_dir):

    csv_out = os.path.join(out_dir, 'puma_results.csv')

    with open(csv_out,'a') as out:
        out_file = csv.writer(out)

        for value in annotations:
            if value == 'genome':
                out_file.writerows([[annotations['accession'],
                                   'CG',"",str(annotations[value]).lower()]])
            elif value == 'URR':
                try:
                    out_file.writerows([[annotations['accession'],
                                       value, 'join(' + str(annotations[value][0])
                                       + '..' +
                                      str(annotations[value][1]) + '+' +
                                      str(annotations[value][2]) + ".." +
                                      str(annotations[value][3]) +')',
                                      annotations[value][4]]])
                except IndexError:
                    out_file.writerows([[annotations['accession'],
                                       value, str(annotations[value][0]) + '..' +
                                       str(annotations[value][1]),
                                        annotations[value][2]]])
            elif value == 'E1BS':
                try:
                    out_file.writerows([[annotations['accession'],
                                       value, 'join(' + str(annotations[value][0]) + '..'
                                       + str(annotations[value][1]) + '+' +
                                       str(annotations[value][2]) + ".." +
                                       str(annotations[value][3]) + ')',
                                       annotations[value][4]]])
                except IndexError:
                    out_file.writerows([[annotations['accession'],
                                        value, str(annotations[value][0]) + '..' +
                                        str(annotations[value][1]),
                                        annotations[value][2]]])
            elif value == 'E2BS':
                for i in range(0,len(annotations[value]),1):
                                  out_file.writerows([[annotations['accession'], value,
                                   str(annotations[value][i]) + '..' +
                                   str(annotations[value][i] +11),str(annotations['genome']
                                   [annotations[value][i]-1:annotations[value][
                                    i]+11]).lower()]])
            elif value == 'E1^E4':
                out_file.writerows([[annotations['accession'], value,
                                    'join(' + str(annotations[value][0]) + '..' +
                                     str(annotations[value][1]) + '+' +
                                    str(annotations[value][2]) + ".."+
                                    str(annotations[value][3]) + ')',
                                    annotations[value][4],annotations[value][5]]])
            elif value == 'E8^E2':
                out_file.writerows([[annotations['accession'], value,
                                      'join(' + str(annotations[value][0]) + '..' +
                                      str(annotations[value][1]) + '+' +
                                      str(annotations[value][2]) + ".."
                                      + str(annotations[value][3]) + ')'
                                      ,annotations[value][4],annotations[value][5]]])
            elif value == 'name':
                pass
            elif value == 'accession':
                pass
            else:
                out_file.writerows([[annotations['accession'],
                                    value, str(annotations[value][0]) + '..' +
                                    str(annotations[value][1]),
                                     annotations[value][2], annotations[value][3]]])

    return
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
def get_args():
    args = sys.argv
    bin = os.path.dirname(args[0])
    parser = argparse.ArgumentParser(
        description='Displays identified protein '
        'information within a given papillomavirus '
        'genome.')

    parser.add_argument(
        '-i',
        '--input',
        metavar='FILE',
        help='Path to a  fasta formatted file that ' +
        'contains a papillomavirus genome.',
        required=True)

    parser.add_argument(
        '-f',
        '--format',
        metavar='FORMAT',
        default='fasta',
        help='FASTA')

    parser.add_argument(
        '-d',
        '--data_dir',
        metavar='DIR',
        type=str,
        default=os.path.join(bin, 'data_dir'),
        help='Directory that has all database files for use')

    parser.add_argument(
        '-o',
        '--outdir',
        metavar='DIR',
        type=str,
        default=os.path.join(bin, 'puma-out'),
        help='Output directory (default: %(default)s)')

    parser.add_argument(
        '-g',
        '--gff3',
        metavar='STR',
        type=str,
        default='y',
        help='Outputs information in gff3 format (y|n)')
    parser.add_argument(
        '-c',
        '--csv',
        metavar='STR',
        type=str,
        default='y',
        help='Outputs information in csv format (y|n)')
    parser.add_argument(
        '-e',
        '--evalue',
        metavar='FLOAT',
        type=float,
        default=0.00001,
        help='BLAST evalue')

    parser.add_argument(
        '-m',
        '--min_prot_len',
        metavar='NUM',
        type=int,
        default=25,
        help='Minimum protein length')

    parser.add_argument(
        '-s',
        '--sites',
        metavar='STR',
        type=str,
        default='ALL',
        help='Comma-separated string of L1 L2 E1 E2 E4 E5 E5_delta '
        'E5_zeta E5_epsilon E6 E7 '
        'E9 E10 '
        'E2BS E1BS URR ALL')

    return parser.parse_known_args()


# -----------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------
def main():
    """main"""
    warnings.simplefilter('ignore', BiopythonWarning)

    args, _ = get_args()
    sites = re.split(r'\s*,\s*', args.sites.upper())
    input_file = args.input
    out_dir = args.outdir
    data_dir = args.data_dir
    input_format = args.format.lower()
    min_prot_len = args.min_prot_len
    evalue = args.evalue
    gff3 = args.gff3
    csv = args.csv
    blastE1E8_dir = os.path.join(out_dir, 'blastE1E8')
    if not os.path.isdir(blastE1E8_dir):
        os.makedirs(blastE1E8_dir)


    valid_sites = set(
        'L1 L2 E1 E2 E4 E5 E5_delta E5_zeta E5_epsilon E6 E7 E9 E10 E2BS '
        'E1BS '
        'URR '
        'ALL'.split())
    if not sites:
        print('--sites is required')
        sys.exit(1)

    if not input_file:
        print("--i is required")
        sys.exit(1)

    if not os.path.isfile(input_file):
        print('--input "{}" is not a file.'.format(input_file))
        sys.exit(1)

    if not os.path.isdir(data_dir):
        print('--data_dir "{}" is not a directory.'.format(data_dir))
        sys.exit(1)

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    bad_sites = list(filter(lambda s: s not in valid_sites, sites))

    if bad_sites:
        print('Invalid site(s): {}'.format(', '.join(bad_sites)))
        sys.exit(1)

    valid_format = set(['fasta', 'genbank'])

    if not input_format in valid_format:
        msg = 'Invalid format ({}), please choose from {}'
        print(msg.format(input_format, ', '.join(valid_format)))
        sys.exit(1)

    startStop = []
    URR = {}
    virus = {}

    blasted = {}
    for seq_record in SeqIO.parse("{}".format(input_file), input_format):
        Origseq = seq_record.seq
        name = seq_record.description

    full_name = name.split("|")[1]
    ID = name.split("|")[0]
    virus['name'] = full_name
    virus['accession'] = ID
    virus['genome'] = Origseq

    print("\nThis is the gene information for {}:".format(virus['name']))

    blasted.update(
        blast_proteins(Origseq, min_prot_len, evalue, data_dir, out_dir))

    virus.update(blasted)
    for protein in blasted:
        if protein == 'L1':
            startStop.append((blasted[protein][1]))
            URRstart = blasted[protein][1]
        else:
            startStop.append(blasted[protein][0])

    startStop = sorted(startStop)

    for numbers in startStop:
        if numbers == URRstart:
            if numbers == startStop[-1]:
                URRstop = startStop[0]
            else:
                position = startStop.index(numbers)
                URRstop = startStop[position + 1]

    URRstart = int(URRstart) + 1
    URRstop = int(URRstop) - 1

    genomelen = int(len(Origseq))


    if URRstop == 0:
        URRstop = genomelen
    if URRstop > URRstart:
        URRfound = str(Origseq[URRstart - 1:URRstop]).lower()
    else:
        URRfound = str(Origseq[URRstart - 1:] + Origseq[:URRstop]).lower()
    if URRstop > URRstart:
        URR['URR'] = [int(URRstart), int(URRstop), URRfound]
    else:
        URR['URR'] = [int(URRstart), int(genomelen), 1, int(URRstop), URRfound]
    virus.update(URR)


    E2BS = find_E2BS(Origseq, URRfound, URRstart, ID, data_dir,out_dir)
    virus.update(E2BS)

    E1BS = find_E1BS(Origseq, URRfound, URRstart, ID,data_dir,
                     out_dir)

    virus.update(E1BS)

    start_splice_site = find_splice_acceptor(virus['E2'], ID, Origseq,
                                             data_dir, out_dir)
    blastResult = find_blastResults_E1E8(virus['E1'], ID, out_dir, data_dir)

    E1_E4 = find_E1E4(virus['E1'], virus['E2'], ID, Origseq, start_splice_site,
                      blastE1E8_dir,blastResult, data_dir)

    E8_E2 = find_E8E2(virus['E1'], virus['E2'], ID, Origseq, start_splice_site,
                      blastE1E8_dir,blastResult,data_dir)

    if E8_E2['E8^E2'] == False:
        pass
    else:
        virus.update(E8_E2)
    if E1_E4['E1^E4'] == False:
        pass
    else:
        virus.update(E1_E4)



    if csv == 'y':
        export_to_csv(virus,out_dir)
    else:
        pass

    if sites[0] == 'ALL':
        sites = {}
        sites.update(virus)
        del sites['name']
        del sites['genome']
        del sites['accession']

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
                print('Line 275:{}'.format(virus[name]))
                print('\n{} start and stop position:\n{},{}\n'.format(
                    name, virus[name][0], virus[name][1]))
                print('{} sequence:\n{}\n'.format(name, virus[name][2]))
                if name != 'URR':
                    print('{} translated seqeunce:\n{}\n'.format(
                        name, virus[name][3][:-1]))


    if gff3 == 'y':
        to_gff3(virus,genomelen,out_dir)
    else:
        pass

    return


# -----------------------------------------------------------------------------------------

if __name__ == '__main__':
    main()

