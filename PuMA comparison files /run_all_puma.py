"""
File to run PuMA agaisnt 655 genomes. This file helps create databases without the input sequence for more accurate and true results.
KVD Lab, University of Arizona
Author: Dr. Koenraad Van Doorslaer
"""


import csv,re,os
from Bio.Seq import Seq

accessions=[]
for r in csv.reader(open("data_dir/PaVE.csv","r")):
    if (r[0]) not in accessions:
        accessions.append(r[0])

#accessions=["KT626573","KX954132","MH049343","NC_001523","NC_001524","NC_001676","NC_001789","NC_011051","NC_001523","NC_001524","JQ692938","KC858263","MH049343","NC_024300"]

#accessions = ["NC_006563"]


with open("combined.csv","w") as out:
    for a in accessions:
        os.system("rm -r puma_out")
        with open("temp.fas","w") as temp, open("./data_dir/blast_E6_updated.fa","w") as E6, open("./data_dir/E1E8_blast.fa","w") as E1E8, open("./data_dir/splice_acceptor_blast_new.fa","w") as splice, open("./data_dir/blast_subject_all.fa","w") as subject:
        #with open("temp.fas","w") as temp, open("data_dir/blast_E6_updated.fa","w") as E6, open("data_dir/E1E8_blast.fa","w") as E1E8, open("data_dir/splice_acceptor_blast_new.fa","w") as splice:

            for r in csv.reader(open("data_dir/PaVE.csv","r")):
                if r[0] == a:
                    if r[1] == 'CG':
                        print (f">{r[0]}|{r[0]}\n{r[3]}", file = temp)
                else:
                    if r[1] == 'E6':
                        print (f">{r[0]}\n{str(Seq(r[3]).translate())}", file = E6)
                    elif r[1] == 'E1':
                        print (f">{r[0]}\n{str(Seq(r[3]).translate())}", file = E1E8)
                    elif r[1] ==  'E2':
                        print (f">{r[0]}\n{str(Seq(r[3]).translate())}", file = splice)
                    if r[4] != "":
                        print (f">{r[0]}\n{r[4]}", file = subject)
                        
                        
                        
        os.system("python scripts/run_puma.py -i temp.fas -d data_dir")
        for r in csv.reader(open("./scripts/puma_out/"+a+"/for_user/"+a+".csv","r")):
             print (",".join(r),file=out)
