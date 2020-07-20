"""
File to run VAPiD for PuMA vs VAPiD comparison

KVD Lab, University of Arizona

Author: Dr. Koenraad Van Doorslaer
"""



import os, csv, re
from Bio import SeqIO

#get accession numbers for the analysis
accessions=[]
for r in csv.reader(open('PaVE.csv',"r")):
    if r[0] not in accessions:
        accessions.append(r[0])



with open("vapid_results.csv","w") as out:
    for a in accessions:
        
        #write mock metdata file
        with open("example_metadata.csv","w") as meta:
            print (f"strain,collection-date,country,coverage\n{a},2017,USA,42.5",file=meta)
        #create blast database based on PaVE without the genome (a) under investigation
        with open("temp.fas","w") as db:
            for r in SeqIO.parse("PaVE.gb","genbank"):
                    if r.id != a:
                        print (f">{r.id}\n{r.seq}", file=db)
                    else:
                        #write the genome to be annotated to a fasta formatted file
                        print (a, r.id)
                        with open("example.fasta","w") as t:
                            print (f">{r.id}\n{r.seq}", file=t)
                
        #make the blast db
        os.system("makeblastdb -in temp.fas -dbtype nucl")
        #make the folder to store the data
        os.system("mkdir "+a)
        #execute VAPiD
        os.system("python2.7 vapid.py example.fasta example.sbt --metadata_loc example_metadata.csv --db temp.fas >"+a+"/screen.out")
        
        #parse VAPiD output
        if os.path.isfile(a+"/screen.out"):
            with open(a+"/screen.out","r") as f:
                lines=[]
                for l in f:
                    l=l.strip()
                    m= re.search("(.*) \['(\d*)', '(\d*)'\]",l)
                    if m:
                        g,s,e = m.groups()[0],str(int(m.groups()[1])),m.groups()[2]
                        print (",".join([a,g,s,e]),file=out)
                        print (",".join([a,g,s,e]))

                #     if "Feature" not in l and "codon" not in l:
                #         lines.append (l)
                # for line in (lines[10:-5]):
                #     #print (l)
                #     m = re.search("(.*) \['(.*)', '(.*)'\]",line)
                #     if m:
                #         gene = m.groups()[0]
                #         start = (m.groups()[1])
                #         end = (m.groups()[2])
                #     print (",".join([a, gene, start, end]), file=out)

        os.system("rm -rf "+a)
        os.system("rm temp.*")
        os.system("rm example_*")
