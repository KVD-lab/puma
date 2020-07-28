import csv
from Bio.Align.Applications import MuscleCommandline

wrong={}
for w in csv.reader(open("wrong.csv","r")):
    wrong[w[0]] = w[1]

for w in (wrong):
    with open(w+"_"+wrong[w]+".fas", "w") as out:
        for p in csv.reader(open("data_dir/PaVE.csv","r")):
            if w == p[0] and wrong[w]==p[1]:
                print (f">{w}_{p[1]}_PaVE\n{p[3]}", file=out)
        for c in csv.reader(open("combined.csv","r")):
            if w == c[0] and wrong[w]==c[1]:
                print (f">{w}_{c[1]}_PuMA\n{c[3]}",file=out)
    
    unaligned = w+"_"+wrong[w]+".fas"
    aligned = w+"_"+wrong[w]+".muscle.fas"
    

    cline = MuscleCommandline(input=unaligned, out=aligned, verbose=False)
    stdout, stderr = cline()