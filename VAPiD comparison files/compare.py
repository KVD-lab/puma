import re, csv


v_accession= [] #store accession numbers from PuMA analysis
r_accession=[] #store accession numbers downloaded from PaVE


exempt = ["CG", "E6*", "E1BS", "E2BS", "URR"]#"E1^E4","E8^E2"] #features excluded from analysis



#parse PaVE data and obtain all unique accessions
for r in csv.reader(open("data_dir/PaVE.csv","r")): 
    if r[0] not in r_accession:
        r_accession.append (r[0])
 
#parse VAPiD data and obtain all unique accessions
for r in csv.reader(open("combined.csv","r")): 
    if r[0].split("|")[0] not in v_accession:
        v_accession.append (r[0].split("|")[0])
        





#count total features in PaVE

count = 0
 
for r in csv.reader(open("data_dir/PaVE.csv","r")):
    if r[1] not in exempt:
        count = count + 1
 
print ("total genes in PaVE: "+str(count))




PuMA_dict = {}
PaVE_dict = {}


#construct a dictionary that contains the accession, gene, and start..stop. Only use those that VAPiD processed
for a in v_accession:
    temp = {}
    for r in csv.reader(open("data_dir/PaVE.csv","r")):
        if a == r[0]:
            if r[1] not in exempt:
                temp[r[1]] = r[3].upper()
    PaVE_dict[a] = temp

#construct a dictionary that contains the accession, gene, and start..stop. Only use those that VAPiD processed
for a in v_accession:
    temp = {}
    for r in csv.reader(open("combined.csv","r")):
        if a == r[0].split("|")[0]:
            if r[1] not in exempt:
                temp[r[1]] = r[3].upper()
    PuMA_dict[a] = temp


#print (PuMA_dict)
#print (PaVE_dict)


wrong = []
missing=[]

for v in PuMA_dict:
    for k in PuMA_dict[v]:
        if k in PaVE_dict[v]:
            if PuMA_dict[v][k] != PaVE_dict[v][k]:
                wrong.append([v,k])
        else:
            wrong.append([v,k])
print ("wrong: "+str(len(wrong)))
print (wrong)

for v in PaVE_dict:
    for k in PaVE_dict[v]:
        if k not in PuMA_dict[v]:
            missing.append([v,k])
print ("missing: "+str(len(missing)))
print (missing)
        