"""
File to compare results of annnoations with  PuMA and VAPiD 
KVD Lab, University of Arizona
Author: Dr. Koenraad Van Doorslaer
"""



import re, csv


v_accession= [] #store accession numbers from PuMA analysis
r_accession=[] #store accession numbers downloaded from PaVE


exempt = ["CG", "E6*", "E1BS", "E2BS", "URR","E1^E4","E8^E2"] #features excluded from analysis



#parse PaVE data and obtain all unique accessions
for r in csv.reader(open("PaVE.csv","r")): 
    if r[0] not in r_accession:
        r_accession.append (r[0])
 
#parse VAPiD data and obtain all unique accessions
for r in csv.reader(open("vapid_results.csv","r")): 
    if r[0].split("|")[0] not in v_accession:
        v_accession.append (r[0].split("|")[0])
        





#count total features in PaVE

count = 0
 
for r in csv.reader(open("PaVE.csv","r")):
    if r[1] not in exempt:
        count = count + 1
 
print ("total genes in PaVE: "+str(count))




vapid_dict = {}
PaVE_dict = {}


#construct a dictionary that contains the accession, gene, and start..stop. Only use those that VAPiD processed
for a in v_accession:
    temp = {}
    for r in csv.reader(open("PaVE.csv","r")):
        if a == r[0]:
            if r[1] not in exempt:
                temp[r[1]] = r[2]
    PaVE_dict[a] = temp

#construct a dictionary that contains the accession, gene, and start..stop. Only use those that VAPiD processed
for a in v_accession:
    temp = {}
    for r in csv.reader(open("vapid_results.csv","r")):
        if a == r[0].split("|")[0]:
            if r[1] not in exempt:
                temp[r[1]] = r[2]+".."+r[3]
    vapid_dict[a] = temp





wrong = []
missing=[]

for v in vapid_dict:
    for k in vapid_dict[v]:
        if k in PaVE_dict[v]:
            if vapid_dict[v][k] != PaVE_dict[v][k]:
                wrong.append([v,k])
        else:
            wrong.append([v,k])
print ("wrong: "+str(len(wrong)))
#print (wrong)

for v in PaVE_dict:
    for k in PaVE_dict[v]:
        if k not in vapid_dict[v]:
            missing.append([v,k])
print ("missing: "+str(len(missing)))
#print (missing)
        
