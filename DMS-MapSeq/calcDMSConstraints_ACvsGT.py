# ----------------------------------------------------------------
# Input: 
# Output: TSV file of DMS reactivities, in SHAPE reactivity format
# Reactivity calculation:
#  1) Mutation rate at each position = # of mutations/# of reads
#     AC and GT mutation rates are calculated and stored separately.
#  2) Average mutation rate across each replicate
#  3) Normalize data 0 to 1
#
#
# ----------------------------------------------------------------
#!/usr/bin/env python3

import numpy as np
import pickle
import sys
from Bio import SeqIO
import os
import time
import csv

start = time.time()

group = "noDMS" #denatured, DMS or noDMS
path = os.getcwd()
#print(path)

pool_names = []
seqs = []

with open(path+"/New_DART_reduce.fa","r") as fasta_h:
    for i in SeqIO.parse(fasta_h,"fasta"):
        pool_names.append(i.id)
        seqs.append(str(i.seq))

#print(len(seqs))
#print(seqs[1])

timepoint0 = time.time()
muts = []
covs = []
for rep in [f'{group}1', f'{group}2', f'{group}3']:
    with open(path+"/"+rep+"__barcode_Bar2_1_mut_cov_PE.pkl","rb") as f:
        mut, cov = pickle.load(f)
    muts.append(mut)
    covs.append(cov)

print(len(muts[0])) # should be 40191
print(len(covs[0]))
print(len(muts)) #should be 3
print(len(covs)) #should be 3

timepoint1 = time.time()

low_cov_seq_count = 0
for i in range(len(seqs)): # Loop through all sequences in the pool

    threshold = 400
    seq = seqs[i][20:290]
    pass_threshold = [True]*len(seq)
    #print(seq)

    # Loop through replicates, create array of positions based on whether position passes cov threshold in all reps
    cover_3rep = []
    for rep in np.arange(len(covs)):
        cover = covs[rep][i][20:290]
        cover_3rep.append(cover)

        for j in np.arange(len(seq)):
            if cover[j] < threshold:
                pass_threshold[j] = False
   # print(pass_threshold)
    #print(cover_3rep)
    if np.average(cover_3rep) < threshold:
        low_cov_seq_count += 1
        continue

    # Loop through replicates, calculate raw reactivities for all replicates, setting low coverage and G/T positions as nan
    all_reactsAC = []
    all_reactsGT = []

    for rep in np.arange(len(muts)):
        
        mut = muts[rep][i][20:290]
        cov = covs[rep][i][20:290]
        
        reactsAC = []
        reactsGT = []
        for j in np.arange(len(seq)):
            if seq[j] in ["A","a","C","c"]:
                reactsGT.append(np.nan)
                if not pass_threshold[j]:
                    reactsAC.append(np.nan)
                else:
                    raw_react = mut[j]/cov[j]
                    reactsAC.append(raw_react)
            else:
                reactsAC.append(np.nan)
                if not pass_threshold[j]:
                    reactsGT.append(np.nan)
                else:
                    raw_react = mut[j]/cov[j]
                    reactsGT.append(raw_react)
        all_reactsAC.append(reactsAC)
        all_reactsGT.append(reactsGT)
    #print(all_reactsAC[0])
    #print(all_reactsGT[0])

    # Average all reactivities across 3 replicates, ignoring nan values
    avg_reactsAC = np.nanmean(all_reactsAC,axis = 0)
    avg_reactsGT = np.nanmean(all_reactsGT,axis = 0)
    #print("Average AC reactivities across 3 replicates:"+str(avg_reactsAC))
    #print("Average GT reactivities across 3 replicates:"+str(avg_reactsGT))

    # Normalize data 0 to 1
    #norm_reacts = np.ma.array(avg_reacts,mask=np.isnan(avg_reacts))
    #norm_reacts = (norm_reacts - np.min(norm_reacts)) / (np.max(norm_reacts) - np.min(norm_reacts))

    # Convert existing nanas to -999 for SHAPE reactivity format
    #norm_reacts[np.isnan(norm_reacts)] = -999
    #print("Final reactivity data: "+str(norm_reacts))
    
    name = pool_names[i]

   # print(react_dict)
    
    # Write positions and reactivities to a csv file (SHAPE reactivity format)
    namecol = [name]*len(seq)
    pos = np.arange(1, len(seq)+1)

    c = open(path+f"/raw_reacts/{group}_rawreactsAC.tsv","a")
    writer = csv.writer(c, delimiter="\t")
    writer.writerows(zip(namecol,pos,avg_reactsAC))
    c.close()
    
    d = open(path+f"/raw_reacts/{group}_rawreactsGT.tsv","a")
    writer = csv.writer(d, delimiter="\t")
    writer.writerows(zip(namecol,pos,avg_reactsGT))
    d.close()

    # Create fasta file with single sequence
    #fa = open(path+f"/DMS_fa_files/{name}_FL.fa","w")
    #fa.write(">"+name+"\n")
    #fa.write(seq+"\n")
    #fa.close()
    
    
    #break

timepoint2 = time.time()

print("Time taken to parse fasta:"+str(timepoint0-start))
print("Time taken to unpickle:"+str(timepoint1-timepoint0))
print("Time taken to loop all sequences:" +str(timepoint2-timepoint1))
print("low coverage sequence count:" + str(low_cov_seq_count))
        

