# ----------------------------------------------------------------
# Input: 
# Output: CSV file of DMS reactivities, in SHAPE reactivity format
# Reactivity calculation:
#  1) Mutation rate at each position = # of mutations/# of reads
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
for rep in ['DMS1', 'DMS2', 'DMS3']:
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
    all_reacts = []

    for rep in np.arange(len(muts)):
        
        mut = muts[rep][i][20:290]
        cov = covs[rep][i][20:290]
        
        reacts = []
        for j in np.arange(len(seq)):
            if seq[j] in ["A","a","C","c"]:
                if not pass_threshold[j]:
                    reacts.append(np.nan)
                else:
                    raw_react = mut[j]/cov[j]
                    reacts.append(raw_react)
            else:
                reacts.append(np.nan)
        all_reacts.append(reacts)
    #print(all_reacts[0])

    # Average all reactivities across 3 replicates, ignoring nan values
    avg_reacts = np.nanmean(all_reacts,axis = 0)
    #print("Average reactivities across 3 replicates:"+str(avg_reacts))

    # Normalize data 0 to 1
    norm_reacts = np.ma.array(avg_reacts,mask=np.isnan(avg_reacts))
    norm_reacts = (norm_reacts - np.min(norm_reacts)) / (np.max(norm_reacts) - np.min(norm_reacts))

    # Convert existing nanas to -999 for SHAPE reactivity format
    norm_reacts[np.isnan(norm_reacts)] = -999
    #print("Final reactivity data: "+str(norm_reacts))
    
    name = pool_names[i]

   # print(react_dict)
    
    # Write positions and reactivities to a csv file (SHAPE reactivity format)
    pos = np.arange(1, len(seq)+1)
    c = open(path+f"/DMS_fa_files/{name}_DMS_FL.csv","w")
    writer = csv.writer(c, delimiter="\t")
    writer.writerows(zip(pos,norm_reacts))
    c.close()
    
    # Create fasta file with single sequence
    fa = open(path+f"/DMS_fa_files/{name}_FL.fa","w")
    fa.write(">"+name+"\n")
    fa.write(seq+"\n")
    fa.close()
    
    
    #break

timepoint2 = time.time()

print("Time taken to parse fasta:"+str(timepoint0-start))
print("Time taken to unpickle:"+str(timepoint1-timepoint0))
print("Time taken to loop all sequences:" +str(timepoint2-timepoint1))
print("low coverage sequence count:" + str(low_cov_seq_count))
        

