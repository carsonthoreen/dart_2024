#!/usr/bin/env python3

# Extract coverage from DMS .pkl files, count and calculate percent sequences that pass a given threshold provided in command.

import numpy as np
import pickle
import os
import sys

path = os.getcwd()

group = "noDMS"

covs = []
for rep in [f"{group}1",f"{group}2",f"{group}3"]:
    with open(path+"/"+rep+"__barcode_Bar2_1_mut_cov_PE.pkl", "rb") as f:
        mut, cov = pickle.load(f)
    covs.append(cov)

print(len(covs[0])) #should be total number of sequences (40191)
print(len(covs)) #should be total number of replicates (3)

low_cov_seq_count = 0
pass_count = 0
threshold = int(sys.argv[1])

for i in range(len(covs[0])): # loop through all sequences in the pool

    cover_3rep = []
    for rep in np.arange(len(covs)): #loop through replicates
        cover = covs[rep][i][20:290]
        cover_3rep.append(cover)
    
    avg_cov = np.average(cover_3rep)
    #print(avg_cov)
    #print(type(avg_cov))
    if avg_cov < threshold:
        low_cov_seq_count += 1
    else:
        pass_count += 1
    
    percentpass = pass_count/(pass_count+low_cov_seq_count)

print(f"Sequences pass coverage threshold threshold {threshold}: {pass_count}")
print(f"Sequences with low coverage: {low_cov_seq_count}")
print(f"Percent pass threshold: {percentpass:.2%}")



