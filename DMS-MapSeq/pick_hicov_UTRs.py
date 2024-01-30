#!/usr/bin/env python3

# DJ 11/17/2023 
# Input: .csv file from DEseq analysis for Nsp1 DART (a list of Nsp1-resistant UTRs). Check DMS MaP seq converage and write results to a new CSV file.

import csv
import pandas as pd
from Bio import SeqIO
import numpy as np
import pickle

# Read .csv file and extract sequence IDs
df = pd.read_csv("80S_spiked_lrt_filtered_pos_sig.csv", sep = ",")
colname = df.columns[1]
IDs = df[colname]
#print(colname) # Should be column name "ID"
#print(IDs) # Should be the ID column without header
print(IDs[0])
print(len(IDs[0]))

# Parse fasta file
pool_names=[]
with open("New_DART_reduce.fa","r") as fasta_h:
    for i in SeqIO.parse(fasta_h,"fasta"):
        pool_names.append(i.id)

# Extract coverage of each ID of interest
covs = [] # covs is the 3*40191 array that stores coverage of all seqs in the 3 replicates
for rep in ['DMS1', 'DMS2', 'DMS3']:
    with open(rep+"__barcode_Bar2_1_mut_cov_PE.pkl","rb") as f:
        mut, cov = pickle.load(f)
    covs.append(cov)
print(len(covs)) #should be 3

# Loop through all UTRs of interest, find index in parsed fasta file (same as in pkl file), then extract average and compare with threshold, store True or False
pass_threshold = []
Threshold = 400
for id in IDs:
    index = pool_names.index(id)
   # print(index)

    cover_3rep = []
    for rep in np.arange(len(covs)):
        cover = covs[rep][index][20:290]
        cover_3rep.append(cover)
    #print(np.shape(cover_3rep)) #Should be (3,270)
    #print(np.average(cover_3rep))
    
    pass_threshold.append(np.average(cover_3rep)>=Threshold)

print(pass_threshold)
pass_count = pass_threshold.count(True)
print(pass_count)

fw = open("80S_spiked_lrt_filtered_pos_sig_hicov.tsv","w")
writer = csv.writer(fw,delimiter = "\t")
writer.writerows(zip(IDs,pass_threshold))
fw.close()
