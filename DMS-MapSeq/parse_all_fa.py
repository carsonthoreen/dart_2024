#!/usr/bin/env python3

# ----------------------------------------------------------------
# file name: parse_all_fa.py
# Description: take a multi-fasta file (reference sequences for DMS-MaP-seq), trim first 20 and last 10 nt (contains T7 and adaptor sequences), then creat individual fasta files for every sequence. These fasta files will be used for RNAfold secondary structure prediction without DMS reactivity constraints.
#
# ----------------------------------------------------------------
# Required modules
# ------------------

from Bio import SeqIO
import os

path = os.getcwd()
#print(path)

pool_names = []
seqs = []

with open(path+"/New_DART_reduce.fa","r") as fasta_h:
    for i in SeqIO.parse(fasta_h,"fasta"):
        pool_names.append(i.id)
        seqs.append(str(i.seq))

for i in range(len(seqs)): # Loop through all sequences in the pool

    seq = seqs[i][20:290]

    name = pool_names[i]

    # Create fasta file with single sequence
    fa = open(path+f"/DMS_fa_files_all/{name}_FL.fa","w")
    fa.write(">"+name+"\n")
    fa.write(seq+"\n")
    fa.close()
    
