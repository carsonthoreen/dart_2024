#!/usr/bin/env python3

'''
Go through list of sequence IDs, extract sequences and pairing probabilities from corresponding dp.ps files (RNAfold output).

Input: a txt file with a list of sequence IDs, generated from sed command on RNAfold output files.
Output:
One .csv file:
    column 1 = seq_name
    column 2 = pos
    column 3 = ppair

Note: put this script in the same directory with all dp.ps files and the txt file before executing.
'''

#import numpy as np
import csv

def parse_structure(name, out_file):

    print(name)
    ps_file = name+'_FL_0001_dp.ps'

    with open(ps_file) as f:
        # Iterate over lines until you find the sequence, then store it
        while True:
            line = f.readline()
            if line.startswith("/sequence"):
                break
        # sequence spans two lines; combine the lines
        seq_line1 = f.readline().strip()[:-1]
        seq_line2 = f.readline().strip()[:-1]
        sequence = seq_line1 + seq_line2
        print(sequence)
        print('sequence length: '+str(len(sequence)))

        # Continue iterating until you find the beginning of the base
        # pair data, then initialize a matrix to store it
        while True:
            line = f.readline()
            if line.startswith("%start of base pair"):
                break
        structure = [ [0.0 for i in range(len(sequence))]
                     for j in range(len(sequence)) ]
        # Iterate over basepair data, only storing the base pairs
        # from the ubox
        while True:
            line = f.readline()
            try:
                i,j,prob,box = line.strip().split()
            except ValueError:
                break
            if box == "lbox":
                continue
            i,j = int(i)-1, int(j)-1
            # reported value is sqrt(prob), square it to get prob
            prob = float(prob)**2
            structure[i][j] = prob
            structure[j][i] = prob

    ppair = [ sum(row) for row in structure]

    # make seq_name, pos columns
    seq_name = [name] * len(sequence)
    #pos = np.arange(1, len(sequence)+1)
    pos = [i+1 for i in range(len(sequence))]

    # export all columns to new file
    with open(out_file, 'a') as c:
        writer = csv.writer(c, delimiter='\t')
        writer.writerows(zip(seq_name, pos, ppair))
        c.close()

# main function -------------------------
for name in open('names_noFluc.txt', 'r').readlines():
    name = name.strip('\n')
    parse_structure(name, out_file='ppairs.csv')
