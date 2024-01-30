#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt

def plot_histogram(tsv_file1,tsv_file2,png_file):

    #Read the TSV file into a dataframe
    df1 = pd.read_csv(tsv_file1, header=None, sep="\t")
    df2 = pd.read_csv(tsv_file2, header=None, sep="\t")

    #Extract the 3rd column that contains pairing probability data
    colname1 = df1.columns[2]
    mutrate1 = df1[colname1]
    colname2 = df2.columns[2]
    mutrate2 = df2[colname2]

    #Remove all nan values
    mutrate_nonan1 = mutrate1.dropna()
    mutrate_nonan2 = mutrate2.dropna()

    print(mutrate_nonan1)
    print(mutrate_nonan2)

    #Plot combined histogram
    fig,ax = plt.subplots(figsize=(8,6))
    ax.hist([mutrate_nonan1,mutrate_nonan2],bins = 1000,label = ["AC","GT"])
    ax.set_xlabel("Per-base mutation rate")
    ax.set_ylabel("Number of bases")
    ax.set_xlim([0,0.15])
    ax.set_ylim([0,150000])
    ax.set_title("Denatured") # CHANGE TITLE AS NEEDED!
    ax.legend()

    #Write histogram to a png graphics file. CHANGE INPUT OUTPUT FILE NAMES AS NEEDED!
    fig.savefig(png_file)
 
plot_histogram("raw_reacts/denatured_rawreactsAC.tsv","raw_reacts/denatured_rawreactsGT.tsv","raw_reacts/rawreacts_denatured.png")
