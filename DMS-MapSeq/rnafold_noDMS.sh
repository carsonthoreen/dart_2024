#!/bin/bash

#SBATCH --partition=ycga
#SBATCH --job-name=calcNFold
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=10GB
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=4

#module load Biopython

#echo CALC REACTIVITIES, SUBSET FASTA
#python3 calcDMSConstraints_allseqs.py

module load ViennaRNA

for i in *FL.fa; do

	fasta_FL=${i}
	name=${i:0:-5}

	echo ${fasta_FL}
	echo ${name}

	RNAfold -p -T 37 --infile=${fasta_FL} --outfile=${name}FL.fold --id-prefix=${name}FL

done
