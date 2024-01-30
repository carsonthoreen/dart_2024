#!/bin/bash
#Pipeline to analyze DART sequencing results.

#v2: 1. updated directories for Mccleary. 2. Add Reference for New pools. 3. Restrict max mismatch to 1. 

module load BBMap
module load STAR

LibType=$1
if [[ $LibType == "Long" ]] 
	then RefDirectory="/gpfs/gibbs/project/thoreen/lhx3/Reference/STAR_reference/LongandVar/STARgenome"
	elif [[ $LibType == "Short" ]]
		then RefDirectory="/gpfs/gibbs/project/thoreen/lhx3/Reference/STAR_reference/Short/STARgenome"
	elif [[ $LibType == "New" ]]
		then RefDirectory="/gpfs/gibbs/project/thoreen/lhx3/Reference/STAR_reference/New/STARindex"
	else echo "Please specify library type: Long, Short or New"
		exit 1
fi



file_R1=$2
file_R2=$3

echo $file_R1 $file_R2

echo "Adaptor trimming"
bbduk.sh in1=${file_R1} in2=${file_R2} out1=${file_R1}.trimmed out2=${file_R2}.trimmed \
ref=/gpfs/gibbs/project/thoreen/lhx3/Documents/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo

echo "Collapse duplicates"
clumpify.sh in1=${file_R1}.trimmed in2=${file_R2}.trimmed out1=${file_R1}.collapsed out2=${file_R2}.collapsed dedupe subs=0

echo "UMI removal"
bbduk.sh in1=${file_R1}.collapsed in2=${file_R2}.collapsed out1=${file_R1}.clean out2=${file_R2}.clean ftl=10

echo "STAR mapping"
STAR --runThreadN 4 --genomeDir $RefDirectory --readFilesIn ${file_R1}.clean ${file_R2}.clean --outFileNamePrefix CountFile. \
--soloStrand Forward --alignSJoverhangMin 999 --alignIntronMax 1 --alignIntronMin 999 --outFilterMismatchNmax 1  

echo "Counting"
pileup.sh in=CountFile.Aligned.out.sam out=Count_STAR_coverage.txt

rm -rf *.trimmed *.collapsed *.clean
