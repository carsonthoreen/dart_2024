#!/bin/bash
#Demultiplexing.sh

module load miniconda
conda activate Flexbar_env

DataLeft=$1
DataRight=$2
echo "Left: $DataLeft"
echo "Right: $DataRight"

flexbar -r $DataRight -p $DataLeft -b /home/lhx3/project/Documents/Barcode.fa -bt LTAIL -be 0.2  --target ${DataLeft%_*} -O ${DataLeft%_*}.Log.out
