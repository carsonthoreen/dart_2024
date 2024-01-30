#!/bin/bash

#SBATCH --partition=ycga
#SBATCH --job-name=VARNA
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=10GB
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=4

module purge

file="sequences.txt"

while read -r line; do
	name=${line}
	echo ${name}

	sequence_FL=$(sed '2q;d' "${name}"_FL.fold)
	echo ${sequence_FL}

	struct_FL=$(sed '3q;d' "${name}"_FL.fold)
	struct_FL=${struct_FL:0:-9}
	echo ${struct_FL}

	dms_FL=${name}_DMS_FL.csv
	reacts_FL=$(cut -f2 "${dms_FL}" | tr -s '\r\n' ';')
	echo ${reacts_FL:0:-1}

	outfile=${name}_FL.svg


	VARNA='java -cp /home/dj448/pi_wgilbert/users_data/dannij/software/VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd'

	$VARNA -sequenceDBN ${sequence_FL} -structureDBN ${struct_FL} -o ${outfile} -algorithm naview -colorMapMax '1.0' -colorMapMin '0' -colorMap ${reacts_FL:0:-1} -colorMapStyle '0:#FFFFFF;1:#ff0000' -bp '#000000' -baseOutline '#000000' -spaceBetweenBases '0.8'

done <$file
