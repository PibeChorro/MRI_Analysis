#!/usr/bin/bash

# define the directory where the raw data is stored
RAWDATA=$HOME/Documents/Master_Thesis/DATA/MRI/rawdata
# loop over every subject in RAWDATA and call qsub with docker_fs.sh to reconstruct the subject
for subject in ${RAWDATA}/sub-*
do
	x=$(basename "$subject")
	echo $x
	qsub -l nodes=1:ppn=1 -F "$x" docker_fs.sh -q shared
done
