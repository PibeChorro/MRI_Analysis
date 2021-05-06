#!/usr/bin/bash

# define the directory where the raw data is stored
FSDATA=$HOME/Documents/Master_Thesis/DATA/MRI/derivatives/freesurfer
# loop over every subject in RAWDATA and call qsub with docker_fs.sh to reconstruct the subject
for subject in ${FSDATA}/sub-*
do
	x=$(basename "$subject")
	echo $x
	qsub -l nodes=1:ppn=1 -F "$x" assign_ROI_overlaps.py -q shared
done
