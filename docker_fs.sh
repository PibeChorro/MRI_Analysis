#!/usr/bin/bash

# Executable bash file to queue on the headnode to run the freesurfer docker container for reconstruction
# VARIABLES
PROJ_DIR=$HOME/Documents/Master_Thesis/DATA/MRI
CONTAINER_NAME=my_freesurfer1

# RUN THE DOCKER IMAGE
docker run \
-d \
--name $CONTAINER_NAME \
--rm \
--user "$(id -u):$(id -g)" \
-v "${PROJ_DIR}:/home" \
-it \
freesurfer/freesurfer:7.1.1

# RUN docker exec TO CREATE FOLDERS IN WHICH DATA CREATED BY
docker exec $CONTAINER_NAME bash -c 'export FS_LICENSE=/home/derivatives/freesurfer/.license; \
export SUBJECTS_DIR=/home/derivatives/freesurfer; \
cd $SUBJECTS_DIR; \
mkdir $SUBJECTS_DIR/sub-01; \
mkdir $SUBJECTS_DIR/sub-01/mri; \
mri_convert /home/rawdata/sub-01/anat/sub-01_T1w.nii /home/derivatives/freesurfer/sub-01/mri/001.mgz; \
recon-all -autorecon-all -subjid sub-01'

docker stop $CONTAINER_NAME
