#!/usr/bin/bash

# Executable bash file to queue on the headnode to run the freesurfer docker container for reconstruction
# VARIABLES
PROJ_DIR=$HOME/Documents/Master_Thesis/DATA/MRI
CONTAINER_NAME=agbartels_freesurfer_${1}

# RUN THE DOCKER IMAGE
docker run \
--detach \
--name $CONTAINER_NAME \
--rm \
--user "$(id -u):$(id -g)" \
--volume "${PROJ_DIR}:/home" \
--interactive \
--tty \
freesurfer/freesurfer:7.1.1
#docker run: executes docker container create and docker container start
#--detach: container runs in background. If this flag is not set you are stuck in the container terminal and all the following commands are not executed
#--name: A unique name for your container you can use to refer to it (start, stop, restart etc)
#--rm: removes container when stopped
#--user: IMPORTANT - docker runs in root by default. If this flag is not set, every file you create is owned by root and not by yourself
#--volume: kind of mounts a directory into the container
#--interactive: this flag is needed so you can send commands to your container
#--tty: opens a terminal in your container

# run docker exec to send the container the actual commands you want to execute
# some commands (like mkdir) can be given to 'docker exec', others however (like cd) need to be given to bash.
# Additionally we perform multiple commands like 'export'. A exec command does not know anything about a previously performed exec command.
docker exec -e SUB=$1 $CONTAINER_NAME bash -c 'export FS_LICENSE=/home/derivatives/freesurfer/.license; \
export SUBJECTS_DIR=/home/derivatives/freesurfer; \
cd $SUBJECTS_DIR; \
mkdir $SUBJECTS_DIR/${SUB}; \
mkdir $SUBJECTS_DIR/${SUB}/mri; \
mri_convert /home/rawdata/${SUB}/anat/${SUB}_T1w.nii /home/derivatives/freesurfer/${SUB}/mri/001.mgz; \
recon-all -autorecon-all -subjid ${SUB}'

# change the fs licence file path. The default fs path lies within root and there is no licence file given
# change the default fs SUBJECTS_DIR to where you want to store your data
# change into your SUBJECTS_DIR
# first create a folder for your subject
# then create a subfolder called 'mri'
# convert your raw anatomical nifti into an mgz file (which needs to be called 001.mgz)
# execute the recon-all command with the given subject

# stop your docker container, which is then automatically removed
docker stop $CONTAINER_NAME
