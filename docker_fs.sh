#!/usr/bin/bash

# Executable bash file to queue on the headnode to run the freesurfer docker container for reconstruction
# VARIABLES
SUB=${1}
PROJ_DIR=$HOME/Documents/Master_Thesis/DATA/MRI
FS_CONTAINER_NAME=agbartels_freesurfer_${SUB}
SUBJECTS_DIR=/home/derivatives/freesurfer
COREGISTERED_DIR=/home/derivatives/spm12/spm12-preproc/coregistered
# Those are the ROI names you just need to know them -- V1v=1, V1d=2 in the label images (... I guess ???)
# (https://hub.docker.com/r/nben/occipital_atlas) "older version" of neuropythy)
roiname_array=("V1v" "V1d" "V2v" "V2d" "V3v" "V3d" "hV4" "VO1" "VO2" "PHC1" "PHC2" \
"TO2" "TO1" "LO2" "LO1" "V3B" "V3A" "IPS0" "IPS1" "IPS2" "IPS3" "IPS4" \
"IPS5" "SPL1" "FEF")
# ROIs we want to merge
rois_to_merge=("V1" "V2" "V3")
# all ROIs combined -- this strange syntax is needed
allROIs=( "${roiname_array[@]}" "${rois_to_merge[@]}" )

# RUN THE DOCKER IMAGE
docker run \
--detach \
--name $FS_CONTAINER_NAME \
--rm \
--user "$(id -u):$(id -g)" \
--volume "${PROJ_DIR}:/home" \
--env FS_LICENSE=/home/derivatives/freesurfer/.license \
--env SUBJECTS_DIR=$SUBJECTS_DIR \
--interactive \
--tty \
freesurfer/freesurfer:7.1.1
#docker run: executes docker container create and docker container start
#--detach: container runs in background. If this flag is not set you are stuck in the container terminal and all the following commands are not executed
#--name: A unique name for your container you can use to refer to it (start, stop, restart etc)
#--rm: removes container when stopped
#--user: IMPORTANT - docker runs in root by default. If this flag is not set, every file you create is owned by root and not by yourself
#--volume: kind of mounts a directory into the container
#--env: exports/overwrites an environment variable - here where FS should look for the license file and subjects
#--interactive: this flag is needed so you can send commands to your container
#--tty: opens a terminal in your container

# run docker exec to send the container the actual commands you want to execute
# some commands (like mkdir) can be given to 'docker exec', others however (like cd) need to be given to bash (docker exec bash -c 'bla bla bla').

docker exec $FS_CONTAINER_NAME mkdir $SUBJECTS_DIR/${SUB}
# create a folder for your subject
docker exec $FS_CONTAINER_NAME mkdir $SUBJECTS_DIR/${SUB}/mri
# create a subfolder called 'mri'
docker exec $FS_CONTAINER_NAME mri_convert /home/rawdata/${SUB}/anat/${SUB}_T1w.nii /home/derivatives/freesurfer/${SUB}/mri/001.mgz
# convert your raw anatomical nifti into an mgz file (which needs to be called 001.mgz)
docker exec --workdir $SUBJECTS_DIR $FS_CONTAINER_NAME recon-all -autorecon-all -subjid ${SUB}
# execute the recon-all command with the given subject

# perform parcelation using noah bensons neuropythy 
python -m neuropythy atlas --verbose $PROJ_DIR/derivatives/freesurfer2/${SUB}

docker stop $FS_CONTAINER_NAME
