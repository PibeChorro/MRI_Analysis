#!/gpfs01/bartels/user/vplikat/anaconda3/bin/python

##########
# HEADER #
##########
# Before you run this script you have to have run the 
# 'searchlight_group_analysis.py' script.
# The underlying rationale and methods for the analysis in this script are 
# explained in detail in Stelzer, Chen & Turner (2012)
# The data for this script is derived from an fMRI experiment in which
# subjects viewed different videos. Either magic, control or surprise videos. 
# The experiment was devided into 3 blocks. Each block consisted of 4 
# experimental runs. In each run subjects viewed 24 videos (each video is 
# considered a trial).
# The videos in each block were associated with one object (Balls, Cards and 
# Sticks) and there were 3 magic effects (Appear, Change and Vanish). For each
# magic effect and object there are two trick versions (i.e. Appear1, Appear2,
# Change1,...). This resulted in 6 magic videos per object = 18 magic videos 
# and for every magic video there was a corresponding control video showing
# the same movements without the magical effect. Additionally per object there
# were 3 surprise videos showing unusual surprising actions performed with the
# objects (e.g. eating a playing card).
# After the second run in each block the underlying method behind each magic 
# trick was presented.

#                       TIME
#   ---------------------------------->
#   OBJECT1 R1  R2  Revelation  R3  R4  |
#   OBJECT2 R1  R2  Revelation  R3  R4  |   TIME
#   OBJECT3 R1  R2  Revelation  R3  R4  v

# RUNS: 2*Appear1 Magic 2*Appear2 Magic 2*Appear1 Control 2*Appear2 Control
#       2*Vanish1 Magic 2*Vanish2 Magic 2*Vanish1 Control 2*Vanish2 Control
#       2*Change1 Magic 2*Change2 Magic 2*Change1 Control 2*Change2 Control
#       2*Surprise1     2*Surprise2     2*Surprise13
#       = 24 Videos

# The aim of the experiment was to find neural correlates of surprise and in 
# particular of surprising events we consider "impossible". 
# The data used are beta estimate NIfTI images derived from a GLM using SPM12 
# in MATLAB and is in MNI space. 

##########################
# PURPOSE OF THIS SCRIPT #
##########################
# The purpose of this script is to calculate a t-test against chance level for
# the searchlight analysis performed on all subjects in MNI space

#############
# LIBRARIES #
#############

# interact with the operating system 
import os
import argparse
from pathlib import Path
import git
import glob
# data structuration and calculations
import numpy as np   # most important numerical calculations
# library for neuroimaging
import nibabel as nib
from nilearn.image import new_img_like,smooth_img
from nipype.interfaces.fsl.preprocess import SUSAN
# optimize time performance
import time
import scipy as sc

################################
# Command line input arguments #
################################

# create a parser object which handles the input variables
parser = argparse.ArgumentParser()

# add all the input arguments
parser.add_argument("--data",       "-d",   nargs="?",  const='pre',    
                    default='pre',  type=str)
parser.add_argument("--analyzed",   nargs='?', const='moment',  
                    default='moment',   type=str)

# parse the arguments to a parse-list(???)
ARGS = parser.parse_args()
# assign values 
DATA     = ARGS.data
ANALYZED = ARGS.analyzed

# variables for path selection and data access
HOME            = str(Path.home())
PROJ_DIR        = os.path.join(HOME, 'Documents/Master_Thesis/DATA/MRI')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')

RAWDATA_DIR     = os.path.join(PROJ_DIR, 'rawdata')
ANALYSIS        = 'ROI-analysis'
if DATA == 'pre':
    DATA_TO_USE = 'decode_effect_on_premagic'
    NUM_LABELS  = 3
elif DATA == 'post':
    DATA_TO_USE = 'decode_effect_on_postmagic'
    NUM_LABELS  = 3
elif DATA == 'all':
    DATA_TO_USE = 'decode_effect_on_allmagic'
    NUM_LABELS  = 3
elif DATA == 'pre-post':
    DATA_TO_USE = 'decode_pre_vs_post'
    NUM_LABELS  = 2
elif DATA == 'mag-nomag':
    DATA_TO_USE = 'magic_vs_nomagic'
    NUM_LABELS  = 2
else:
    raise
    
if ANALYZED == 'moment':
    data_analyzed = 'SpecialMoment'
elif ANALYZED == 'video':
    data_analyzed = 'WholeVideo'
else:
    raise
    
DATA_DIR        = os.path.join(DERIVATIVES_DIR, 'decoding', 'decoding_magic', 
                               DATA_TO_USE, data_analyzed,'SearchLight','LDA')
RESULTS_DIR     = os.path.join(DATA_DIR,'group-statistics')

SMOOTHING_SIZE = 4

if not os.path.isdir(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

SUBJECTS = glob.glob(os.path.join(DATA_DIR,'sub-*'))
SUBJECTS.sort()

smooth = SUSAN()
smooth.inputs.fwhm        = SMOOTHING_SIZE
smooth.inputs.brightness_threshold = 0.1
smooth.inputs.output_type = 'NIFTI'

for sub in SUBJECTS:
    smooth.inputs.in_file  = os.path.join(sub,'searchlight_results.nii') 
    smooth.inputs.out_file = os.path.join(sub,'s'+str(SMOOTHING_SIZE) +'searchlight_results.nii')
    smooth.run()
    
subs_accuracy_dirs = [os.path.join(sub,'s'+str(SMOOTHING_SIZE) +'searchlight_results.nii') 
                   for sub in SUBJECTS]
subs_accuracy_images = smooth_img(subs_accuracy_dirs, None)
subs_accuracy_map= []
for img in subs_accuracy_images:
    img_data    = img.get_fdata()                  # get data from NIfTI image
    img.uncache() 
    subs_accuracy_map.append(img_data)
    
# turn the accuracy map list into a 4D array (subjects * x * y * z)
# Get the dimensions, reshape it to get a 2D array (subjects * voxels)
# Calculate a one sided t-test against chance level for each voxel.
# Reshape the T and p-value maps back into a 3D matrix and save as nifti images
subs_accuracy_map = np.array(subs_accuracy_map)
image_shape = subs_accuracy_map.shape
subs_accuracy_data = subs_accuracy_map.reshape((image_shape[0],
                                                image_shape[1]*image_shape[2]*image_shape[3]))

# run t-tests on the data
Ts, ps = sc.stats.ttest_1samp(a=subs_accuracy_data, popmean=1/NUM_LABELS)
# correct p-values
ps[Ts>0] /= 2
ps[Ts<=0] = 1.0-(ps[Ts<=0]/2)

significant_voxels = (ps<=0.01) & (ps > 0)

T_map = Ts.reshape(image_shape[1],image_shape[2],image_shape[3])
p_map = ps.reshape(image_shape[1],image_shape[2],image_shape[3])
significant_map = significant_voxels.reshape(image_shape[1],image_shape[2],image_shape[3])
significant_map = significant_map.astype(float)

results = new_img_like(ref_niimg=subs_accuracy_images[0],data=T_map)
nib.save(results,os.path.join(DATA_DIR,'T_map.nii'))

results = new_img_like(ref_niimg=subs_accuracy_images[0],data=p_map)
nib.save(results,os.path.join(DATA_DIR,'uncorrected_p_map.nii'))

results = new_img_like(ref_niimg=subs_accuracy_images[0],data=significant_map)
nib.save(results,os.path.join(DATA_DIR,'uncorrected_significant_voxels.nii'))

##################
# WRITE LOG FILE #
##################
# We want to save all important information of the script execution
# To get the git hash we have to check if the script was run locally or on the
# cluster. If it is run on the cluster we want to get the $PBS_O_WORKDIR 
# variable, which preserves the location from which the job was started. 
# If it is run locally we want to get the current working directory.

try:
    script_file_directory = os.environ["PBS_O_WORKDIR"]
except KeyError:
    script_file_directory = os.getcwd()
    
try:
    rep = git.Repo(script_file_directory, search_parent_directories=True)
    git_hash = rep.head.object.hexsha
except git.InvalidGitRepositoryError:
    git_hash = 'not-found'
    
# create a log file, that saves some information about the run script
with open(os.path.join(RESULTS_DIR,'parametric_analysis-logfile.txt'), 'w+') as writer:
    writer.write('Codeversion: {} \n'.format(git_hash))
