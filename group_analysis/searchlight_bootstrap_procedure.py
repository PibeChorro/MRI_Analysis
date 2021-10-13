#!/gpfs01/bartels/user/vplikat/anaconda3/bin/python

##########
# HEADER #
##########
# Before you run this script you have to have run the 
# 'decode_magic_effects_seachlight.py' script.
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
# The purpose of this script is to generate a number of bootstrapped mean 
# chance accuracy maps, based on permutation maps from multiple subjects.

################################
# FUNCTIONALITY OF THIS SCRIPT #
################################
# FIRST STEP 
# Get all important path information and information about the analysis from 
# command line input.
# SECOND STEP
# From all subjects select one random permutation accuracy map and calculate 
# the mean accuracy for every voxel over the drawn permutation accuracy maps.
# Save the created mean accuracy map and repeat as often as it was stated in 
# the command line flag (10 000 in the original paper).

######################
# COMMAND LINE FLAGS #
######################
# --data: data from which previously calculated permutation tests should be 
# used. Accepted values are pre (magic effect decoding pre revelation), post
# (magic effect decoding post revelation), all (magic effect decoding using 
# pre and post revelation data), pre-post (decoding videos pre vs. post 
# revelation) and mag-nomag (decoding magic vs nomagic).
# --bootstraps: how many chance distribution maps should be created from the 
# 100 permutation maps per subject
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
from nilearn.image import new_img_like
from nipype.interfaces.fsl.preprocess import SUSAN
# optimize time performance
import time

# get start time
T_START = time.time()

################################
# Command line input arguments #
################################

# create a parser object which handles the input variables
parser = argparse.ArgumentParser()

# add all the input arguments
parser.add_argument("--data",       "-d",   nargs="?",  const='pre',    
                    default='pre',  type=str)
parser.add_argument("--analyzed",           nargs='?', const='moment',  
                    default='moment',   type=str)
parser.add_argument("--bootstraps", "-b",   nargs="?",  const=1000,     
                    default=10000,   type=int)   # how many bootstrapping draws are done for the null statistic

# parse the arguments to a parse-list(???)
ARGS = parser.parse_args()
# assign values 
DATA        = ARGS.data
BOOTSTRAPPS = ARGS.bootstraps
ANALYZED    = ARGS.analyzed
    
HOME            = str(Path.home())
PROJ_DIR        = os.path.join(HOME, 'Documents/Master_Thesis/DATA/MRI')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')

RAWDATA_DIR     = os.path.join(PROJ_DIR, 'rawdata')
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
    
SMOOTHING_SIZE = 4

GLM_DATA_DIR    = 'mnispace' 
FLA_DIR         = os.path.join(DERIVATIVES_DIR,'spm12',
                           'spm12-fla','WholeBrain',
                           'EveryVideo',GLM_DATA_DIR,
                           data_analyzed)
    
MASK_DIR    = os.path.join(FLA_DIR, 'group_mask.nii')
MASK_IMG    = nib.load(MASK_DIR)
MASK        = MASK_IMG.get_fdata()
MASK        = np.array(MASK,dtype=bool)
    
DATA_DIR        = os.path.join(DERIVATIVES_DIR, 'decoding', 'decoding_magic', 
                               DATA_TO_USE, data_analyzed,'SearchLight','LDA')
RESULTS_DIR     = os.path.join(DATA_DIR, 'group-statistics')
if not os.path.isdir(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)
    
# set random seed generator for reproducable results
rnd_seed = 1234
np.random.seed(rnd_seed)

null_distribution = []

SUBJECTS = glob.glob(os.path.join(DATA_DIR,'sub-*'))
SUBJECTS.sort()

smooth = SUSAN()
smooth.inputs.fwhm        = SMOOTHING_SIZE
smooth.inputs.brightness_threshold = 0.1
smooth.inputs.output_type = 'NIFTI'

# get the real mean decoding accuracy 
mean_accuracy_map = []
for sub in SUBJECTS:
    img         = nib.load(os.path.join(sub,'searchlight_results.nii'))
    img_data    = img.get_fdata()                  # get data from NIfTI image
    img_data[~MASK] = 0
    img.uncache() 
    mean_accuracy_map.append(img_data)
    
# convert list to array and calculate the mean
mean_accuracy_map = np.array(mean_accuracy_map)
mean_accuracy_map = mean_accuracy_map.mean(axis=0)
mean_accuracy_map[~MASK] = 0
# save as NIfTI file
results = new_img_like(ref_niimg=img, data=mean_accuracy_map)
nib.save(results,os.path.join(RESULTS_DIR, 'mean_accuracy.nii'))
smooth.inputs.in_file     = os.path.join(RESULTS_DIR,'mean_accuracy.nii')
smooth.inputs.out_file = os.path.join(RESULTS_DIR,'s'+str(SMOOTHING_SIZE) +'mean_accuracy.nii')
smooth.run()
os.remove(os.path.join(RESULTS_DIR, 'mean_accuracy.nii'))
# Using the inference method described in Stelzer (2013)
# From each subject randomly select one permutation image and create a mean
# over subjects. Do this several times (10^5 was used in Stelzer) and get a 
# distribution
for draw in range(BOOTSTRAPPS):
    # empty list to store permuted maps into
    perm_classification_images = []
    for s, sub in enumerate(SUBJECTS):
        # get all permuted accuracy maps of current subject and sort them
        perm_images     = glob.glob(os.path.join(sub,'perm*.nii'))
        perm_images.sort()
        # get a random number between 0 and number of permuted maps
        rnd_perm_map    = np.random.randint(0,len(perm_images))
        # read in the randomly selected accuracy map, transform into numpy array and add it to the list
        img             = nib.load(perm_images[rnd_perm_map])
        img_data        = img.get_fdata()                  # get data from NIfTI image
        img.uncache() 
        perm_classification_images.append(img_data)
        
    # the same as with the right accuracy maps: turn list into numpy array, calculate mean and save as NIfTI file
    perm_classification_images  = np.array(perm_classification_images)
    perm_mean_accuracy_map      = perm_classification_images.mean(axis=0)
    perm_mean_accuracy_map[~MASK] = 0
    results = new_img_like(ref_niimg=img,data=perm_mean_accuracy_map)
    
    nib.save(results,
             os.path.join(RESULTS_DIR,
                          'bootstrapped_{:04d}mean_accuracy.nii'.format(draw)))
    smooth.inputs.in_file     = os.path.join(RESULTS_DIR,
                          'bootstrapped_{:04d}mean_accuracy.nii'.format(draw))
    smooth.inputs.out_file = os.path.join(RESULTS_DIR,
                          's'+str(SMOOTHING_SIZE) +'bootstrapped_{:04d}mean_accuracy.nii'.format(draw))
    smooth.run()
    os.remove(os.path.join(RESULTS_DIR,
                          'bootstrapped_{:04d}mean_accuracy.nii'.format(draw)))

rep = git.Repo(os.getcwd(),search_parent_directories=True)
git_hash = rep.head.object.hexsha

# create a log file, that saves some information about the run script
with open(os.path.join(RESULTS_DIR,'bootstrapping-logfile.txt'), 'w+') as writer:
    writer.write('Codeversion: {} \n'.format(git_hash))
    writer.write('Random seed generator: {}\n'.format(rnd_seed))
    writer.write('Time for computation: {}h'.format(str((time.time() - T_START)/3600)))