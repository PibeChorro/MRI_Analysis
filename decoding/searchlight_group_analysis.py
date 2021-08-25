#!/gpfs01/bartels/user/vplikat/anaconda3/bin/python

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
parser.add_argument("--runs",   "-r",   nargs="?",  const='pre',    default='pre',  type=str)
parser.add_argument("--draws",  "-d",   nargs="?",  const=100,    default=100,  type=int)   # how many bootstrapping draws are done for the null statistic

# parse the arguments to a parse-list(???)
ARGS = parser.parse_args()
# assign values 
RUNS_TO_USE     = ARGS.runs
BOOTSTRAPPS     = ARGS.draws

HOME            = str(Path.home())
PROJ_DIR        = os.path.join(HOME, 'Documents/Master_Thesis/DATA/MRI')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')

RAWDATA_DIR     = os.path.join(PROJ_DIR, 'rawdata')
ANALYSIS        = 'ROI-analysis'
DATA_DIR        = os.path.join(DERIVATIVES_DIR, 'decoding', 'decoding_magic', 
                               'decode_effect_on_'+RUNS_TO_USE+'magic',
                               'SpecialMoment','SearchLight_old','LDA')
RESULTS_DIR     = os.path.join(DATA_DIR,'group-statistics')
if not os.path.isdir(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

null_distribution = []

SUBJECTS = glob.glob(os.path.join(DATA_DIR,'sub-*'))
SUBJECTS.sort()

# get the real mean decoding accuracy 
mean_accuracy_map = []
for sub in SUBJECTS:
    img         = nib.load(os.path.join(sub,'searchlight_results.nii'))
    img_data    = img.get_fdata()                  # get data from NIfTI image
    img.uncache() 
    mean_accuracy_map.append(img_data)
    
# convert list to array and calculate the mean
mean_accuracy_map = np.array(mean_accuracy_map)
mean_accuracy_map = mean_accuracy_map.mean(axis=0)
# save as NIfTI file
results = new_img_like(ref_niimg=img,data=mean_accuracy_map)
nib.save(results,os.path.join(RESULTS_DIR,'mean_accuracy.nii'))

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
    results = new_img_like(ref_niimg=img,data=perm_mean_accuracy_map)
    
    nib.save(results,os.path.join(RESULTS_DIR,
                                  'perm_{:04d}mean_accuracy.nii'.format(draw)))
    
print ('Done')