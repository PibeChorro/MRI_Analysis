#!/gpfs01/bartels/user/vplikat/anaconda3/bin/python

##########
# HEADER #
##########
# Searchlight decoding analysis with permutations 
# The data used are beta estimate NIfTI images derived from a GLM using SPM12 
# in MATLAB. For group analysis purpose we use data normalized to MNI space.
# The experiment was devided into 3 blocks. Each block consisted of 4 
# experimental runs. In each run subjects viewed 24 videos (each video is 
# considered a trial) of three different categories: Magic, Control and 
# Surprise.
# The videos in each block were associated with one object (Balls, Cards and 
# Sticks) and there were 3 magic effects (Appear, Change and Vanish). For each
# magic effect and object there are two trick versions (i.e. Appear1, Appear2,
# Change1,...).
# There were 6 Magic videos (presented twice each), the other 12 videos are of
# no interest here. 
# After the second run in each block the underlying method behind each magic 
# trick was presented.
# Here I run a searchlight decoding analysis on subject level using either
# Support Vector Machine (SVM) or Linear Discriminant Analysis (LDA) 
# classification. The aim is to predict magic effect trained on two objects
# and test it on the remaining third.

#############
# LIBRARIES #
#############

# interact with the operating system 
import os
import argparse
from pathlib import Path
import git
# data structuration and calculations
import pandas as pd  # to create data frames
import numpy as np   # most important numerical calculations
# read in mat files
import readmat
# needed to extract the run number out of the parentesis of the string in the SPM.mat file
import re
# library for neuroimaging
import nibabel as nib
from nilearn.decoding import SearchLight
from nilearn.image import smooth_img, new_img_like
# machine learning algorithms and stuff
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.svm import SVC
from sklearn.model_selection import PredefinedSplit
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
parser.add_argument("--sub",        "-s",                               default='sub-01')           # subject
parser.add_argument("--smooth",             nargs='?',  const=0,        default=0,      type=int)   # what data should be used
parser.add_argument("--algorythm",  "-a",   nargs='?',  const='LDA',    default='LDA',  type=str)
parser.add_argument("--kernels",    "-k",   nargs='?',  const=12,       default=12,     type=int)   # how many processes should be run in parallel
parser.add_argument("--runs",       "-r",   nargs="?",  const='pre',    default='pre',  type=str)
parser.add_argument("--perms",      "-p",   nargs="?",  const=10,       default=10,     type=int)   # how many permutations
# parse the arguments to a parse-list(???)
ARGS = parser.parse_args()
# assign values 
SUB             = ARGS.sub
N_PROC          = ARGS.kernels
SMOOTHING_SIZE  = ARGS.smooth
DECODER         = ARGS.algorythm
RUNS_TO_USE     = ARGS.runs
N_PERMS         = ARGS.perms

# initiate the decoder based on command line input (either LDA (default) or SVM)
if DECODER =='LDA':
    my_decoder          = LDA(solver='lsqr', shrinkage='auto')
elif DECODER == 'SVM':
    SVM_C = 1
    my_decoder          = SVC(kernel='linear', C=SVM_C)
else:
    raise
        
# decide what data (from what runs) shall be used based on command line input
# 'pre' takes data before revealing the method behind the magic tricks, 
# 'post' takes data after revealing the method behind the magic tricks and
# 'all' takes all data pre and post together
if RUNS_TO_USE == 'pre':
    runs_of_interest = [1,2,5,6,9,10]
elif RUNS_TO_USE == 'post':
    runs_of_interest = [3,4,7,8,11,12] 
elif RUNS_TO_USE == 'all':
    runs_of_interest = [1,2,3,4,5,6,7,8,9,10,11,12] 
else:
    raise
        

# variables for path selection and data access
HOME            = str(Path.home())
PROJ_DIR        = os.path.join(HOME, 'Documents/Master_Thesis/DATA/MRI')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')
if SMOOTHING_SIZE > 0:
    GLM_DATA_DIR    = str(SMOOTHING_SIZE)+'mm-smoothed-mnispace' 
    FLA_DIR         = os.path.join(DERIVATIVES_DIR,'spm12',
                               'spm12-fla','WholeBrain',
                               'EveryVideo',GLM_DATA_DIR,
                               'WholeVideo',SUB)
else:
    GLM_DATA_DIR    = 'mnispace' 
    FLA_DIR         = os.path.join(DERIVATIVES_DIR,'spm12',
                               'spm12-fla','WholeBrain',
                               'EveryVideo',GLM_DATA_DIR,
                               'WholeVideo',SUB)
FREESURFER_DIR  = os.path.join(DERIVATIVES_DIR, 'freesurfer')
RAWDATA_DIR     = os.path.join(PROJ_DIR, 'rawdata')
MASK_DIR        = os.path.join(FLA_DIR, 'mask.nii')
ROI_DIR         = os.path.join(FREESURFER_DIR,SUB,'corrected_ROIs')
SPM_MAT_DIR     = os.path.join(FLA_DIR, 'SPM.mat')
ANALYSIS        = 'SearchLight'
RESULTS_DIR     = os.path.join(DERIVATIVES_DIR, 'decoding', 'decoding_magic', 
                               'decode_effect_on_'+RUNS_TO_USE+'magic', 
                               'SpecialMoment', ANALYSIS, DECODER, SUB)
OUTPUT_DIR = os.path.join(RESULTS_DIR, 'searchlight_results.nii')   # where to store the results
if not os.path.isdir(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

# The labels to be predicted from the decoder
LABEL_NAMES = [
    'Appear',
    'Change',
    'Vanish'
]
    
print ('Analysing subject: {}'.format(SUB))
print ('Getting ROIs from:	 {}'.format(FREESURFER_DIR))
print ('Saving data at:	 {}'.format(RESULTS_DIR))

########################################
# reading in the necessary information #
########################################

# From the previously created SPM.mat file we read in the information we need
# The filenames of our beta images
SPM_BETADICT = readmat.load(SPM_MAT_DIR, isStruct=True)['SPM']['Vbeta']
BETA_DIRS = [f['fname'] for f in SPM_BETADICT]
# The corresponding Regressor names - are in form of 'Sn(<run-number>) <Regressor-Name>*bf(1)'
SPM_REGRESSORS = readmat.load(SPM_MAT_DIR,isStruct=True)['SPM']['xX']['name']
# again a complex process to throw out regressors of no interest (like realignment)
REGRESSORS_OF_INTEREST  = [True if 'Magic' in n else False for n in SPM_REGRESSORS]

# store beta filenames and regressornames in a dictionary
DATA_DICT = {
    'Regressors': SPM_REGRESSORS,
    'BetaNames': BETA_DIRS
}

# convert dictionary into a pandas DataFrame for further analysis
label_df = pd.DataFrame(DATA_DICT, columns=DATA_DICT.keys())

# This complex list comprehensions are necessary to get the run number out of the regressor name
x       = [' '.join(re.findall(r"\((\d+)\)",string)) for string in label_df.Regressors]
runs    = [int(s_filter.split()[0]) for s_filter in x]

# add further data to DataFrame
label_df['Runs']    = runs                  # In which run
label_df['Chunks']  = (label_df.Runs-1)//4  # The chunks (needed for cross validation)
label_df['Labels']  = np.nan                # Labels will be filled later 

# throw out all rows of regressors of no interest
label_df = label_df.iloc[REGRESSORS_OF_INTEREST]
label_df = label_df[label_df.Runs.isin(runs_of_interest)]

# Check for every entry in Regressors if it contains one of the label names. If so, assign the label name
for l in LABEL_NAMES:
    label_df.Labels = np.where(label_df.Regressors.str.contains(l),l,label_df.Labels)

#######################
# THE ACTUAL DECODING #
#######################
targets                 = np.asarray(label_df.Labels)   # get labels as numpy array from pandas dataframe
chunks                  = np.asarray(label_df.Chunks)    # get chunks for cross validation as numpy array from data frame
runs_for_permutation    = np.asarray(label_df.Runs)

# initialize the searchlight decoding object
MY_SEARCH_LIGHT = SearchLight(mask_img=MASK_DIR,
                             radius=4.0,
                             estimator=my_decoder,
                             n_jobs=N_PROC,
                             cv=PredefinedSplit(chunks),
                             verbose=3)
betas = smooth_img(FLA_DIR+'/'+label_df.BetaNames,fwhm=None)
# fit the decoding object based on the previously loaded betas and labes
# perform cross validation over objects
MY_SEARCH_LIGHT.fit(imgs=betas,y=targets,groups=chunks)

# Form results into a NIfTI 
results = new_img_like(ref_niimg=betas[0],data=MY_SEARCH_LIGHT.scores_)
nib.save(results,OUTPUT_DIR)

################
# PERMUTATIONS #
################
# shuffle the labels within each run and start a new decoding based on the same
# data. Store data again as a NIfTI image
for i in range(N_PERMS):
    perm_results_dir    = os.path.join(RESULTS_DIR, 
                                       'perm_{:04d}_searchlight_results.nii'.format(i))   # where to store the results
    permed_targets      = []    # empty list which will be filled with permuted labels
    for r in runs_of_interest:
        tmp = targets[runs_for_permutation==r]  # get labels of one run
        permed_targets.extend(np.random.permutation(tmp))   # permute labels and add to list
    
    # searchlight decoding with permuted labels
    MY_SEARCH_LIGHT.fit(imgs=betas,y=permed_targets,groups=chunks)
    results = new_img_like(ref_niimg=betas[0],data=MY_SEARCH_LIGHT.scores_)
    nib.save(results,perm_results_dir)
    
del label_df
del betas


rep = git.Repo(os.getcwd(),search_parent_directories=True)
git_hash = rep.head.object.hexsha

# create a log file, that saves some information about the run script
with open(os.path.join(RESULTS_DIR,'serchlight-logfile.txt'), 'w+') as writer:
    writer.write('Codeversion: {} \n'.format(git_hash))
    writer.write('Number of permutations: {}\n'.format(str(N_PERMS)))
    writer.write('Number of kernels used: {}\n'.format(str(N_PROC)))
    writer.write('Time for computation: {}h'.format(str((time.time() - T_START)/3600)))

# print time the whole processe took
print ((time.time() - T_START)/3600)