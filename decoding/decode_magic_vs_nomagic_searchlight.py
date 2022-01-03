#!/gpfs01/bartels/user/vplikat/anaconda3/bin/python

##########
# HEADER #
##########
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
# The purpose of this script is to use a machine learning algorithm to predict
# the magical effects performed on one object, based on training data from the
# other two objects. To be able to test for statistical significance, 
# permutation testing is applied.
# All this is done for the whole brain, but on a single subject level.
# In order to analize the whole brain a small sphere is used as a mask. This
# sphere moves through the whole brain and the decoding accuracy value obtained
# is assigned to the voxel in the center of the sphere.

################################
# FUNCTIONALITY OF THIS SCRIPT #
################################
# FIRST STEP 
# Get all important path information and names of magic EFFECTS (Appear, 
# Change, Vanish) and information about the analysis from command line input.
# SECOND STEP 
# Read in the SPM.mat file created by SPM12 when the GLM is estimated.
# From this SPM.mat file we read out the names of the beta NIfTI images and
# the names of the regressors, that correspond to the beta image.
# From the Regressor names the run number (1-12) is extracted and added to the
# DataFrame. Then all regressors of NO interest (realignment, controll and 
# surprise videos, etc.) are removed, as well as the data we do not want to
# test (when analyzing pre revelation data, the post revelation data is removed
# and vise versa).
# Again from the Regressor name the label (=the effect) is extracted and added
# to the dataframe. Finally the chunks are defined. Chunks are needed for cross
# validation (meaning the whole dataset is sperated in chunks and each chunk is
# used for testing once).
# THIRD STEP
# Perform one Searchlight analysis using the 'real' labels and save the result
# in a NIfTI file. 
# Then loop over the number of permutations provided by the command line input
# and every time shuffle the labels WITHIN THE RUNS (!!!) and again save the
# result as perm_XXXX_searchlight_result in a NIfTI file.

######################
# COMMAND LINE FLAGS #
######################
# --sub: the subject ID that shall be analyzed
# --smooth: if the script should read in beta images that are the result of a
# GLM based on smoothed functional images
# --algorythm: which algorythm should be used. Currently implemented SVM and LDA
# --kernels: How many kernels should be used to parallize the permutation testing
# --runs: Which data should be used. Either pre, post revelation or all data
# together
# --perms: How many permutations are applied

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

# FIRST STEP
################################
# Command line input arguments #
################################

# create a parser object which handles the input variables
parser = argparse.ArgumentParser()

# add all the input arguments
parser.add_argument("--sub",        "-s",                               
                    default='sub-01')           # subject
parser.add_argument("--over",       "-o",   nargs='?',  const='objects',    
                    default='objects')
parser.add_argument("--smooth",             nargs='?',  const=0,        
                    default=0,      type=int)   # what data should be used
parser.add_argument("--algorythm",  "-a",   nargs='?',  const='LDA',    
                    default='LDA',  type=str)
parser.add_argument("--kernels",    "-k",   nargs='?',  const=12,       
                    default=12,     type=int)   # how many processes should be run in parallel
parser.add_argument("--runs",       "-r",   nargs="?",  const='pre',    
                    default='pre',  type=str)
parser.add_argument("--analyzed",           nargs='?', const='moment',  
                    default='moment',   type=str)
parser.add_argument("--perms",      "-p",   nargs="?",  const=10,       
                    default=10,     type=int)   # how many permutations
parser.add_argument("--radius",             nargs="?", const=4.0,
                    default=4.0,    type=float)
# parse the arguments to a parse-list(???)
ARGS = parser.parse_args()
# assign values 
SUB             = ARGS.sub
OVER            = ARGS.over
SMOOTHING_SIZE  = ARGS.smooth
DECODER         = ARGS.algorythm
N_PROC          = ARGS.kernels
RUNS_TO_USE     = ARGS.runs
ANALYZED        = ARGS.analyzed
N_PERMS         = ARGS.perms
SEARCHLIG_RAD   = ARGS.radius

# based on command line flag decide what decoding algorithm should be used
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
        
if ANALYZED == 'moment':
    data_analyzed = 'SpecialMoment'
elif ANALYZED == 'video':
    data_analyzed = 'WholeVideo'
else:
    raise
    
################################################
# VARIABLES FOR PATH SELECTION AND DATA ACCESS #
################################################
HOME            = str(Path.home())
PROJ_DIR        = os.path.join(HOME, 'Documents/Magic_fMRI/DATA/MRI')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')
# where to look for the beta images
if SMOOTHING_SIZE > 0:
    GLM_DATA_DIR    = str(SMOOTHING_SIZE)+'mm-smoothed-mnispace' 
    FLA_DIR         = os.path.join(DERIVATIVES_DIR,'spm12',
                               'spm12-fla','WholeBrain',
                               'EveryVideo',GLM_DATA_DIR,
                               data_analyzed)
else:
    GLM_DATA_DIR    = 'mnispace' 
    FLA_DIR         = os.path.join(DERIVATIVES_DIR,'spm12',
                               'spm12-fla','WholeBrain',
                               'EveryVideo',GLM_DATA_DIR,
                               data_analyzed)
# where the brain mask for the subjects can be found
MASK_DIR        = os.path.join(FLA_DIR, 'group_mask.nii')
# where the .mat file can be found created by SPM12
SPM_MAT_DIR     = os.path.join(FLA_DIR, SUB, 'SPM.mat')
# wheret to store the results
ANALYSIS        = 'SearchLight'
RESULTS_DIR     = os.path.join(DERIVATIVES_DIR, 'decoding', 'decoding_magic', 
                               'decode_magic_vs_nomagic',''+RUNS_TO_USE+'_videos',  
                               'over_' + OVER, data_analyzed , ANALYSIS, 
                               DECODER, SUB)
OUTPUT_DIR = os.path.join(RESULTS_DIR, 'searchlight_results.nii')   # where to store the results
if not os.path.isdir(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)


EFFECT_NAMES = [
    'Appear',
    'Change',
    'Vanish'
]
    
print ('Analysing subject: {}'.format(SUB))
print ('Saving data at:	 {}'.format(RESULTS_DIR))

# SECOND STEP 
########################################
# reading in the necessary information #
########################################

# From the previously created SPM.mat file we read in the information we need
# The filenames of our beta images
SPM_BETADICT = readmat.load(SPM_MAT_DIR, isStruct=True)['SPM']['Vbeta']
BETA_DIRS = [f['fname'] for f in SPM_BETADICT]
# The corresponding Regressor names - are in form of 'Sn(<run-number>) <Regressor-Name>*bf(1)'
SPM_REGRESSORS = readmat.load(SPM_MAT_DIR,isStruct=True)['SPM']['xX']['name']

# store beta filenames and regressornames in a dictionary
DATA_DICT = {
    'Regressors': SPM_REGRESSORS,
    'BetaNames': BETA_DIRS
}

# convert dictionary into a pandas DataFrame for further analysis
label_df = pd.DataFrame(DATA_DICT, columns=DATA_DICT.keys())

# This complex loop is necessary to get the run number out of the regressor name
x       = [' '.join(re.findall(r"\((\d+)\)",string)) for string in label_df.Regressors]
runs    = [int(s_filter.split()[0]) for s_filter in x]

# add further data to DataFrame
label_df['Runs']    = runs                  # In which run
label_df['Chunks']  = (label_df.Runs-1)//4  # The chunks (needed for cross validation)
label_df['Labels']  = np.nan                # Labels

# again a complex process to throw out regressors of no interest (like realignment)
regressors_of_interest  = [True if ('Magic' in n) 
                           or ('Control' in n) 
                           or ('Surprise' in n) 
                           else False for n in SPM_REGRESSORS]
# throw out all rows of regressors of no interest
label_df = label_df.iloc[regressors_of_interest]
label_df = label_df[label_df.Runs.isin(runs_of_interest)]
# Check for every entry in Regressors if it contains one of the label names. 
# If so, assign the label name
label_df.Labels = np.where(label_df.Regressors.str.contains('Magic'),'Magic',label_df.Labels)
label_df.Labels = np.where(label_df.Regressors.str.contains('Control'),'NoMagic',label_df.Labels)
label_df.Labels = np.where(label_df.Regressors.str.contains('Surprise'),'NoMagic',label_df.Labels)

VideoNames  = [''.join(re.search(" (.+?)\*",string).group(1)) 
            for string in label_df.Regressors]
label_df['VideoNames'] = VideoNames

# depending on whether we want to decode over objects or trick versions the 
# chunks change
if OVER == 'objects':
    label_df['Chunks']  = (label_df.Runs-1)//4
    chunks              = np.asarray(label_df.Chunks)
    ps = PredefinedSplit(chunks)
elif OVER == 'effects':
    # decode over effects resulting in a 3 fold cross classification
    # train on the videos of two effects and test on the remaining third
    # in order to have the same number of datapoints in each label we add 
    # surprise video 1 to effect 1, surprise 2 to effect 2 and surprise 3 to 
    # effect 3
    label_df['Effects']  = np.nan                # Labels
    # Check for every entry in Regressors if it contains one of the label names. 
    for e,effect in enumerate(EFFECT_NAMES):
        label_df.Effects = np.where(label_df.Regressors.str.contains(effect),e,label_df.Effects)
        label_df.Effects = np.where(label_df.Regressors.str.contains('Surprise'+str(e+1)),e,label_df.Effects)
    
    chunks = np.asarray(label_df.Effects)
    ps = PredefinedSplit(chunks)

elif OVER == 'tricks':
    # training on all tricks of verion 1 and test on all tricks of version 2 
    # caused a shifted null distribution below chance, hence resulting in a 
    # below chance max-statistic null distribution for multiple comparison. 
    # To resolve the problem we train on tricks of version 1 in odd trials 
    # and test on tricks of version 2 in even trials. 
    # DISADVANTAGE: less data
    # ADVANTAGE: four, instead of two validation folds (version 1/2 x odd/even)
    
    # stange and not beautiful solution (that I don't really understand) taken
    # from here:
    # (https://stackoverflow.com/questions/28837633/pandas-get-position-of-a-given-index-in-dataframe)
    
    label_df['Version'] = (label_df.VideoNames.str.contains('1')) 
    train1  = label_df.index.get_indexer_for((label_df[(label_df.Runs%2==0) & (label_df.Version==0)].index))
    test1   = label_df.index.get_indexer_for((label_df[(label_df.Runs%2==1) & (label_df.Version==1)].index))
    train2  = label_df.index.get_indexer_for((label_df[(label_df.Runs%2==1) & (label_df.Version==0)].index))
    test2   = label_df.index.get_indexer_for((label_df[(label_df.Runs%2==0) & (label_df.Version==1)].index))
    train3  = label_df.index.get_indexer_for((label_df[(label_df.Runs%2==0) & (label_df.Version==1)].index))
    test3   = label_df.index.get_indexer_for((label_df[(label_df.Runs%2==1) & (label_df.Version==0)].index))
    train4  = label_df.index.get_indexer_for((label_df[(label_df.Runs%2==1) & (label_df.Version==1)].index))
    test4   = label_df.index.get_indexer_for((label_df[(label_df.Runs%2==0) & (label_df.Version==0)].index))
    ps = [[train1, test1], [train2,test2],
          [train3, test3], [train4,test4]]
    chunks = np.asarray(label_df.Version)    # get chunks for cross validation as numpy array from data frame
else:
    raise

# THIRD STEP
#######################
# THE ACTUAL DECODING #
#######################
targets                 = np.asarray(label_df.Labels)   # get labels as numpy array from pandas dataframe
runs_for_permutation    = np.asarray(label_df.Runs)

# initialize the searchlight decoding object
MY_SEARCH_LIGHT = SearchLight(mask_img=MASK_DIR,
                             radius=SEARCHLIG_RAD,
                             estimator=my_decoder,
                             n_jobs=N_PROC,
                             cv=ps,
                             verbose=3)
betas = smooth_img(FLA_DIR+ '/' + SUB +'/'+label_df.BetaNames,fwhm=None)
# fit the decoding object based on the previously loaded betas and labes
# perform cross validation over objects
MY_SEARCH_LIGHT.fit(imgs=betas,y=targets,groups=chunks)

# Form results into a NIfTI 
results = new_img_like(ref_niimg=MASK_DIR,data=MY_SEARCH_LIGHT.scores_)
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
with open(os.path.join(RESULTS_DIR,'serchlight-logfile.txt'), 'w+') as writer:
    writer.write('Codeversion: {} \n'.format(git_hash))
    writer.write('Decoder used: {}\n'.format(DECODER))
    writer.write('Smoothing kernel of data: {}\n'.format(str(SMOOTHING_SIZE)))
    writer.write('Number of permutations: {}\n'.format(str(N_PERMS)))
    writer.write('Searchlight radius: {}\n'.format(str(SEARCHLIG_RAD)))
    writer.write('Number of kernels used: {}\n'.format(str(N_PROC)))
    writer.write('Time for computation: {}h'.format(str((time.time() - T_START)/3600)))