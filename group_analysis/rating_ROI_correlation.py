#!/gpfs01/bartels/user/vplikat/anaconda3/bin/python

##########
# HEADER #
##########
# This script is the first one to run before doing any ROI based univariate
# analysis!
# The data for this script is derived from an fMRI experiment in which subjects
# viewed different videos. Either magic, control or surprise videos. 
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
# in MATLAB. 

##########################
# PURPOSE OF THIS SCRIPT #
##########################
# The purpose of this script is to read in the data frame containing beta 
# values of a set of ROIs and surprise ratings (created by 'create_Rating_Betas_DF.py').
# and perform a set of analyses on these data:
# Test if there are main effects of the different experimental conditions
# --> Pre-Post, Magic-Control-Surprise, Object (or Run), Effect
# Test for interaction effects (how to interprete)

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
from scipy import stats
import pingouin as pg
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
parser.add_argument("--smooth",     nargs='?', const=0,         
                    default=0,          type=int)
parser.add_argument("--analyzed",   nargs='?', const='moment',  
                    default='moment',   type=str)
# parse the arguments to a parse-list(???)
ARGS = parser.parse_args()
# assign values 
SMOOTHING_SIZE  = ARGS.smooth
ANALYZED        = ARGS.analyzed

if ANALYZED == 'moment':
    data_analyzed = 'SpecialMoment'
elif ANALYZED == 'video':
    data_analyzed = 'WholeVideo'
else:
    raise

# variables for path selection and data access
HOME            = str(Path.home())
PROJ_DIR        = os.path.join(HOME, 'Documents/Magic_fMRI/DATA/MRI')
RAWDATA_DIR     = os.path.join(PROJ_DIR, 'rawdata')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')
DATA_DIR        = os.path.join(DERIVATIVES_DIR, 'univariate-ROI',
                               data_analyzed,'EveryVideo')

data_frame = pd.read_csv(filepath_or_buffer=os.path.join(DATA_DIR,'ratings_df.csv'))

# define ROIs
ROIS = [
        'V1', 'V2', 'V3', 'hV4', 
        'V3A', 'V3B', 
        'LO', 'VO', 
        'FEF', 'IPS',
        'ACC', 'PCC', 
        'IFG', 'aINSULA', 
        'IFJ', 'PHT', 'PF'
      ]

VIDEO_TYPE = [
    'Magic',
    'Control',
    'Surprise'
    ]

EFFECTS = [
    'Appear',
    'Change',
    'Vanish'
    ]

############
# ANALYSES #
############
subject_ids = data_frame.ids.unique()
for roi in ROIS:
    sub_correlations = []
    
    for sub in subject_ids:
        ratings = data_frame.Ratings[data_frame.ids == sub].values
        betas   = data_frame[roi][data_frame.ids == sub].values
        
        [r, p] = stats.spearmanr(a=ratings, b=betas, nan_policy='omit')
        #print ('Pearson correlation between rating and {} activity r={} (p={})'.format(roi,r,p))
        
        sub_correlations.append(r)
    
    # fisher z-transform correlations
    sub_fisher_correlations = np.arctanh(sub_correlations)
    res=pg.ttest(x=sub_fisher_correlations,y=0)
    p_val = res['p-val'].values[0]
    t_val = res['T'].values[0]
    
    print ('Spearman correlation between rating and {} tested against 0: t={} (p={})'.format(roi,t_val,p_val))   
    
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
with open(os.path.join(DATA_DIR,'rating_ROI_correlation-logfile.txt'), 'w+') as writer:
    writer.write('Codeversion: {} \n'.format(git_hash))
    writer.write('Time for computation: {}h'.format(str((time.time() - T_START)/3600)))