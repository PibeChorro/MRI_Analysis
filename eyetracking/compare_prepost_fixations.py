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
# in MATLAB. 

# FIRST STEP
# interact with the operating system 
import os
from pathlib import Path
import git
import glob
# data structuration and calculations
import pandas as pd  # to create data frames
import numpy as np   # most important numerical calculations
import pingouin as pg
# read in mat files
import readmat
# needed to extract the run number out of the parentesis of the string in the SPM.mat file
import re
# plotting
import matplotlib.pyplot as plt

################################################
# VARIABLES FOR PATH SELECTION AND DATA ACCESS #
################################################
HOME            = str(Path.home())
# DATA
STIM_DIR = os.path.join(HOME,'Documents/Master_Thesis/Stimuli')
PROJ_DIR        = os.path.join(HOME, 'Documents/Magic_fMRI/DATA/MRI')
RAW_DIR         = os.path.join(PROJ_DIR, 'rawdata')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')
DATA_DIR = os.path.join(DERIVATIVES_DIR, 'eyetracking')
RESULTS_DIR = os.path.join(DATA_DIR, 'group_analysis')
if not os.path.isdir(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)
SUBJECS = glob.glob(os.path.join(DATA_DIR, 'sub*'))
SUBJECS.sort()
# ANALYSIS
MAGIC_MOMENT_MAT = os.path.join(HOME,'Documents/Master_Thesis/MRI_Analysis/glm/info_MagicMoment.mat')

info_magic_moment = readmat.load(MAGIC_MOMENT_MAT,isStruct=True)['do']
videos = [v for v in info_magic_moment['ListOfVideos'] if 'Magic' in v and 'excluded' not in v]
special_frame = np.array(info_magic_moment['all_frames_of_effect'])
special_location = np.array(info_magic_moment['all_pixels_of_effect'])

fps = 25
spf = 1/fps

all_runs = np.linspace(1,12,12,dtype=int)
pre_post_runs = (all_runs-1)//2
pre_runs = all_runs[pre_post_runs%2==0]
post_runs = all_runs[pre_post_runs%2==1]

flipping_runs = [1,4,5,8,9,12]
ps = []

slack_time = 0.25
pres_frame_width = 1600

for vid in videos:
    matrices = []
    fisher_matrices = []
    # loop over subjects
    for s, sub in enumerate(SUBJECS):
        event_files = glob.glob(os.path.join(RAW_DIR,os.path.basename(sub),'func','*_events.tsv'))
        ET_files = glob.glob(os.path.join(DATA_DIR,os.path.basename(sub),'*_recording-eyetracking_physio_preprocessed.tsv'))
    
        event_files.sort()
        ET_files.sort()
    
        if len(event_files) != len(ET_files):
            #print (sub +' is missing some files')
            continue
    
        positions_x = []
        positions_y = []
        for event, ET in zip(event_files, ET_files):
            run = list(map(int, re.findall(r'\d+', os.path.basename(event))))[1]
    
            event_df = pd.read_csv(event,sep='\t')
            ET_data = pd.read_csv(ET, sep='\t')
    
            indices_of_interest = [i for i,trial in enumerate(event_df.trial_type) if vid in trial]
            for index in indices_of_interest:
                pos_x = ET_data.X_Coord[(ET_data.TimeStamp>=event_df.onset[index]) &
                                       (ET_data.TimeStamp<=event_df.rating_onset[index])].values
                pos_y = ET_data.Y_Coord[(ET_data.TimeStamp>=event_df.onset[index]) &
                                       (ET_data.TimeStamp<=event_df.rating_onset[index])].values
                pos = ET_data[['X_Coord','Y_Coord']][(ET_data.TimeStamp>=event_df.onset[index]) &
                                       (ET_data.TimeStamp<=event_df.rating_onset[index])].values
    
                flip = ('_F' in event_df.trial_type[index] and run in flipping_runs) or (
                    '_F' not in event_df.trial_type[index] and run not in flipping_runs)
                if flip:
                    pos_x = abs(pos_x - pres_frame_width)
                positions_x.append(pos_x)
                positions_y.append(pos_y)
    
        df_x=pd.DataFrame(list(map(np.ravel, positions_x))).T
        df_y=pd.DataFrame(list(map(np.ravel, positions_y))).T
    
        corr_matrix = (df_x.corr()+df_y.corr())/2
        # do fisher transformation and append the correlation matrix to the array
        matrices.append(corr_matrix)
        fisher_matrices.append(np.arctanh(corr_matrix))
        
    fig = plt.figure(figsize=(15,15))
    plt.imshow(np.mean(fisher_matrices,axis=0), cmap='hot', vmin=0, vmax=2)
    plt.colorbar()
    plt.savefig(os.path.join(RESULTS_DIR,vid+'.png'))
    
    pre_fisher = [np.mean([fisher_matrices[f].iloc[0,1],
                           fisher_matrices[f].iloc[0,2],
                           fisher_matrices[f].iloc[0,3],
                           fisher_matrices[f].iloc[1,2],
                           fisher_matrices[f].iloc[1,3],
                           fisher_matrices[f].iloc[2,3]]) 
                  for f in range(len(fisher_matrices))]
    
    post_fisher = [np.mean([fisher_matrices[f].iloc[4,5],
                            fisher_matrices[f].iloc[4,6],
                            fisher_matrices[f].iloc[4,7],
                            fisher_matrices[f].iloc[5,6],
                            fisher_matrices[f].iloc[5,7],
                            fisher_matrices[f].iloc[6,7]]) 
                   for f in range(len(fisher_matrices))]
    
    ttest_res = pg.ttest(pre_fisher,post_fisher,paired=True)
    
    print ('{}: p={}'.format(vid,ttest_res['p-val'].values[0]))
    ps.append(ttest_res['p-val'].values[0])
        
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
with open(os.path.join('fixation_comparison-log.txt'), 'w+') as writer:
    writer.write('Codeversion: {} \n'.format(git_hash))