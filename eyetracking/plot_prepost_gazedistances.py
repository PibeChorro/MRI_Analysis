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

##########################
# PURPOSE OF THIS SCRIPT #
##########################
# The purpose of this script is to look into the eyetracking data and check 
# if subjects fixated at the same position during the special moment pre vs post
# revelation

# FIRST STEP
# interact with the operating system 
import os
from pathlib import Path
import git
import glob
# data structuration and calculations
import pandas as pd  # to create data frames
import numpy as np   # most important numerical calculations
from scipy.spatial import distance
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
PROJ_DIR        = os.path.join(HOME, 'Documents/Magic_fMRI/DATA/MRI')
RAW_DIR         = os.path.join(PROJ_DIR, 'rawdata')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')
DATA_DIR = os.path.join(DERIVATIVES_DIR, 'eyetracking')
SUBJECS = glob.glob(os.path.join(DATA_DIR, 'sub*'))
SUBJECS.sort()
# ANALYSIS
MAGIC_MOMENT_MAT = os.path.join(HOME,'Documents/Master_Thesis/MRI_Analysis/glm/info_MagicMoment.mat')

info_magic_moment = readmat.load(MAGIC_MOMENT_MAT,isStruct=True)['do']
videos = np.array(info_magic_moment['ListOfVideos'])
special_frame = np.array(info_magic_moment['all_frames_of_effect'])
special_location = np.array(info_magic_moment['all_pixels_of_effect'])

fps = 25
spf = 1/fps

all_runs = np.linspace(1,12,12,dtype=int)
pre_post_runs = (all_runs-1)//2
pre_runs = all_runs[pre_post_runs%2==0]
post_runs = all_runs[pre_post_runs%2==1]

flipping_runs = [1,4,5,8,9,12]

slack_time = 0.25
frame_width = 1600

for s, sub in enumerate(SUBJECS):
    event_files = glob.glob(os.path.join(RAW_DIR,os.path.basename(sub),
                                         'func','*_events.tsv'))
    ET_files = glob.glob(os.path.join(sub,'*_recording-eyetracking_physio_preprocessed.tsv'))
    
    event_files.sort()
    ET_files.sort()
    
    if len(event_files) != len(ET_files):
        print (sub +' is missing some files')
        continue
        
    sub_dict = {
        'vids': [],
        'pre_revelation': [],
        'runs': [],
        'type': [],
        'x_positions': [],
        'y_positions': []
        }
    
    for event, ET in zip(event_files, ET_files):
        run = list(map(int, re.findall(r'\d+', os.path.basename(event))))[1]
        
        event_df = pd.read_csv(event,sep='\t')
        ET_data = pd.read_csv(ET, sep='\t')
        
        # iterate over videos
        for index,row in event_df.iterrows():
            vid = row.trial_type
            flip = ('_F' in vid and run in flipping_runs) or (
                '_F' not in vid and run not in flipping_runs)
            
            if '_F' in vid:
                vid = vid[:-2]
                
            special_moment = special_frame[videos==vid][0]*spf + row.onset
            
            x_pos = ET_data.X_Coord[(ET_data.TimeStamp>=special_moment-slack_time) & 
                                    (ET_data.TimeStamp<=special_moment+slack_time)].values
            y_pos = ET_data.Y_Coord[(ET_data.TimeStamp>=special_moment-slack_time) & 
                                    (ET_data.TimeStamp<=special_moment+slack_time)].values
            
            x = np.nanmean(x_pos)
            if flip:
                x = abs(x-frame_width)
    
            y = np.nanmean(y_pos)
            
            sub_dict['x_positions'].append(x)
            sub_dict['y_positions'].append(y)
            sub_dict['vids'].append(vid)
            sub_dict['pre_revelation'].append(run in pre_runs)
            sub_dict['runs'].append(run)
            if 'Magic' in vid:
                sub_dict['type'].append('Magic')
            elif 'Control' in vid:
                sub_dict['type'].append('Control')
            elif 'Surprise' in vid:
                sub_dict['type'].append('Surprise')
            else:
                raise
            
    sub_df = pd.DataFrame(sub_dict,columns=sub_dict.keys())
    
    euc_distances = []

    for vid in sub_df[sub_df.type=='Magic'].vids.unique():
        x_pre = sub_df.x_positions[(sub_df.vids==vid)&
                                   (sub_df.pre_revelation==True)].mean()
        x_post = sub_df.x_positions[(sub_df.vids==vid)&
                                    (sub_df.pre_revelation==False)].mean()
        
        y_pre = sub_df.y_positions[(sub_df.vids==vid)&
                                   (sub_df.pre_revelation==True)].mean()
        y_post = sub_df.y_positions[(sub_df.vids==vid)&
                                    (sub_df.pre_revelation==False)].mean()
        
        dist = distance.euclidean([x_pre,y_pre],[x_post,y_post])
        euc_distances.append(dist)
        if dist>200:
            print(os.path.basename(sub))
            print(vid)
            print(dist)
    plt.hist(euc_distances);
    
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
with open(os.path.join(DATA_DIR,'dist_plot-logfile.txt'), 'w+') as writer:
    writer.write('Codeversion: {} \n'.format(git_hash))