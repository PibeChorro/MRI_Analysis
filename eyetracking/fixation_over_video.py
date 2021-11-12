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
# The purpose of this script is plot subject fixations over an example video
# One has to choose the subject and the video

# FIRST STEP
# interact with the operating system 
import os
from pathlib import Path
import git
import glob
# data structuration and calculations
import pandas as pd  # to create data frames
import numpy as np   # most important numerical calculations
# needed to extract the run number out of the parentesis of the string in the SPM.mat file
import re
# video processing
import cv2


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
SUBJECS = glob.glob(os.path.join(DATA_DIR, 'sub*'))
SUBJECS.sort()

all_runs = np.linspace(1,12,12,dtype=int)
pre_post_runs = (all_runs-1)//2
pre_runs = all_runs[pre_post_runs%2==0]
post_runs = all_runs[pre_post_runs%2==1]

# The runs in which videos with _F as suffix were flipped. In the rest of the 
# runs the videos WITHOUT _F as suffix were flipped
flipping_runs = [1,4,5,8,9,12]

# The width we used for video presentation (we cut off some parts of the frame
# left and right)
pres_frame_width = 1600

# choose a video
vid = 'Stick_Vanish2_Magic'

# choose a subjects
s = 5
sub = os.path.basename(SUBJECS[s])

# get event and ET files of chosen subject
event_files = glob.glob(os.path.join(RAW_DIR,sub,'func','*_events.tsv'))
ET_files = glob.glob(os.path.join(DATA_DIR,sub,'*_recording-eyetracking_physio_preprocessed.tsv'))

event_files.sort()
ET_files.sort()

# if we do not have the same number of event and ET files, the subject might
# miss some ET files (happened to the first few subjects) we stop here!
if len(event_files) != len(ET_files):
    print (sub +' is missing some files')
    raise
    
# empty arrays saving x and y coordinates of fixations
positions_x = []
positions_y = []
for event, ET in zip(event_files, ET_files):
    # get the run number out of the event file name (over complicated way, but works)
    run = list(map(int, re.findall(r'\d+', os.path.basename(event))))[1]
    
    # read in current event and ET data frame
    event_df = pd.read_csv(event,sep='\t')
    ET_data = pd.read_csv(ET, sep='\t')
    
    # find indices of trials that showed the chosen video
    indices_of_interest = [i for i,trial in enumerate(event_df.trial_type) if vid in trial]
    
    for index in indices_of_interest:
        # get x and y position during the trial
        pos_x = ET_data.X_Coord[(ET_data.TimeStamp>=event_df.onset[index]) &
                               (ET_data.TimeStamp<=event_df.rating_onset[index])].values
        pos_y = ET_data.Y_Coord[(ET_data.TimeStamp>=event_df.onset[index]) &
                               (ET_data.TimeStamp<=event_df.rating_onset[index])].values
        
        # the video was flipped if it contained _F in the name AND was in one
        # of the flip runs OR when it did NOT contain _F in the name AND was
        # NOT in the flip runs
        flip = ('_F' in event_df.trial_type[index] and run in flipping_runs) or (
            '_F' not in event_df.trial_type[index] and run not in flipping_runs)
        if flip:
            print('flix x-axis')
            pos_x = abs(pos_x - pres_frame_width)
        positions_x.append(pos_x)
        positions_y.append(pos_y)
        
##################
# VIDEO CREATION #
##################
# define video codec - mp4v for mp4 video format
fourcc = cv2.VideoWriter_fourcc(*'mp4v')
# read in video and get fps and framesize
cap = cv2.VideoCapture(os.path.join(STIM_DIR,vid+'.mp4'))
fps = cap.get(cv2.CAP_PROP_FPS)
frame_width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
frame_heigt = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
# create video writer with defined codec, fps and framesize
out = cv2.VideoWriter('output.mp4', fourcc, fps, (frame_width,frame_heigt))
idx = 0
# got through video frame by frame and change pixel values of fixations 
# (+- 5 pixels)
while True:
    ret,frame = cap.read()
    if not ret:
        break
        
    for run_idx, (p_x, p_y) in enumerate(zip(positions_x,positions_y)):
        if np.isnan(p_x[idx]) or np.isnan(p_y[idx]):
            continue
        x = int(p_x[idx]*frame_width/pres_frame_width)
        y = int(p_y[idx])
        if x<0: x=0
        if x>frame_width: x=frame_width
        if y<0: y=0
        if y>frame_heigt: y=frame_heigt
        # pre runs are green
        if run_idx<3:
            frame[y-5:y+5,x-5:x+5,:] = [255,0,0]
        # post runs are blue
        else:
            frame[y-5:y+5,x-5:x+5,:] = [0,255,0]
    idx+=1
    out.write(frame)
    
cap.release()
out.release()
cv2.destroyAllWindows()

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
with open(os.path.join('fixation_video-log.txt'), 'w+') as writer:
    writer.write('Codeversion: {} \n'.format(git_hash))