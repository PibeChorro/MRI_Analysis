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
# The purpose of this script is to preprocess the eyetracking data and save 
# the created clean data in a seperate location

################################
# FUNCTIONALITY OF THIS SCRIPT #
################################
# Additionally move the edf files to sourcedata
# FIRST STEP 
# Import all libraries needed and get all important path information
# SECOND STEP 
# Eyetracking data preprocessing:
# a) remove data before and after blinks
# b) interpolate missing data (blinks and other)
# c) low pass and high pass filter of diameter
# d) demean pupil diameter per condition
# e) divide by standard deviation per session (run in this case)
# f) resample - all data 
# THIRD STEP
# Save preprocessed data in 'derivatives' folder

# FIRST STEP
# interact with the operating system 
import os
import argparse
from pathlib import Path
import git
import glob
# data structuration and calculations
import pandas as pd  # to create data frames
import numpy as np   # most important numerical calculations
from scipy import signal
# needed to extract the run number out of the parentesis of the string in the SPM.mat file
import re
from my_ET_functions import eyedata2pandasframe

# create a parser object which handles the input variables
parser = argparse.ArgumentParser()

# add all the input arguments
parser.add_argument("--sub", "-s", default='sub-01')         # subject

# parse the arguments to a parse-list(???)
ARGS = parser.parse_args()
# assign values 
SUB             = ARGS.sub

################################################
# VARIABLES FOR PATH SELECTION AND DATA ACCESS #
################################################
HOME            = str(Path.home())
# DATA
PROJ_DIR        = os.path.join(HOME, 'Documents/Magic_fMRI/DATA/MRI')
RAW_DIR         = os.path.join(PROJ_DIR, 'rawdata')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')
RESULT_DIR = os.path.join(DERIVATIVES_DIR, 'eyetracking')

# ANALYSIS 
fps = 25    # frames per second of the videos
spf = 1/fps # inverse of fps

# ET preprocessing following Knappen et al 2016
# removing data from 0.15 seconds before and after blinks
margin              = 0.15
interpolation_limit = margin*2 + 0.2

#filter
butter_order    = 3
high_pass_freq  = 0.02
low_pass_freq   = 4

# raw unprocessed ET data files
ET_files = glob.glob(os.path.join(RAW_DIR, SUB, 'func',
                                  '*_recording-eyetracking_physio.asc'))
ET_files.sort()

# event files containing trial onset, duration video name and surprise response

for ET in ET_files:
    # get run number. We do not want to enumerate ET_files because some subs
    # miss some data files
    run = list(map(int, re.findall(r'\d+', os.path.basename(ET))))[1]
    
    # read in the coresponding event file depending on run number, because
    # some subjects miss some ET data, so the nth ET file does not necessarily
    # need to be the coresponding event file
    event_file = glob.glob(os.path.join(RAW_DIR,SUB,'func',
                                        '*run-{:02d}_events.tsv'.format(run)))
    try:
        event_df = pd.read_csv(filepath_or_buffer=event_file[0],sep='\t')
    except:
        print('could not read event file')
        raise
    
    # call own function that reads out infos from the ET ascii file
    ET_data, fixations, saccades, blinks, trials, sample_rate, frame_dim = eyedata2pandasframe(ET)
    # Cut off all data, that exceeds the last trial
    ET_data = ET_data[ET_data.TimeStamp<event_df.onset.iloc[-1]+
                      event_df.duration.iloc[-1]]
    frame_width = frame_dim[2]
    frame_height = frame_dim[3]
    
    # set all ET data values during blinks (+- 150ms) to np.nan
    for b_s,b_e in zip(blinks.Start,blinks.End):
        ET_data.Diameter[(ET_data.TimeStamp>b_s-margin) & (ET_data.TimeStamp<b_e+margin)] = np.nan
        ET_data.X_Coord[(ET_data.TimeStamp>b_s-margin) & (ET_data.TimeStamp<b_e+margin)] = np.nan
        ET_data.Y_Coord[(ET_data.TimeStamp>b_s-margin) & (ET_data.TimeStamp<b_e+margin)] = np.nan

    # now interpolate them
    ET_data.Diameter.interpolate(method ='linear',
                                 limit=int(interpolation_limit*sample_rate), 
                                 limit_direction='both',
                                 inplace=True)
    ET_data.X_Coord.interpolate(method ='linear',
                                 limit=int(interpolation_limit*sample_rate), 
                                 limit_direction='both',
                                 inplace=True)
    ET_data.Y_Coord.interpolate(method ='linear',
                                 limit=int(interpolation_limit*sample_rate), 
                                 limit_direction='both',
                                 inplace=True)
    
    # filter stuff
    high_pass_cof_sample    = high_pass_freq / (sample_rate/2)
    low_pass_cof_sample     = low_pass_freq / (sample_rate/2)
    
    bbp, abp = signal.butter(butter_order, 
                             [high_pass_cof_sample, low_pass_cof_sample], 
                             btype='bandpass')

    diameter_filtered = []
    
    # iterate over trials and filter data within a trial
    for idx,tr in event_df.iterrows():
        tmp_diameter = ET_data.Diameter[(ET_data.TimeStamp>=tr.onset) & 
                                        (ET_data.TimeStamp<tr.onset+tr.duration)]
        tmp_diameter_bp = signal.filtfilt(bbp, abp, tmp_diameter)
        diameter_filtered.extend(tmp_diameter_bp)

    ET_data['Diameter_filtered'] = diameter_filtered
    
    # normalization - mean over trials, standard deviation over session (run)
    diameter_std = ET_data.Diameter_filtered.std() 
    diameter_normalized = []
    for idx,tr in event_df.iterrows():
        tmp_diameter = ET_data.Diameter_filtered[(ET_data.TimeStamp>=tr.onset) & 
                                        (ET_data.TimeStamp<tr.onset+tr.duration)]
        tmp_diameter -= np.nanmean(tmp_diameter)
        tmp_diameter /= diameter_std
        diameter_normalized.extend(tmp_diameter)
    
    ET_data['Diameter_normalized'] = diameter_normalized
    
    # set to borders or to np.nan ???
    #ET_data.X_Coord[ET_data.X_Coord > frame_width] = frame_width
    #ET_data.X_Coord[ET_data.X_Coord < 0] = 0
    
    #ET_data.Y_Coord[ET_data.Y_Coord > frame_height] = frame_height
    #ET_data.Y_Coord[ET_data.Y_Coord < 0] = 0
    
    #Resample data --> one data point per frame (use mean)
    ET_data = ET_data.resample(str(spf*1000)+'ms').mean()
    
    # save preprocessed data 
    file_name = SUB+'_task-magic_run-{:02d}_recording-eyetracking_physio_preprocessed.tsv'.format(run)
    save_dir = os.path.join(RESULT_DIR, SUB)
    if not os.path.isdir(save_dir):
        os.makedirs(save_dir)
    ET_data.to_csv(path_or_buf=os.path.join(save_dir,file_name), 
                   sep='\t',
                   na_rep = 'n/a',
                   index = False)
    
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
with open(os.path.join(save_dir,'ET_preproc-logfile.txt'), 'w+') as writer:
    writer.write('Codeversion: {} \n'.format(git_hash))