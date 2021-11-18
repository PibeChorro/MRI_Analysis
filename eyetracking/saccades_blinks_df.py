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
import argparse
from pathlib import Path
import git
import glob
# data structuration and calculations
import pandas as pd  # to create data frames
import numpy as np   # most important numerical calculations
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
RESULT_DIR      = os.path.join(DERIVATIVES_DIR, 
                               'eyetracking', 
                               'group_analysis')

subjects = glob.glob(os.path.join(RAW_DIR, 'sub*'))
subjects.sort()
subjects = [os.path.basename(sub) for sub in subjects]

data_dict = {
    'num_saccades': [],
    'num_blinks': [],
    'pre_post': [],
    'video_type': [],
    'subject': [],
    'run': []
    }

all_runs        = np.linspace(1,12,12,dtype=int)
pre_post_runs   = (all_runs-1)//2
pre_runs        = all_runs[pre_post_runs%2==0]
post_runs       = all_runs[pre_post_runs%2==1]

for s, sub in enumerate(subjects):

    # raw unprocessed ET data files
    ET_files = glob.glob(os.path.join(RAW_DIR, sub, 'func',
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
            print('Could not read event file  {}of {}'.format(run, sub))
            continue
        
        try:
            # call own function that reads out infos from the ET ascii file
            _, _, saccades, blinks, _, _, _ = eyedata2pandasframe(ET)
        except:
            print('Could not read ET file {} of {}'.format(run, sub))
            continue    
        # Cut off all data, that exceeds the last trial
        saccades = saccades[saccades.Start<event_df.onset.iloc[-1]+
                          event_df.duration.iloc[-1]]
        
        saccades['duration'] = saccades.End - saccades.Start
        
        blinks = blinks[blinks.Start<event_df.onset.iloc[-1]+
                          event_df.duration.iloc[-1]]
        
        blinks['duration'] = blinks.End - blinks.Start
        
        # calculate blink duration mean and std to exclude blinks that exceed 
        # mean+std
        blink_duration_mean = np.mean(blinks.duration)
        blink_duration_std = np.std(blinks.duration)
        
        blinks = blinks[blinks.duration<blink_duration_mean+blink_duration_std]
        
        num_blinks_mag = 0
        num_saccs_mag = 0
        num_blinks_con = 0
        num_saccs_con = 0
        num_blinks_sur = 0
        num_saccs_sur = 0
        
        for idx, row in event_df.iterrows():
            tr_blinks = len (blinks[(blinks.End>row.onset) & 
                                    (blinks.Start<row.rating_onset)])
            tr_sacs = len (saccades[(saccades.End>row.onset) & 
                                    (saccades.Start<row.rating_onset)])
            
            if 'Magic' in row.trial_type:
                num_blinks_mag += tr_blinks
                num_saccs_mag += tr_sacs
            elif 'Control' in row.trial_type:
                num_blinks_con += tr_blinks
                num_saccs_con += tr_sacs
            elif 'Surprise' in row.trial_type:
                num_blinks_sur += tr_blinks
                num_saccs_sur += tr_sacs
            else:
                raise
        
        data_dict['num_saccades'].extend([num_saccs_mag,num_saccs_con,num_saccs_sur])
        data_dict['num_blinks'].extend([num_blinks_mag,num_blinks_con,num_blinks_sur])
        data_dict['video_type'].extend(['Magic','Control','Surprise'])
        data_dict['subject'].extend([sub, sub, sub])
        data_dict['run'].extend([run, run, run])
        if run in pre_runs:
            data_dict['pre_post'].extend(['pre', 'pre', 'pre'])
        else:
            data_dict['pre_post'].extend(['post', 'post', 'post'])
        
data_df = pd.DataFrame(data=data_dict,columns=data_dict.keys())
data_df.to_csv(os.path.join(RESULT_DIR,'blink_saccade_df.csv'))

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
with open(os.path.join(RESULT_DIR,'saccade_blink_comparison-logfile.txt'), 'w+') as writer:
    writer.write('Codeversion: {} \n'.format(git_hash))