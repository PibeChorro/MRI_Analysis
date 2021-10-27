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
from statsmodels.stats.anova import AnovaRM
import pingouin as pg
import matplotlib.pyplot as plt
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
PROJ_DIR        = os.path.join(HOME, 'Documents/Master_Thesis/DATA/MRI')
RAWDATA_DIR     = os.path.join(PROJ_DIR, 'rawdata')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')
DATA_DIR        = os.path.join(DERIVATIVES_DIR, 'univariate-ROI',
                               data_analyzed,'EveryVideo')

data_frame = pd.read_csv(filepath_or_buffer=os.path.join(DATA_DIR,'ratings_df.csv'))

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
# repeated measures ANOVA of pre-post x Video type

na_removed = pg.remove_rm_na(data=data_frame,
                             dv='Ratings', 
                             subject='ids',
                             within=['PrePost', 'Type'],
                             aggregate='mean')

spher = pg.sphericity(data=na_removed,
                      dv='Ratings', 
                      subject='ids',
                      within=['PrePost', 'Type'])
    
aov_res = pg.rm_anova(data=na_removed,
                      dv='Ratings', 
                      subject='ids',
                      within=['PrePost', 'Type'],
                      correction = not spher.spher,
                      detailed=True)

fig = plt.figure()
for vid_type in VIDEO_TYPE:
    rating_pre = na_removed.loc[:,'Ratings'][(na_removed.Type==vid_type)&
                                             (na_removed.PrePost==0)]
    rating_post = na_removed.loc[:,'Ratings'][(na_removed.Type==vid_type)&
                                              (na_removed.PrePost==1)]
    pre_post_means  = [rating_pre.mean(),rating_post.mean()]
    pre_post_means  = np.array(pre_post_means)
    pre_post_sems   = [rating_pre.sem(),rating_post.sem()]
    pre_post_sems   = np.array(pre_post_sems)
    upper = pre_post_means+pre_post_sems
    lower = pre_post_means-pre_post_sems
    plt.plot(pre_post_means)
    plt.fill_between([0,1], upper, lower, alpha=0.2)

# # if the interaction effect is significant
# if aov_res['p-GG-corr'][2]<0.05:
#     Magic_df = na_removed[na_removed.Type == 'Magic']
#     magic_ttest = pg.ttest(x=Magic_df.Ratings[Magic_df.PrePost==0],
#                            y=Magic_df.Ratings[Magic_df.PrePost==1],
#                            paired=True, 
#                            alternative='greater')
    
#     Control_df = na_removed[na_removed.Type == 'Control']
#     control_ttest = pg.ttest(x=Control_df.Ratings[Control_df.PrePost==0],
#                              y=Control_df.Ratings[Control_df.PrePost==1],
#                              paired=True, 
#                              alternative='greater')
    
#     Surprise_df = na_removed[na_removed.Type == 'Surprise']
#     surprise_ttest = pg.ttest(x=Surprise_df.Ratings[Surprise_df.PrePost==0],
#                               y=Surprise_df.Ratings[Surprise_df.PrePost==1],
#                               paired=True, 
#                               alternative='greater')
    
# repeated measures ANOVA for effects x objects (only magic videos)
magic_df = data_frame[data_frame.Type=='Magic']
na_removed = pg.remove_rm_na(data=magic_df,
                             dv='Ratings', 
                             subject='ids',
                             within=['Objects', 'Effect'],
                             aggregate='mean')

rm_aov = AnovaRM(data=na_removed, 
                 depvar='Ratings', 
                 subject='ids',
                 within=['Objects', 'Effect'])
res = rm_aov.fit()

fig = plt.figure()
eff_values  = []
eff_sems    = []
for eff in EFFECTS:
    eff_ratings = na_removed.loc[:,'Ratings'][(na_removed.Effect==eff)]
    eff_values.append(eff_ratings.mean())
    eff_sems.append(eff_ratings.sem())
    
eff_values  = np.array(eff_values)
eff_sems    = np.array(eff_sems)
upper = eff_values+eff_sems
lower = eff_values-eff_sems
plt.plot(eff_values)
plt.fill_between([0,1,2], upper, lower, alpha=0.2)


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
with open(os.path.join(DATA_DIR,'rating_analysis-logfile.txt'), 'w+') as writer:
    writer.write('Codeversion: {} \n'.format(git_hash))
    writer.write('Time for computation: {}h'.format(str((time.time() - T_START)/3600)))