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
# data structuration and calculations
import pandas as pd  # to create data frames
import pingouin as pg

################################################
# VARIABLES FOR PATH SELECTION AND DATA ACCESS #
################################################
HOME            = str(Path.home())
# DATA
PROJ_DIR        = os.path.join(HOME, 'Documents/Magic_fMRI/DATA/MRI')
RAW_DIR         = os.path.join(PROJ_DIR, 'rawdata')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')
DATA_DIR        = os.path.join(DERIVATIVES_DIR, 
                               'eyetracking', 
                               'group_analysis')

BLINK_SAC_DF = pd.read_csv(os.path.join(DATA_DIR,'blink_saccade_df.csv'))
# normalize number of blinks and saccades by number of video presentations
BLINK_SAC_DF.num_saccades[BLINK_SAC_DF.video_type=='Magic'] /= 12
BLINK_SAC_DF.num_saccades[BLINK_SAC_DF.video_type=='Control'] /= 6
BLINK_SAC_DF.num_saccades[BLINK_SAC_DF.video_type=='Surprise'] /= 6

BLINK_SAC_DF.num_blinks[BLINK_SAC_DF.video_type=='Magic'] /= 12
BLINK_SAC_DF.num_blinks[BLINK_SAC_DF.video_type=='Control'] /= 6
BLINK_SAC_DF.num_blinks[BLINK_SAC_DF.video_type=='Surprise'] /= 6

# first remove missing values and average over conditions
aggregated_df = pg.remove_rm_na(data=BLINK_SAC_DF,
                                subject='subject',
                                within=['pre_post', 'video_type'],
                                aggregate='mean')

# remove subject 4, due to missing some files --> only pre revelation files
aggregated_df = aggregated_df[aggregated_df.subject!='sub-04']

spher = pg.sphericity(data=aggregated_df,
                      dv='num_saccades', 
                      subject='subject',
                      within=['pre_post', 'video_type'])

aov_res = pg.rm_anova(data=aggregated_df,
                      dv='num_saccades', 
                      subject='subject',
                      within=['pre_post', 'video_type'],
                      correction = not spher.spher,
                      detailed=True)

if aov_res['p-unc'][0]<0.05:
    x_var = aggregated_df.num_saccades[(aggregated_df.video_type=='Magic')&
                                       (aggregated_df.pre_post=='pre')]
    y_var = aggregated_df.num_saccades[(aggregated_df.video_type=='Magic')&
                                       (aggregated_df.pre_post=='post')]
    ttest_res_mag = pg.ttest(x=x_var, y=y_var, paired=True)
    
    x_var = aggregated_df.num_saccades[(aggregated_df.video_type=='Control')&
                                       (aggregated_df.pre_post=='pre')]
    y_var = aggregated_df.num_saccades[(aggregated_df.video_type=='Control')&
                                       (aggregated_df.pre_post=='post')]
    ttest_res_con = pg.ttest(x=x_var, y=y_var, paired=True)
    
    x_var = aggregated_df.num_saccades[(aggregated_df.video_type=='Surprise')&
                                       (aggregated_df.pre_post=='pre')]
    y_var = aggregated_df.num_saccades[(aggregated_df.video_type=='Surprise')&
                                       (aggregated_df.pre_post=='post')]
    ttest_res_sur = pg.ttest(x=x_var, y=y_var, paired=True)
    
    post_hoc_prepost_sacc_res=pd.concat([ttest_res_mag,ttest_res_con,ttest_res_sur])
    post_hoc_prepost_sacc_res['VideoType'] = ['Magic', 'Control', 'Surprise']
    post_hoc_prepost_sacc_res.to_csv(path_or_buf=os.path.join(DATA_DIR,'post_hoc_prepost_sacc_results.csv'),
                        index=False)
    
spher = pg.sphericity(data=aggregated_df,
                      dv='num_blinks', 
                      subject='subject',
                      within=['pre_post', 'video_type'])

aov_res = pg.rm_anova(data=aggregated_df,
                      dv='num_blinks', 
                      subject='subject',
                      within=['pre_post', 'video_type'],
                      correction = not spher.spher,
                      detailed=True)

if aov_res['p-unc'][0]<0.05:
    x_var = aggregated_df.num_blinks[(aggregated_df.video_type=='Magic')&
                                       (aggregated_df.pre_post=='pre')]
    y_var = aggregated_df.num_blinks[(aggregated_df.video_type=='Magic')&
                                       (aggregated_df.pre_post=='post')]
    ttest_res_mag = pg.ttest(x=x_var, y=y_var, paired=True)
    
    x_var = aggregated_df.num_blinks[(aggregated_df.video_type=='Control')&
                                       (aggregated_df.pre_post=='pre')]
    y_var = aggregated_df.num_blinks[(aggregated_df.video_type=='Control')&
                                       (aggregated_df.pre_post=='post')]
    ttest_res_con = pg.ttest(x=x_var, y=y_var, paired=True)
    
    x_var = aggregated_df.num_blinks[(aggregated_df.video_type=='Surprise')&
                                       (aggregated_df.pre_post=='pre')]
    y_var = aggregated_df.num_blinks[(aggregated_df.video_type=='Surprise')&
                                       (aggregated_df.pre_post=='post')]
    ttest_res_sur = pg.ttest(x=x_var, y=y_var, paired=True)
    
    post_hoc_prepost_blink_res=pd.concat([ttest_res_mag,ttest_res_con,ttest_res_sur])
    post_hoc_prepost_blink_res['VideoType'] = ['Magic', 'Control', 'Surprise']
    post_hoc_prepost_blink_res.to_csv(path_or_buf=os.path.join(DATA_DIR,'post_hoc_prepost_sacc_results.csv'),
                        index=False)
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