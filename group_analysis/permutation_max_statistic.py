#!/gpfs01/bartels/user/vplikat/anaconda3/bin/python

##########
# HEADER #
##########
# Befor you run this script you have to have previously run any of the ROI 
# based decoding scripts.
# The underlying rationale and methods for the analysis in this script are 
# explained in detail in Nichols & Holmes (2002)
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
# The purpose of this script is to create one group max statistic based on the
# single subject ROI decoding results. This allows for a multiple comparison 
# corrected p-value for the ROIs.

################################
# FUNCTIONALITY OF THIS SCRIPT #
################################
# FIRST STEP 
# Get all important path information and names of ROIs.
# SECOND STEP 
# Iterate over ROIs and read in the null distribution for the current ROI from 
# all subjects. From this huge 2D array (num_permutations x num_subjects) 
# calculate the mean over subjects (=1D array with length num_permutations).
# Cast together all mean null distributions (resulting again in a 2D array 
# num_permutations x num_ROIs) and get the max value for each permutation along
# the ROIs (=1D array with length num_permutations).
# The mean accuracy of the ROIs over subjects is also calculated and saved.
# This mean is later tested against the max statistic distribution.
# In between all null distributions of all subjects in one ROI are plotted in 
# one figure to make sure that it is centered around chance level. 
# THIRD STEP
# Calculate corrected p-value for each ROI: Look how many accuracy values in 
# the max statistic distribution are higher than the observed mean accuracy
# value (called C)
# p_corrected = (C+1)/(num_permutation+1)
# FOURTH STEP
# Plot the p-value for each ROI (first figure)
# For each ROI calculate the 95% confidence interval and plot mean + CI next
# to the max statistic. 
# Plot the max statistic

######################
# COMMAND LINE FLAGS #
######################
# --data: data from which previously calculated permutation tests should be 
# used. Accepted values are pre (magic effect decoding pre revelation), post
# (magic effect decoding post revelation), all (magic effect decoding using 
# pre and post revelation data), pre-post (decoding videos pre vs. post 
# revelation) and mag-nomag (decoding magic vs nomagic).

#############
# LIBRARIES #
#############

# interact with the operating system 
import os
import argparse
from pathlib import Path
import glob
import git
# import/export data
import h5py
# data structuration and calculations
import numpy as np   # most important numerical calculations
import pandas as pd
import scipy.stats as st
# optimize time performance
import time
# plotting
import matplotlib.pyplot as plt

def mean_confidence_interval(data, confidence=0.95):
    '''mean_confidence_interval: A function that calculates the confidence
    interval (CI) based on the normal distribution using numpy and scipy.stats
    data: array like object
    confidence: float between 0 and 1
    returns: mean and CI of data'''
    m, se = np.mean(data), st.sem(data)
    h = se * st.norm.ppf((1 + confidence) / 2.)
    return m, m-h, m+h

# get start time
T_START = time.time()

# FIRST STEP
################################
# Command line input arguments #
################################

# create a parser object which handles the input variables
parser = argparse.ArgumentParser()

# add all the input arguments
parser.add_argument("--data", "-d", nargs="?", const='pre', default='pre', 
                    type=str)
parser.add_argument("--over",  "-o",   nargs='?',  const='objects', 
                    default='objects')
# parse the arguments to a parse-list(???)
ARGS = parser.parse_args()
# assign values 
DATA = ARGS.data
OVER = ARGS.over

#############################
# ALL IMPORTANT DIRECTORIES #
#############################
HOME            = str(Path.home())
PROJ_DIR        = os.path.join(HOME, 'Documents/Magic_fMRI/DATA/MRI')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')

RAWDATA_DIR     = os.path.join(PROJ_DIR, 'rawdata')
ANALYSIS        = 'ROI-analysis'
if DATA == 'pre':
    DATA_TO_USE = 'decode_effect_on_premagic'
    NUM_LABELS  = 3
elif DATA == 'post':
    DATA_TO_USE = 'decode_effect_on_postmagic'
    NUM_LABELS  = 3
elif DATA == 'all':
    DATA_TO_USE = 'decode_effect_on_allmagic'
    NUM_LABELS  = 3
elif DATA == 'pre-post':
    DATA_TO_USE = 'decode_pre_vs_post'
    NUM_LABELS  = 2
elif DATA == 'mag-nomag':
    DATA_TO_USE = 'magic_vs_nomagic'
    NUM_LABELS  = 2
else:
    raise
DATA_DIR        = os.path.join(DERIVATIVES_DIR, 'decoding', 'decoding_magic', 
                               DATA_TO_USE, 'over_' + OVER, 'SpecialMoment',
                               'ROI-analysis')
RESULTS_DIR     = os.path.join(DATA_DIR,'group-statistics')
if not os.path.isdir(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

SUBJECTS = glob.glob(os.path.join(DATA_DIR,'sub-*'))
SUBJECTS.sort()

# define ROIs
ROIS = [
        'V1', 'V2', 'V3', 'hV4', 
        'V3A', 'V3B', 
        'LO', 'VO', 
        'PHT', 'PF', 
        'FEF', 'IPS',
        'ACC', 'PCC', 
        'IFG', 'aINSULA', 
        'IFJ'
      ]

# SECOND STEP 
# empty list, that stores the full accuracy arrays per ROI, so that a CI can 
# be calculated
roi_accuracies = []

# empty lists where the mean accuracies and null distributions for ROIs will 
# be stored
accuracy_mean_over_subs             = []
null_distribution_mean_over_subs    = []

for r, roi in enumerate(ROIS):
    # empty lists for accuracies and null distributions of subjects for the
    # current ROI
    accuracies          = []
    null_distributions  = []
    
    # create a figure. Here we plot all null distributions to make sure they
    # are symetrical and centured around chance level
    roi_fig = plt.figure()

    # inner loop - iterating over mask (=ROIs)
    for s, sub in enumerate(SUBJECTS):
        # read in the hdf5 data file for ROI[r] and SUBJECT[s]
        roi_file = os.path.join(sub,roi + '.hdf5')
        res = h5py.File(roi_file, 'r')
        
        # read out the accuracy and null distribution
        accuracies.append(res['accuracy'][()])
        null_distributions.append(res['null_distribution'][()])
        
        # plot null distribution of subject for current ROI
        plt.hist(res['null_distribution'][()],
                  bins=30,
                  color='blue',
                  alpha=1/len(SUBJECTS))
        
    # append the accuracies list to the roi_accuracies list to later plot their
    # confidence intervals
    roi_accuracies.append(accuracies)
    # calculate mean accuracy and null distribution of current ROI over 
    # subjects. Append both to outer lists
    mean_accuracy           = np.mean(accuracies)
    mean_null_distribution  = np.mean(null_distributions,axis = 0)
    accuracy_mean_over_subs.append(mean_accuracy)
    null_distribution_mean_over_subs.append(mean_null_distribution)
    
    # plot mean null distribution for ROI over subjects
    plt.hist(mean_null_distribution,
             bins=30,
             color='red',
             alpha=0.2)
    plt.axvline(1/3)
    roi_fig.savefig(os.path.join(RESULTS_DIR,roi+'_sub_null_distributions.png'))
    
# create a dataframe containing the accuracies of every subject for every roi
accuracies_df = pd.DataFrame(list(map(np.ravel, roi_accuracies))).T
# set the columns correct
accuracies_df.columns = ROIS
accuracies_df.to_csv(os.path.join(RESULTS_DIR, 'accuracies.csv'))
    
# THIRD STEP
# after getting all mean null distributions get max-null statistic
null_distribution_mean_over_subs    = np.asarray(null_distribution_mean_over_subs)
max_statistic_null_distribution     = null_distribution_mean_over_subs.max(axis=0)

# calculate p-values for each ROI by getting the number of accuracies larger
# than the 'real' accuracy
ps = []
for acc in accuracy_mean_over_subs:
    ps.append(sum(max_statistic_null_distribution>acc))
    
ps = np.asarray(ps)
ps = (ps+1)/(max_statistic_null_distribution.shape[0]+1)

result_dict = {
    'ROIs': ROIS,
    'mean_accuracies': accuracy_mean_over_subs,
    'p_values': ps
    }
    
# FOURTH STEP
# plot confidence intervalls for all ROIs and max statistic distribution in
# one figure and draw significance threshold and theoretical chance level
CI_fig, ax= plt.subplots(1,2,sharey=True,figsize=(20,10))
ax[0].hist(max_statistic_null_distribution,
           bins=50, 
           orientation='horizontal',
           label='max null distribution')
# draw the significance threshold (alpha=0.05 -> 95% percentile)
ax[0].axhline(np.percentile(max_statistic_null_distribution,95),
              color='green', 
              label='alpha=0.05 threshold')
ax[1].axhline(np.percentile(max_statistic_null_distribution,95),
              color='green', 
              label='alpha=0.05 threshold')

# draw the chance level
ax[0].axhline(1/NUM_LABELS,
              color='red',
              linestyle='--', 
              label='Chance level')
ax[1].axhline(1/NUM_LABELS,
              color='red',
              linestyle='--', 
              label='Chance level')

# empty lists for CI-bounds
lows = []
highs = []
# for all ROIs plot the mean and CI 
for a, accs in enumerate(roi_accuracies):
    mean, low, high = mean_confidence_interval(accs)
    lows.append(low)
    highs.append(high)
    ax[1].plot([a,a],[low,high],'b')
    ax[1].plot(a,mean,'b*')

# figure make up
plt.xticks(np.arange(len(ROIS)),ROIS,rotation=45)
ax[0].legend()
ax[1].legend()
CI_fig.savefig(os.path.join(RESULTS_DIR,'accuracy_CI.png'))    

# plot max statistic null distribution
max_stat_fig = plt.figure()
plt.hist(max_statistic_null_distribution,bins=30)
plt.axvline(1/NUM_LABELS)
max_stat_fig.savefig(os.path.join(RESULTS_DIR,'max_statistic_null_distribution.png'))

# save the mean accuracies for the ROIs, the max statistic null distribution
# and the p-values for decoding accuracies for the ROIs
result_dict['CI_lower'] = lows
result_dict['CI_higher'] = highs

result_df = pd.DataFrame(data=result_dict,columns=result_dict.keys())
result_df.to_csv(path_or_buf=os.path.join(RESULTS_DIR,'max-statistic-results.csv'),
                 index=False)
with h5py.File(os.path.join(RESULTS_DIR,'permutation-max-statistic.hdf5'), 'w') as f:
    f.create_dataset('max_null_distribution', data=max_statistic_null_distribution)
    f.create_dataset('crit_acc_value', data=np.quantile(max_statistic_null_distribution,0.95))

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
with open(os.path.join(RESULTS_DIR,'max_stat-logfile.txt'), 'w+') as writer:
    writer.write('Codeversion: {} \n'.format(git_hash))
    writer.write('Time for computation: {}h'.format(str((time.time() - T_START)/3600)))
