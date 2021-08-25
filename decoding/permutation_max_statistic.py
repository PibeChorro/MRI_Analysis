#!/gpfs01/bartels/user/vplikat/anaconda3/bin/python

#############
# LIBRARIES #
#############

# interact with the operating system 
import os
import argparse
from pathlib import Path
import glob
# import/export data
import h5py
# data structuration and calculations
import numpy as np   # most important numerical calculations
# optimize time performance
import time
# plotting
import matplotlib.pyplot as plt

# get start time
T_START = time.time()

################################
# Command line input arguments #
################################

# create a parser object which handles the input variables
parser = argparse.ArgumentParser()

# add all the input arguments
parser.add_argument("--runs",       "-r",   nargs="?",  const='pre',    default='pre',  type=str)
# parse the arguments to a parse-list(???)
ARGS = parser.parse_args()
# assign values 
RUNS_TO_USE     = ARGS.runs

#############################
# ALL IMPORTANT DIRECTORIES #
#############################
HOME            = str(Path.home())
PROJ_DIR        = os.path.join(HOME, 'Documents/Master_Thesis/DATA/MRI')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')

RAWDATA_DIR     = os.path.join(PROJ_DIR, 'rawdata')
ANALYSIS        = 'ROI-analysis'
DATA_DIR        = os.path.join(DERIVATIVES_DIR, 'decoding', 'decoding_magic', 
                               'decode_effect_on_'+RUNS_TO_USE+'magic',
                               'SpecialMoment','ROI-analysis','LDA')
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
        'FEF', 'IPS',
        'ACC', 'PCC', 
        'IFG', 'aINSULA', 
        'IFJ', 'PHT', 'PF'
      ]

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
    fig = plt.figure()

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
    fig.savefig(os.path.join(RESULTS_DIR,roi+'_sub_null_distributions.png'))
    
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

with h5py.File(os.path.join(RESULTS_DIR,'permutation-max-statistic.hdf5'), 'w') as f:
    f.create_dataset('mean_accuracy', data=mean_accuracy)
    f.create_dataset('max_null_distribution', data=max_statistic_null_distribution)
    f.create_dataset('p_values', data=ps)
    
fig = plt.figure()
plt.bar(range(len(ps)), ps)
plt.ylim((0,0.1))
plt.axhline(0.05,color='green')
plt.xticks(np.arange(len(ROIS)),ROIS,rotation=45)
fig.savefig(os.path.join(RESULTS_DIR,'p_values.png'))

fig = plt.figure()
plt.hist(max_statistic_null_distribution,bins=30)
plt.axvline(1/3)
fig.savefig(os.path.join(RESULTS_DIR,'max_statistic_null_distribution.png'))