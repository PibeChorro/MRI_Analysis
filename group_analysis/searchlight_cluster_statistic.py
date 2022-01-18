#!/gpfs01/bartels/user/vplikat/anaconda3/bin/python

##########
# HEADER #
##########
# Before you run this script you have to have run the 
# 'searchlight_group_analysis.py' script.
# The underlying rationale and methods for the analysis in this script are 
# explained in detail in Stelzer, Chen & Turner (2012)
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
# in MATLAB and is in MNI space. 

##########################
# PURPOSE OF THIS SCRIPT #
##########################
# The purpose of this script is to calculate the significant cluster size in 
# a mean accuracy map derived from searchlight analysis.

################################
# FUNCTIONALITY OF THIS SCRIPT #
################################
# FIRST STEP 
# Get all important path information and information about the analysis from 
# command line input.
# SECOND STEP
# Read in the mean accuracy map of the observed values. 
# Create a "critical accuracy value map". For each voxel get the accuracy value
# that is significant. 
# Since we cannot read in all 10 000 chance accuracy maps (previously created 
# by 'searchlight_group_analysis.py') we read in slice after slice. 
# THIRD STEP
# Get the observed cluster size distribution. Therefore get the voxels which 
# have an accuracy above the threshold (voxel-wise) and look for voxels that
# share a face. Those are considered a cluster. Save the size of each cluster
# Get the null cluster-size distribution. Do the same with the 10 000 chance
# accuracy distributions.
# FOURTH STEP
# Calculate voxel wize p-values for significant clusters. That is the sum of 
# all accuracies larger than the observed in the chance distribution.
# Save the result. 

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
import sys
import argparse
from pathlib import Path
import git
import glob
# data structuration and calculations
import numpy as np   # most important numerical calculations
# library for neuroimaging
import nibabel as nib
from nilearn.image import smooth_img, new_img_like
# optimize time performance
import time
from multiprocessing import Pool
import matplotlib.pyplot as plt

# A function to check if a given cell
# (row, col) can be included in DFS
# shamelessly stolen from 
# https://www.geeksforgeeks.org/find-length-largest-region-boolean-matrix/
# and modified for 3D matrix
def isSafe(M, row, col, ais, visited):
    '''isSave: Checks if a voxel (one cell in a 3D Matrix) has already been
    checked.
    Input: 
        M: 3D matrix. This one is investigated in the whole process
        row: index of row
        col: index of column
        ais: index of aisle (3rd dimension)
        visited: if the cell was already checked
    Returns:
        If the cell entry is True AND has not been visited before.'''
    ROW = M.shape[0]
    COL = M.shape[1]
    AIS = M.shape[2]
    
    # row number is in range, column number is in
    # range and value is 1 and not yet visited
    return ((row >= 0) and (row < ROW) and
            (col >= 0) and (col < COL) and
            (ais >= 0) and (ais < AIS) and
            (M[row,col,ais] and not visited[row,col,ais]))
 
# A utility function to do DFS for a 3D
# boolean matrix. It only considers
# the 6 neighbours as adjacent faces
def DFS(M, row, col, ais, visited, count, indices):
 
    # These arrays are used to get row and column and aisle
    # numbers of 6 neighbouring faces of a given cell
    rowNbr = [-1,1,0 ,0, 0,0]
    colNbr = [0 ,0,-1,1, 0,0]
    aisNbr = [0 ,0,0 ,0,-1,1]
 
    # Mark this cell as visited
    visited[row,col,ais] = True
    # append indices of visited cell
    indices.append([row,col,ais])
 
    # Recur for all connected neighbours
    for k in range(6):
        if (isSafe(M, row + rowNbr[k],
                   col + colNbr[k], 
                   ais + aisNbr[k], visited)):
 
            # increment region length by one
            count[0] += 1
            DFS(M, row + rowNbr[k],
                col + colNbr[k], ais + aisNbr[k],
                visited, count,indices)
    
 
# The main function that returns largest
# length region of a given boolean 3D matrix
def largestRegion(M):
    ROW = M.shape[0]
    COL = M.shape[1]
    AIS = M.shape[2]
 
    # Make a bool array to mark visited cells.
    # Initially all cells are unvisited
    visited = np.zeros(M.shape)
 
    # Initialize result as 0 and travesle
    # through the all cells of given matrix
    cluster_sizes = []
    all_cluster_indices = []
    for i in range(ROW):
        for j in range(COL):
            for k in range(AIS):
 
                # If a cell with value 1 is not
                if (M[i,j,k] and not visited[i,j,k]):
     
                    # visited yet, then new region found
                    count = [1]
                    clust_indices = []
                    DFS(M, i, j, k, visited, count,clust_indices)
                    all_cluster_indices.append(clust_indices)
     
                    # maximum region
                    cluster_sizes.append(count[0])
    return all_cluster_indices, cluster_sizes

# get start time
T_START = time.time()

################################
# Command line input arguments #
################################

# create a parser object which handles the input variables
parser = argparse.ArgumentParser()

# add all the input arguments
parser.add_argument("--what", "-w", nargs="?",const='effect', default='effect', 
                    type=str)
parser.add_argument("--data", "-d", nargs="?", const='pre', default='pre', 
                    type=str)
parser.add_argument("--over",  "-o",   nargs='?',  const='objects', 
                    default='objects')
parser.add_argument("--analyzed",           nargs='?', const='moment',  
                    default='moment',   type=str)

# parse the arguments to a parse-list(???)
ARGS = parser.parse_args()
# assign values 
WHAT = ARGS.what
DATA = ARGS.data
OVER = ARGS.over
ANALYZED    = ARGS.analyzed

# variables for path selection and data access
HOME            = str(Path.home())
PROJ_DIR        = os.path.join(HOME, 'Documents/Magic_fMRI/DATA/MRI')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')

RAWDATA_DIR     = os.path.join(PROJ_DIR, 'rawdata')
ANALYSIS        = 'ROI-analysis'
if WHAT == 'effect':
    DATA_TO_USE = 'decode_effect'
    NUM_LABELS  = 3
elif WHAT == 'pre-post':
    DATA_TO_USE = 'decode_pre_vs_post'
    NUM_LABELS  = 2
elif WHAT == 'mag-nomag':
    DATA_TO_USE = 'decode_magic_vs_nomagic'
    NUM_LABELS  = 2
else:
    raise

if ANALYZED == 'moment':
    data_analyzed = 'SpecialMoment'
elif ANALYZED == 'video':
    data_analyzed = 'WholeVideo'
else:
    raise
    
DATA_DIR        = os.path.join(DERIVATIVES_DIR, 'decoding', 'decoding_magic', 
                               DATA_TO_USE, DATA+'_videos','over_'+ OVER, 
                               data_analyzed, 'SearchLight','LDA')

RESULTS_DIR     = os.path.join(DATA_DIR,'group-statistics')
BOOTSTRAPPS     = glob.glob(os.path.join(RESULTS_DIR,'s4bootstrapped_*'))

if not os.path.isdir(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

# SECOND STEP    
CRIT_P_VAL          = 0.001
PERCENTILE          = 100 - 100*CRIT_P_VAL
CLUSTER_SIZE_ALPHA  = 0.05  # CHECK FOR THE VALUE IN THE PAPER
SMOOTH_KERNEL_SIZE  = None

# get the real mean decoding accuracy 
mean_accuracy_path = os.path.join(RESULTS_DIR,'s4mean_accuracy.nii')
mean_img           = smooth_img(mean_accuracy_path,fwhm=SMOOTH_KERNEL_SIZE)
mean_accuracy_map  = mean_img.get_fdata()       # get data from NIfTI image
mean_img.uncache() 

# set unreasonable values (below chance to zero)

# get dimensions of mean_accuracy_map to know the dimensions of all used chance
# distribution maps
BRAIN_DIMENSIONS = mean_accuracy_map.shape

# Using the inference method described in Stelzer (2013)
# From each subject randomly select one permutation image and create a mean
# over subjects. Do this several times (10^5 was used in Stelzer) and get a 
# distribution

# Since it would take too much memory to read in all chance distribution maps
# at once to calculate voxel-wise threshold values, we do it slice by slice

# empty list to store permuted slices into
boot_classification_images  = smooth_img(BOOTSTRAPPS,fwhm=SMOOTH_KERNEL_SIZE)
boot_classification_maps    = []
for d,draw in enumerate(boot_classification_images):
    if d%500 == 0:
        print('Reading in img Nr' + str(d))
    img_data        = draw.get_fdata()   # get data from NIfTI image
    boot_classification_maps.append(img_data)
    del img_data

# A map storing the individual accuracy values for each voxel that is 
# considered significant (the percentile is used)
crit_acc_value_map = np.percentile(boot_classification_maps,
                                   PERCENTILE,
                                   axis=0)

# save map consisting of critical values
results = new_img_like(ref_niimg=mean_img,
                       data=crit_acc_value_map)
nib.save(results,os.path.join(RESULTS_DIR, 'crit_acc_value_map.nii'))

# THIRD STEP
# Get cluster_sizes and locations in the mean accuracy maps
# EXTRA STEP: set maximum recursion depth to number of voxels exceeding 
# significant threshold
sig_voxel_map = mean_accuracy_map>crit_acc_value_map
num_sig_voxel = sig_voxel_map.sum(axis=None)
sys.setrecursionlimit(num_sig_voxel)
cluster_indices, decoding_cluster_sizes = largestRegion(sig_voxel_map)

#############################################
# get 'null distribution' of cluster sizes: #
#############################################
# Do so by iterate over all chance distributions and read them in one after
# another
# A bool 4D array giving significant voxels for a specific bootstrap
perm_significant_maps = boot_classification_maps>crit_acc_value_map

# empty list to store cluster sizes in
perm_cluster_sizes = []
# Use multiprocessing.Pool to make it faster
with Pool(25) as pool:
    for cl_size in pool.map(largestRegion,perm_significant_maps):
        perm_cluster_sizes.extend(cl_size[1])
        
del perm_significant_maps
        
# normalize histogram of cluster sizes: 
# (occurence/total number of detected clusters)
cluster_sizes, cluster_counts   = np.unique(perm_cluster_sizes,
                                            return_counts=True)
norm_cluster_hist               = cluster_counts/len(perm_cluster_sizes)
# flip normalized cluster histogram to calculate p-values for cluster sizes
norm_cluster_hist = np.flip(norm_cluster_hist)
# calculate p-values for cluster sizes: p_cluster = sum_{s'>s} H_cluster(s')
# where H_cluster is the normalized histogram 
p_cluster = [sum(norm_cluster_hist[0:i+1]) 
             for i in range(len(norm_cluster_hist))]
p_cluster = np.array(p_cluster)

sig_cluster_sizes_uncor         = cluster_sizes[p_cluster<CLUSTER_SIZE_ALPHA]
cluster_size_threshold_uncor    = min(sig_cluster_sizes_uncor)

# Get the FDR corrected cluster size. Check if the i-th p-value is smaller than
# alpha/(number of tests - i). A soon as it gets larger than that you stop.
for p_idx, p_clu in enumerate(p_cluster):
    test = p_clu<CLUSTER_SIZE_ALPHA/(len(p_cluster)-p_idx)
    if not test:
        break
    
sig_cluster_sizes_cor         = cluster_sizes[p_idx-1:]
cluster_size_threshold_cor    = min(sig_cluster_sizes_cor)

# FOURTH STEP
# create a zero matrix with MNI dimensions to create a map that stores p-values
# of voxels within significant clusters
p_voxel_map = np.zeros(shape=mean_accuracy_map.shape)

# iterate over the clusters in actual decoding map
for c,clust in enumerate(cluster_indices):
    # if the cluster size is larger than the cluster theshold calculate the
    # p-value for each voxel in the cluster
    if decoding_cluster_sizes[c]>cluster_size_threshold_cor:
        # go through the voxels in your significant cluster
        for voxel in clust:
            # calculate voxel-wise p-values for these clusters: 
            # p_voxel = sum_{a'>a} H_voxel(a')
            # H_voxel is normalized chance distribution for this voxel
            # a original accuracy
            # First: calculate normalized chance distribution for the specific
            # voxel in question
            chance_accuracies = [item[voxel[0],voxel[1],voxel[2]]
                                 for item in boot_classification_maps]
            normed_chance_dist = chance_accuracies/sum(chance_accuracies)
            
            # read out the 'real' decoding accuracy
            real_accuracy = mean_accuracy_map[voxel[0],voxel[1],voxel[2]]
            
            # calculate p-value and store it in the p-value map
            tmp_p_val = sum(normed_chance_dist[chance_accuracies>real_accuracy])
            p_voxel_map[voxel[0],voxel[1],voxel[2]] = tmp_p_val

# save map consisting of critical values
sig_voxel_map_corr = p_voxel_map<0.05

results = new_img_like(ref_niimg=mean_img, data=p_voxel_map)
nib.save(results,os.path.join(RESULTS_DIR, 'voxel_p-value_map.nii'))

results = new_img_like(ref_niimg=mean_img, data=sig_voxel_map_corr)
nib.save(results,os.path.join(RESULTS_DIR, 'significant_voxel_map.nii'))

# from generated normalized null distribution create a histogram
fig = plt.figure()
plt.hist(norm_cluster_hist,bins=50)
fig.savefig(os.path.join(RESULTS_DIR,'clustersize_null_distribution.png'))

#--------------------REMOVE ALL BOOTSTRAPPED NIFTI IMAGES--------------------#
for file in BOOTSTRAPPS:
    os.remove(file)

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
with open(os.path.join(RESULTS_DIR,'stelzer_analysis-logfile.txt'), 'w+') as writer:
    writer.write('Codeversion: {} \n'.format(git_hash))