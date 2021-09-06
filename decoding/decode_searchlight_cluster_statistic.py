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
import argparse
from pathlib import Path
import git
import glob
# data structuration and calculations
import numpy as np   # most important numerical calculations
# library for neuroimaging
import nibabel as nib
from nilearn.image import new_img_like
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
    result = []
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
                    result.append(count[0])
    return result, all_cluster_indices

# get start time
T_START = time.time()

################################
# Command line input arguments #
################################

# create a parser object which handles the input variables
parser = argparse.ArgumentParser()

# add all the input arguments
parser.add_argument("--data",       "-d",   nargs="?",  const='pre',    
                    default='pre',  type=str)

# parse the arguments to a parse-list(???)
ARGS = parser.parse_args()
# assign values 
DATA     = ARGS.data


# variables for path selection and data access
HOME            = str(Path.home())
PROJ_DIR        = os.path.join(HOME, 'Documents/Master_Thesis/DATA/MRI')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')

RAWDATA_DIR     = os.path.join(PROJ_DIR, 'rawdata')
ANALYSIS        = 'ROI-analysis'
if DATA == 'pre':
    DATA_TO_USE = 'decode_effect_on_premagic'
elif DATA == 'post':
    DATA_TO_USE = 'decode_effect_on_postmagic'
elif DATA == 'all':
    DATA_TO_USE = 'decode_effect_on_allmagic'
elif DATA == 'pre-post':
    DATA_TO_USE = 'decode_pre_vs_post'
elif DATA == 'mag-nomag':
    DATA_TO_USE = 'magic_vs_nomagic'
else:
    raise
DATA_DIR        = os.path.join(DERIVATIVES_DIR, 'decoding', 'decoding_magic', 
                               DATA_TO_USE, 'SpecialMoment','SearchLight','LDA')
RESULTS_DIR     = os.path.join(DATA_DIR,'group-statistics')
BOOTSTRAPPS     = glob.glob(os.path.join(RESULTS_DIR,'perm_*'))

if not os.path.isdir(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

# SECOND STEP    
CRIT_P_VAL          = 0.001 # FOR NOW! change after larger bootstrap method
CLUSTER_SIZE_ALPHA  = 0.05  # CHECK FOR THE VALUE IN THE PAPER

# p=(C+1)/(N+1)
# p: p-value
# C: Number of accuracies better than the 'original' accuray
# N: Number of permutations
# ==> C = p*(N+1)-1
# C is the index of the bootstrap array which's value depicts the cirtical 
# value that is considered significant
# With a large N we can simplify it to C = p*N

crit_acc_idx = CRIT_P_VAL*len(BOOTSTRAPPS)
crit_acc_idx = round(crit_acc_idx)

# get the real mean decoding accuracy 
mean_accuracy_path = os.path.join(RESULTS_DIR,'mean_accuracy.nii')
img                 = nib.load(mean_accuracy_path)
mean_accuracy_map   = img.get_fdata()                  # get data from NIfTI image
img.uncache() 

# get dimensions of mean_accuracy_map to know the dimensions of all used chance
# distribution maps
BRAIN_DIMENSIONS = mean_accuracy_map.shape

# Using the inference method described in Stelzer (2013)
# From each subject randomly select one permutation image and create a mean
# over subjects. Do this several times (10^5 was used in Stelzer) and get a 
# distribution

# Since it would take too much memory to read in all chance distribution maps
# at once to calculate voxel-wise threshold values, we do it slice by slice

# empty list where the voxel-wise threshold value slices are stored in
crit_acc_value_map = []
for sl in range(BRAIN_DIMENSIONS[0]):

    # empty list to store permuted slices into
    perm_classification_slices = []
    for draw in BOOTSTRAPPS:
        img             = nib.load(draw)
        img_data        = img.get_fdata()   # get data from NIfTI image
        img.uncache() 
        perm_classification_slices.append(img_data[sl])
        
    perm_classification_slices  = np.array(perm_classification_slices)
    sorted_perm_class_images    = perm_classification_slices.copy()
    sorted_perm_class_images.sort(axis=0)
    
    crit_acc_value_map.append(sorted_perm_class_images[-crit_acc_idx])
    
crit_acc_value_map = np.array(crit_acc_value_map)

# save map consisting of critical values
results = new_img_like(ref_niimg=img,data=crit_acc_value_map)
nib.save(results,os.path.join(RESULTS_DIR, 'crit_acc_value_map.nii'))

# delete the copy of perm_classification_images, once it is not needed anymore
# it takes up a bunch of RAM!
del sorted_perm_class_images

# THIRD STEP
# Get cluster_sizes and locations in the mean accuracy maps
decoding_cluster_sizes, cluster_indices = largestRegion(mean_accuracy_map>crit_acc_value_map)

# with Pool(25) as pool:
#     for cl_size,_ in pool.map(largestRegion,tmp_bool):
#         perm_cluster_sizes.extend(cl_size)

#############################################
# get 'null distribution' of cluster sizes: #
#############################################
# Do so by iterate over all chance distributions and read them in one after
# another

# empty list to store cluster sizes in
perm_cluster_sizes = []
# Use multiprocessing.Pool to make it faster
with Pool(25) as pool:
     for img in pool.map(nib.load,BOOTSTRAPPS):
         perm_cl = img.get_fdata()
         tmp_bool = perm_cl>crit_acc_value_map
         cl_size,_ = largestRegion(tmp_bool)
         perm_cluster_sizes.extend(cl_size)
        
img.uncache()
del tmp_bool
# normalize histogram of cluster sizes: 
# (occurence/total number of detected clusters)
cluster_sizes, cluster_counts   = np.unique(perm_cluster_sizes,return_counts=True)
cluster_sizes                   = np.flip(cluster_sizes)
norm_cluster_hist               = cluster_counts/len(perm_cluster_sizes)
# flip normalized cluster histogram to calculate p-values for cluster sizes
norm_cluster_hist = np.flip(norm_cluster_hist)
# calculate p-values for cluster sizes: p_cluster = sum_{s'>s} H_cluster(s')
# where H_cluster is the normalized histogram 
p_cluster = [sum(norm_cluster_hist[0:i]) 
             for i in range(len(norm_cluster_hist))]
p_cluster = np.array(p_cluster)

sig_cluster_sizes       = cluster_sizes[p_cluster<CLUSTER_SIZE_ALPHA]
cluster_size_threshold  = min(sig_cluster_sizes)

# FOURTH STEP
# create a zero matrix with MNI dimensions to create a map that stores p-values
# of voxels within significant clusters
p_voxel_map = np.zeros(shape=mean_accuracy_map.shape)

# iterate over the clusters in actual decoding map
for c,clust in enumerate(cluster_indices):
    # if the cluster size is larger than the cluster theshold calculate the
    # p-value for each voxel in the cluster
    if decoding_cluster_sizes[c]>cluster_size_threshold:
        # go through the voxels in your significant cluster
        for voxel in clust:
            # calculate voxel-wise p-values for these clusters: 
            # p_voxel = sum_{a'>a} H_voxel(a')
            # H_voxel is normalized chance distribution for this voxel
            # a original accuracy
            # First: calculate normalized chance distribution for the specific
            # voxel in question
            perm_classification_voxel = []
            for draw in BOOTSTRAPPS:
                img             = nib.load(draw)
                img_data        = img.get_fdata()   # get data from NIfTI image
                img.uncache() 
                perm_classification_voxel.append(img_data[voxel])
            chance_accuracies = perm_classification_voxel
            normed_chance_dist = chance_accuracies/sum(chance_accuracies)
            
            # read out the 'real' decoding accuracy
            real_accuracy = mean_accuracy_map[voxel[0],voxel[1],voxel[2]]
            
            # calculate p-value and store it in the p-value map
            tmp_p_val = sum(normed_chance_dist[chance_accuracies>real_accuracy])
            p_voxel_map[voxel[0],voxel[1],voxel[2]] = tmp_p_val

# save map consisting of critical values
results = new_img_like(ref_niimg=img,data=p_voxel_map)
nib.save(results,os.path.join(RESULTS_DIR, 'voxel_p-value_map.nii'))

# from generated normalized null distribution create a histogram
fig = plt.figure()
plt.hist(norm_cluster_hist,bins=50)
fig.savefig(os.path.join(RESULTS_DIR,'clustersize_null_distribution.png'))


print('Done')