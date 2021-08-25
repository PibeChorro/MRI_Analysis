#!/gpfs01/bartels/user/vplikat/anaconda3/bin/python

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
import matplotlib.pyplot as plt

# A function to check if a given cell
# (row, col) can be included in DFS
# shamelessly stolen from 
# https://www.geeksforgeeks.org/find-length-largest-region-boolean-matrix/
# and modified for 3D matrix

 
def isSafe(M, row, col, ais, visited):
    ROW = M.shape[0]
    COL = M.shape[1]
    AIS = M.shape[2]
    
    # row number is in range, column number is in
    # range and value is 1 and not yet visited
    return ((row >= 0) and (row < ROW) and
            (col >= 0) and (col < COL) and
            (ais >= 0) and (ais < AIS) and
            (M[row,col,ais] and not visited[row,col,ais]))
 
# A utility function to do DFS for a 2D
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
parser.add_argument("--runs", "-r", nargs="?", const='pre', default='pre', type=str)

# parse the arguments to a parse-list(???)
ARGS = parser.parse_args()
# assign values 
RUNS_TO_USE     = ARGS.runs


# variables for path selection and data access
HOME            = str(Path.home())
PROJ_DIR        = os.path.join(HOME, 'Documents/Master_Thesis/DATA/MRI')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')

RAWDATA_DIR     = os.path.join(PROJ_DIR, 'rawdata')
ANALYSIS        = 'ROI-analysis'
DATA_DIR        = os.path.join(DERIVATIVES_DIR, 'decoding', 'decoding_magic', 
                               'decode_effect_on_'+RUNS_TO_USE+'magic',
                               'SpecialMoment','SearchLight_old','LDA')
RESULTS_DIR     = os.path.join(DATA_DIR,'group-statistics')
BOOTSTRAPPS     = glob.glob(os.path.join(RESULTS_DIR,'perm_*'))

if not os.path.isdir(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)
    
CRIT_P_VAL          = 0.05    # FOR NOW! change after larger bootstrap method
CLUSTER_SIZE_ALPHA  = 0.05

# get the real mean decoding accuracy 
mean_accuracy_path = os.path.join(RESULTS_DIR,'mean_accuracy.nii')
img                 = nib.load(mean_accuracy_path)
mean_accuracy_map   = img.get_fdata()                  # get data from NIfTI image
img.uncache() 

# Using the inference method described in Stelzer (2013)
# From each subject randomly select one permutation image and create a mean
# over subjects. Do this several times (10^5 was used in Stelzer) and get a 
# distribution
# empty list to store permuted maps into
perm_classification_images = []
for draw in BOOTSTRAPPS:
    img             = nib.load(draw)
    img_data        = img.get_fdata()   # get data from NIfTI image
    img.uncache() 
    perm_classification_images.append(img_data)
    
perm_classification_images  = np.array(perm_classification_images)
sorted_perm_class_images    = perm_classification_images.copy()
sorted_perm_class_images.sort(axis=0)

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

crit_acc_value_map = sorted_perm_class_images[-crit_acc_idx]

# save map consisting of critical values
results = new_img_like(ref_niimg=img,data=crit_acc_value_map)
nib.save(results,os.path.join(RESULTS_DIR, 'crit_acc_value_map.nii'))

# delete the copy of perm_classification_images, once it is not needed anymore
# it takes up a bunch of RAM!
del sorted_perm_class_images

# Get cluster_sizes and locations in the mean accuracy maps
decoding_cluster_sizes, cluster_indices = largestRegion(mean_accuracy_map>crit_acc_value_map)

# get 'null distribution' of cluster sizes:
# iterate over bootstrapped accuracy maps and get all existing clusters in the
# map. Add them all up to one huge list and you have the distribution
perm_cluster_sizes = []
for perm in perm_classification_images:
    tmp_cluster_size, _ = largestRegion(perm>crit_acc_value_map)
    perm_cluster_sizes.extend(tmp_cluster_size)
    
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
sig_cluster_sizes = cluster_sizes[p_cluster<CLUSTER_SIZE_ALPHA]
cluster_size_threshold = min(sig_cluster_sizes)

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
            chance_accuracies = perm_classification_images[:,voxel[0],voxel[1],voxel[2]]
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