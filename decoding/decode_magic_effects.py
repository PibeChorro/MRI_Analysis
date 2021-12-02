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
# The purpose of this script is to use a machine learning algorithm to predict
# the magical effects performed on one object, based on training data from the
# other two objects. To be able to test for statistical significance, 
# permutation testing is applied.
# All this is done for a set of ROIs, but on a single subject level

################################
# FUNCTIONALITY OF THIS SCRIPT #
################################
# FIRST STEP 
# Get all important path information and names of ROIs, magic EFFECTS (Appear, 
# Change, Vanish) and information about the analysis from command line input.
# SECOND STEP 
# Read in the SPM.mat file created by SPM12 when the GLM is estimated.
# From this SPM.mat file we read out the names of the beta NIfTI images and
# the names of the regressors, that correspond to the beta image.
# From the Regressor names the run number (1-12) is extracted and added to the
# DataFrame. Then all regressors of NO interest (realignment, controll and 
# surprise videos, etc.) are removed, as well as the data we do not want to
# test (when analyzing pre revelation data, the post revelation data is removed
# and vise versa).
# Again from the Regressor name the label (=the effect) is extracted and added
# to the dataframe. Finally the chunks are defined. Chunks are needed for cross
# validation (meaning the whole dataset is sperated in chunks and each chunk is
# used for testing once).
# THIRD STEP
# Read in all beta images that are left and store them in a huge 2D array 
# (every beta image is flattened turning the 3D image in a 1D array, these
# arrays are then formed into a 2D array)
# Iterate over all ROIs, read in the ROI mask image and apply it on the 2D 
# array. The remaining data can be manipulated according to given flags 
# (scaling, cutoff etc.) and is then fed into the permutation_test_score method
# to calculate the measured accuracy and get a null distribution from the 
# permutation. 
# FOURTH STEP
# Save the result in a hdf5 file and make a barplot with the accuracies for
# each ROI

######################
# COMMAND LINE FLAGS #
######################
# --sub: the subject ID that shall be analyzed
# --smooth: if the script should read in beta images that are the result of a
# GLM based on smoothed functional images
# --algorythm: which algorythm should be used. Currently implemented SVM and LDA
# --scaling: if the beta image data should be scaled. Currently implemented z,
# min0max1 and de-meaning
# --cutoff: If values should be cut off. IMPORTANT makes only sense if scaling
# was applied
# --feature: Feature transformation of the data for dimension reduction. 
# Currently implemented PCA
# --kernels: How many kernels should be used to parallize the permutation testing
# --runs: Which data should be used. Either pre, post revelation or all data
# together
# --perms: How many permutations are applied

#############
# LIBRARIES #
#############

# interact with the operating system 
import os
import argparse
from pathlib import Path
import git
import glob
# import/export data
import h5py
# data structuration and calculations
import pandas as pd  # to create data frames
import numpy as np   # most important numerical calculations
# read in mat files
import readmat
# needed to extract the run number out of the parentesis of the string in the SPM.mat file
import re
# library for neuroimaging
import nibabel as nib
# machine learning algorithms and stuff
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.svm import SVC
from sklearn.model_selection import PredefinedSplit, permutation_test_score
# optimize time performance
import time
from tqdm import tqdm
# plotting
import matplotlib.pyplot as plt

# get start time
T_START = time.time()

# FIRST STEP 
################################
# Command line input arguments #
################################

# create a parser object which handles the input variables
parser = argparse.ArgumentParser()

# add all the input arguments
parser.add_argument("--sub",        "-s",                               
                    default='sub-01')         # subject
parser.add_argument("--over",       "-o",   nargs='?',  const='objects',    
                    default='tricks')
parser.add_argument("--smooth",             nargs='?',  const=0,        
                    default=0,      type=int) # what data should be used
parser.add_argument("--algorythm",    "-a",   nargs='?',  const='LDA',    
                    default='LDA')
parser.add_argument("--scaling",            nargs='?',  const='None',   
                    default='None', type=str)
parser.add_argument("--cutoff",     "-c",   nargs='?',  const=np.inf,   
                    default=np.inf, type=float) # if and with which value (in std) data is cut off 
parser.add_argument("--feature",    "-f",   nargs='?',  const='None',   
                    default='None', type=str)
parser.add_argument("--kernels",    "-k",   nargs='?',  const=12,       
                    default=12,     type=int)   # how many processes should be run in parallel
parser.add_argument("--runs",       "-r",   nargs="?",  const='pre',    
                    default='pre',  type=str)
parser.add_argument("--perms",      "-p",   nargs="?",  const=1000,     
                    default=1000,   type=int)   # how many permutations
# parse the arguments to a parse-list(???)
ARGS = parser.parse_args()
# assign values 
SUB             = ARGS.sub
OVER            = ARGS.over
N_PROC          = ARGS.kernels
SMOOTHING_SIZE  = ARGS.smooth
DECODER         = ARGS.algorythm
CUTOFF          = ARGS.cutoff
FEAT_TRANS      = ARGS.feature
SCALE           = ARGS.scaling
RUNS_TO_USE     = ARGS.runs
N_PERMS         = ARGS.perms

# based on command line flag decide what decoding algorithm should be used
if DECODER =='LDA':
    my_decoder          = LDA(solver='lsqr', shrinkage='auto')
    decoder_parameters  = os.path.join(
        'scale_'+SCALE,
        'cutoff_'+str(CUTOFF),
        'feat_trans_'+FEAT_TRANS)
elif DECODER == 'SVM':
    SVM_C = 1
    decoder_parameters  = os.path.join(
        'scale_'+SCALE,
        'cutoff_'+str(CUTOFF),
        'feat_trans_'+FEAT_TRANS,
        'C_'+str(SVM_C))
    my_decoder          = SVC(kernel='linear', C=SVM_C)
else:
    raise
        
# based on command line flag decide what data should be used. 
# 'pre' = pre revelation data
# 'post' = post revelation data
# 'all' = pre and post data together
if RUNS_TO_USE == 'pre':
    runs_of_interest = [1,2,5,6,9,10]
elif RUNS_TO_USE == 'post':
    runs_of_interest = [3,4,7,8,11,12] 
elif RUNS_TO_USE == 'all':
    runs_of_interest = [1,2,3,4,5,6,7,8,9,10,11,12] 
else:
    raise

################################################
# VARIABLES FOR PATH SELECTION AND DATA ACCESS #
################################################
HOME            = str(Path.home())
PROJ_DIR        = os.path.join(HOME, 'Documents/Magic_fMRI/DATA/MRI')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')
# where to look for the beta images
if SMOOTHING_SIZE > 0:
    GLM_DATA_DIR    = str(SMOOTHING_SIZE)+'mm-smoothed-nativespace' 
    FLA_DIR         = os.path.join(DERIVATIVES_DIR,'spm12',
                               'spm12-fla','WholeBrain',
                               'EveryVideo',GLM_DATA_DIR,
                               'SpecialMoment',SUB)
else:
    GLM_DATA_DIR    = 'nativespace' 
    FLA_DIR         = os.path.join(DERIVATIVES_DIR,'spm12',
                               'spm12-fla','WholeBrain',
                               'EveryVideo',GLM_DATA_DIR,
                               'SpecialMoment',SUB)
# where the ROI masks can be found
FREESURFER_DIR  = os.path.join(DERIVATIVES_DIR, 'freesurfer')
ROI_DIR         = os.path.join(FREESURFER_DIR,SUB,'corrected_ROIs')
# where the .mat file can be found created by SPM12
SPM_MAT_DIR     = os.path.join(FLA_DIR, 'SPM.mat')
# wheret to store the results
ANALYSIS        = 'ROI-analysis'
RESULTS_DIR     = os.path.join(DERIVATIVES_DIR, 'decoding', 'decoding_magic', 
                               'decode_effect',''+RUNS_TO_USE+'_videos', 
                               'over_' + OVER, 'SpecialMoment', ANALYSIS, 
                               SUB)
if not os.path.isdir(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

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

# What should be predicted
LABEL_NAMES = [
    'Appear',
    'Change',
    'Vanish'
]

# empty lists that will be filled with the results to plot after calculation
decode_accuracy = []
decode_p_value = []

# create a 'random' seed number for the permutation based on the subject name
rng_seed = 0
for letter in SUB:
    rng_seed += ord(letter)
    
print ('Analysing subject: {}'.format(SUB))
print ('Getting ROIs from:	 {}'.format(FREESURFER_DIR))
print ('Saving data at:	 {}'.format(RESULTS_DIR))
print ('Number of permutations:	 {}'.format(N_PERMS))

# SECOND STEP 
########################################
# reading in the necessary information #
########################################

# From the previously created SPM.mat file we read in the information we need
# The filenames of our beta images
SPM_BETADICT = readmat.load(SPM_MAT_DIR, isStruct=True)['SPM']['Vbeta']
BETA_DIRS = [f['fname'] for f in SPM_BETADICT]
# The corresponding Regressor names - are in form of 'Sn(<run-number>) <Regressor-Name>*bf(1)'
SPM_REGRESSORS  = readmat.load(SPM_MAT_DIR,isStruct=True)['SPM']['xX']['name']

# store beta filenames and regressornames in a dictionary
data_dict = {
    'Regressors': SPM_REGRESSORS,
    'BetaNames': BETA_DIRS
}

# convert dictionary into a pandas DataFrame for further analysis
label_df = pd.DataFrame(data_dict, columns=data_dict.keys())

# a complex process to throw out regressors of no interest (like realignment)
regressors_of_interest  = [True if 'Magic' in n else False for n in SPM_REGRESSORS]
# throw out all rows of regressors of no interest
label_df = label_df.iloc[regressors_of_interest]

# This complex loop is necessary to get the run number out of the regressor name
run_info    = [' '.join(re.findall(r"\((\d+)\)",string)) for string in label_df.Regressors]
runs        = [int(s_filter.split()[0]) for s_filter in run_info]

# add further data to DataFrame
label_df['Runs'] = runs # In which run
# only keep data from the runs you are interested in
label_df = label_df[label_df.Runs.isin(runs_of_interest)]

VideoNames  = [''.join(re.search(" (.+?)\*",string).group(1)) 
            for string in label_df.Regressors]
label_df['VideoNames'] = VideoNames

# depending on whether we want to decode over objects or trick versions the 
# chunks change
if OVER == 'objects':
    label_df['Chunks']  = (label_df.Runs-1)//4
    chunks              = np.asarray(label_df.Chunks)
    ps = PredefinedSplit(chunks)
elif OVER == 'tricks':
    # training on all tricks of verion 1 and test on all tricks of version 2 
    # caused a shifted null distribution below chance, hence resulting in a 
    # below chance max-statistic null distribution for multiple comparison. 
    # To resolve the problem we train on tricks of version 1 in odd trials 
    # and test on tricks of version 2 in even trials. 
    # DISADVANTAGE: less data
    # ADVANTAGE: four, instead of two validation folds (version 1/2 x odd/even)
    
    # stange and not beautiful solution (that I don't really understand) taken
    # from here:
    # (https://stackoverflow.com/questions/28837633/pandas-get-position-of-a-given-index-in-dataframe)
    
    label_df['Version'] = label_df.VideoNames.str.contains('1')
    train1  = label_df.index.get_indexer_for((label_df[(label_df.Runs%2==0) & (label_df.Version==0)].index))
    test1   = label_df.index.get_indexer_for((label_df[(label_df.Runs%2==1) & (label_df.Version==1)].index))
    train2  = label_df.index.get_indexer_for((label_df[(label_df.Runs%2==1) & (label_df.Version==0)].index))
    test2   = label_df.index.get_indexer_for((label_df[(label_df.Runs%2==0) & (label_df.Version==1)].index))
    train3  = label_df.index.get_indexer_for((label_df[(label_df.Runs%2==0) & (label_df.Version==1)].index))
    test3   = label_df.index.get_indexer_for((label_df[(label_df.Runs%2==1) & (label_df.Version==0)].index))
    train4  = label_df.index.get_indexer_for((label_df[(label_df.Runs%2==1) & (label_df.Version==1)].index))
    test4   = label_df.index.get_indexer_for((label_df[(label_df.Runs%2==0) & (label_df.Version==0)].index))
    ps = [[train1, test1], [train2,test2],
          [train3, test3], [train4,test4]]
else:
    raise

label_df['Labels']  = np.nan                # Labels
# Check for every entry in Regressors if it contains one of the label names. If so, assign the label name
for l in LABEL_NAMES:
    label_df.Labels = np.where(label_df.Regressors.str.contains(l),l,label_df.Labels)

# THIRD STEP
# read in all beta from regressors of interest (flatten and then combine in on
# large 2D array)
betas = []                                              # empty list to store data arrays in
for b, beta in enumerate(label_df.BetaNames):
    beta_nii    = nib.load(os.path.join(FLA_DIR,beta))  # read in beta NIfTI image
    beta_data   = beta_nii.get_fdata()                  # get data from NIfTI image
    beta_nii.uncache()                                  # remove beta image from memory
    beta_data   = beta_data.flatten()                   # convert into one-dimensional array
    betas.append(beta_data)                             # append array on betas list
betas = np.array(betas)

# inner loop - iterating over mask (=ROIs)
for r, roi in tqdm(enumerate(ROIS)):
    output_dir = os.path.join(RESULTS_DIR,roi + '.hdf5')   # where to store the results

    # Get all NifTi files containing the name of your ROI. Read them in and combine them to one ROI
    maskdir_list = glob.glob(os.path.join(ROI_DIR,'*' + roi + '*.nii'))
    masklist = []
    for mask in maskdir_list:
        mask_nii = nib.load(mask)
        mask_img = mask_nii.get_fdata()
        mask_nii.uncache()
        mask_img = np.asarray(mask_img)
        mask_img = mask_img.flatten()
        masklist.append(mask_img)
        
    # sum up the ROI mask list and thus create one ROI
    ROI = np.sum(masklist,axis=0)
    # turn into boolean values
    ROI = ROI>0
    
    # get data you want to use in the ROI analysis
    ROI_data = []
    ROI_data.append([beta[ROI & ~np.isnan(beta)] for beta in betas])
    # the list comprehension wraps the matrix in an additional, unnecessary array
    ROI_data = np.array(ROI_data [0])               
    
    # apply scaling and cutoff according to command line input
    if SCALE == 'z':
        # z-transform data within features
        ROI_data = ROI_data - ROI_data.mean(axis=0)
        ROI_data = ROI_data / ROI_data.std(axis=0)    
        # set outliers to CUTOFF value (only makes sense after z-scaling)
        ROI_data[ROI_data>CUTOFF]   = CUTOFF
        ROI_data[ROI_data<-CUTOFF]  = -CUTOFF
    elif SCALE == 'min0max1':
        # substract the minimum value (set min=0)
        # then divide by the maximum value (set max=1)
        ROI_data = ROI_data - ROI_data.min(axis=0)
        ROI_data = ROI_data / ROI_data.max(axis=0)
    elif SCALE == 'mean':
        # substract the mean (set mean=0)
        ROI_data = ROI_data - ROI_data.mean(axis=0)
    
    # apply feature transformation and dimension reduction according to 
    # command line input
    if FEAT_TRANS == 'PCA':
        n_components = min (ROI_data.shape)
        my_PCA = PCA(n_components=n_components)
        ROI_data = my_PCA.fit_transform(ROI_data)

    
    # the actual decoding
    targets                 = np.asarray(label_df.Labels)
    runs_for_permutation    = np.asarray(label_df.Runs)
    if N_PERMS > 0:
        res = permutation_test_score(
            estimator=my_decoder,
            X=ROI_data,
            y=targets,
            groups=runs_for_permutation,
            cv=ps,
            n_permutations=N_PERMS,
            random_state=rng_seed,
            n_jobs=N_PROC,
            verbose=3)
        accuracy = res[0]
        null_distribution = res[1]
        p_value = res[2]

    decode_accuracy.append(accuracy-1/len(set(targets)))
    decode_p_value.append(p_value)
    plt.hist(null_distribution,bins=50)
    plt.show()

    with h5py.File(output_dir, 'w') as f:
        f.create_dataset('accuracy', data=accuracy)
        if N_PERMS > 0:
            f.create_dataset('null_distribution', data=null_distribution)
            f.create_dataset('p_value', data=p_value)
            
    # let the programm 'sleep' for some time so cpu_usage is calculated correctly
    time.sleep(5)
    
del label_df
del betas

# FOURTH STEP
decode_accuracy = np.array(decode_accuracy)
decode_p_value  = np.array(decode_p_value)

x = np.arange(len(decode_accuracy))
ps = decode_p_value<0.05

fig = plt.figure()
plt.bar(x, decode_accuracy)
plt.plot(x[ps],decode_accuracy[ps]+0.1,'*')
plt.xticks(np.arange(len(decode_accuracy)),ROIS[:len(decode_accuracy)],rotation=45)
d = '_'
d = d.join(decoder_parameters.split(os.sep))
fig.savefig(os.path.join(RESULTS_DIR,SUB + '_' + DECODER + '_' + GLM_DATA_DIR + '_' + d +'.png'))

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
with open(os.path.join(RESULTS_DIR,'logfile.txt'), 'w+') as writer:
    writer.write('Codeversion: {} \n'.format(git_hash))
    writer.write('Random seed: {} \n'.format(rng_seed))
    writer.write('Number of permutations: {}\n'.format(N_PERMS))
    writer.write('Scaling: {}\n'.format(SCALE))
    writer.write('Cutoff: {}\n'.format(str(CUTOFF)))
    writer.write('Decoder used: {}\n'.format(DECODER))
    writer.write('Smoothing kernel of data: {}\n'.format(str(SMOOTHING_SIZE)))
    writer.write('Number of kernels used: {}\n'.format(str(N_PROC)))
    writer.write('Time for computation: {}h'.format(str((time.time() - T_START)/3600)))