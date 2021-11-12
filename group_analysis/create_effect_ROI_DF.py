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
# The purpuse of this script is to create a data frame that contains the mean
# values of all voxels within a set of ROIs.

################################
# FUNCTIONALITY OF THIS SCRIPT #
################################
# In order to do so, we first read in the SPM.mat file created by
# SPM12 when the GLM is estimated. From this SPM.mat file we read out the
# names of the beta NIfTI images and the names of the regressors, that 
# correspond to the beta image.
# This is done for every subject and with these information a pandas DataFrame
# is created looking like that: 
# Idx   Regressors          BetaNames       Subject_ID
# 0     Card_Appear_Magic   beta0001.nii    1
# 1     Card_Change_Control beta0002.nii    1
# .
# .
# .
# X     Constant            beta0180        24

# Then the run number (read out from the Regressors) is added, als well as the
# object, magic effect, if it was pre or post revelation and the video type 
# (Magic, Control, Surprise)
# Before any next step, the regressors of NO interest (realignment, constant)
# are removed
# In the final step we go through all ROIs (e.g. V1,V2 etc.) and within one ROI
# we iterate over all subjects. From the subject we read in the subject 
# specific ROI mask. This ROI mask is applied on all beta images of interest 
# and for every masked beta image a mean value is calculated. This long list 
# is then added to the dataframe.
# Idx   Regressors          BetaNames       Subject_ID  V1      V2      ...
# 0     Card_Appear_Magic   beta0001.nii    1           0.25    0.1
# 1     Card_Change_Control beta0002.nii    1           -0.32   0.7
# .
# .
# .
# X     Ball_Surprise3      beta0180        24          0.11    -2.1
#
# The resulting dataframe is saved in a hdf5/csv file

######################
# COMMAND LINE FLAGS #
######################
# --smooth, const=0, default=0 
# if smoothed data is used and if so what smoothing kernel 
# --analyzed, const='moment',  default='moment'
# The GLM was either calculated using all scans during one trial or depending
# on a specific moment during each video. This flag decides which of the GLMs
# is used

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
import pandas as pd  # to create data frames
import numpy as np   # most important numerical calculations
# read in mat files
import readmat
# needed to extract the run number out of the parentesis of the string in the SPM.mat file
import re
# library for neuroimaging
import nibabel as nib
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
PROJ_DIR        = os.path.join(HOME, 'Documents/Magic_fMRI/DATA/MRI')
RAWDATA_DIR     = os.path.join(PROJ_DIR, 'rawdata')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')
RESULTS_DIR     = os.path.join(DERIVATIVES_DIR, 'univariate-ROI',
                               data_analyzed, 'MagicEffects')
FREESURFER_DIR  = os.path.join(DERIVATIVES_DIR, 'freesurfer')
if SMOOTHING_SIZE > 0:
    GLM_DATA_DIR    = str(SMOOTHING_SIZE)+'mm-smoothed-nativespace' 
    FLA_DIR         = os.path.join(DERIVATIVES_DIR,'spm12',
                               'spm12-fla','WholeBrain',
                               'MagicEffects',GLM_DATA_DIR,
                               data_analyzed)
else:
    GLM_DATA_DIR    = 'nativespace' 
    FLA_DIR         = os.path.join(DERIVATIVES_DIR,'spm12',
                               'spm12-fla','WholeBrain',
                               'MagicEffects',GLM_DATA_DIR,
                               data_analyzed)
if not os.path.isdir(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)
    
SUBJECTS = glob.glob(os.path.join(FLA_DIR,'sub-*'))
SUBJECTS.sort()

# define ROIs
ROIS = [
        'V1', 'V2', 'V3', 'hV4', 
        'V3A', 'V3B', 
        'LO', 'VO', 
        'FEF', 'IPS',
        'ACC', 'PCC', 
        'IFG', 'aINSULA', 
        'IFJ', 'PHT', 'PF',
        'DMN', 'DAN', 'VAN', 'visual'
      ]

# objects used on magic tricks
OBJECTS = [
    'Ball',
    'Card',
    'Stick'
    ]

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

# empty lists to store the regressor names, names of beta images and subjec
# IDs
regressors  = []
beta_dirs   = []
ids         = []

# iterate over subjects to read in ALL the data in one huge data frame
for s, sub in enumerate(SUBJECTS):

    ########################################
    # reading in the necessary information #
    ########################################
    sub_mat_dir = os.path.join(FLA_DIR, sub, 'SPM.mat')

    # From the previously created SPM.mat file we read in the information we need
    # The filenames of our beta images
    sub_beta_dict = readmat.load(sub_mat_dir, isStruct=True)['SPM']['Vbeta']
    sub_beta_names = [f['fname'] for f in sub_beta_dict]
    # The corresponding Regressor names - are in form of 'Sn(<run-number>) <Regressor-Name>*bf(1)'
    sub_regressors = readmat.load(sub_mat_dir,isStruct=True)['SPM']['xX']['name']
    sub_id = [s+1]*len(sub_regressors)
    
    regressors.extend(sub_regressors)
    beta_dirs.extend(sub_beta_names)
    ids.extend(sub_id)

# store beta filenames and regressornames in a dictionary
data_dict = {
    'Regressors': regressors,
    'BetaNames': beta_dirs,
    'Subject_ID': ids
}

# convert dictionary into a pandas DataFrame for further analysis
label_df = pd.DataFrame(data_dict, columns=data_dict.keys())

# This complex loop is necessary to get the run number out of the regressor name
x       = [' '.join(re.findall(r"\((\d+)\)",string)) 
           for string in label_df.Regressors]
runs    = [int(s_filter.split()[0]) for s_filter in x]

# add further data to DataFrame
label_df['Runs']    = runs                  # In which run
pre_post_inblock    = (label_df.Runs-1)//2  # Resulting in 0,0,1,1,2,2,3,3 ...
label_df['PrePost'] = pre_post_inblock%2    # Resulting in 0,0,1,1,0,0,1,1 ...
label_df['Objects'] = np.nan                # Objects used in the videos
label_df['Effect']  = np.nan                # Which magic effect was performed
label_df['Type']    = np.nan                # is it magic, control or surprise

# again a complex process to throw out regressors of no interest (like realignment)
regressors_of_interest  = [True if ('Magic' in n) 
                           or ('Control' in n) 
                           or ('Surprise' in n) 
                           else False for n in label_df.Regressors]
# throw out all rows of regressors of no interest
label_df = label_df.iloc[regressors_of_interest]

# Check for every entry in Regressors if it contains one of the effect names. 
# If so, assign the effect name
for e in EFFECTS:
    label_df.Effect = np.where(label_df.Regressors.str.contains(e),e,label_df.Effect)
    
# Check for every entry in Regressors if it contains one of the object names. 
# If so, assign the object name
for o in OBJECTS:
    label_df.Objects = np.where(label_df.Regressors.str.contains(o),o,label_df.Objects)

# Check for every entry in Regressors if it contains one of the type names. 
# If so, assign the type name
for t in VIDEO_TYPE:
    label_df.Type = np.where(label_df.Regressors.str.contains(t),t,label_df.Type)

# inner loop - iterating over mask (=ROIs)
for r, roi in enumerate(ROIS):
    # empty list, that that will be filled with the mean of all voxels within
    # a ROI for all subjects
    roi_beta_means = []
    
    # iterate over all subjects
    for s, sub in enumerate(SUBJECTS):
        # get the directory where the subject specific ROI mask nifti immage
        # is stored
        roi_dir = os.path.join(FREESURFER_DIR,'sub-{:02d}'.format(s+1),
                               'corrected_ROIs')
    
        # Get all NifTi files containing the name of your ROI. Read them in and 
        # combine them to one ROI
        maskdir_list = glob.glob(os.path.join(roi_dir,'*' + roi + '*.nii'))
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
        
        # All beta Names of the current subject
        tpm_beta_series = label_df.BetaNames[label_df.Subject_ID==s+1]
        
        # iterate over the subjects beta images, read them in, calculate the
        # mean of the ROI in that image and append the mean to the list
        for b, beta in enumerate(tpm_beta_series):
            beta_nii    = nib.load(os.path.join(sub,beta))  # read in beta NIfTI image
            beta_data   = beta_nii.get_fdata()                  # get data from NIfTI image
            beta_nii.uncache()                                  # remove beta image from memory
            beta_data   = beta_data.flatten()                   # convert into one-dimensional array
            beta_roi    = beta_data[ROI & ~np.isnan(beta_data)]
            roi_beta_means.append(beta_roi.mean())
    
    # create a new column in the data frame with the name of the current ROI
    # and fill the column with the mean beta values of this ROI
    label_df[roi] = roi_beta_means
    
# save data as hdf5 and as csv file
label_df.to_hdf(os.path.join(RESULTS_DIR,'data_frame.hdf5'),key='df',mode='w')
label_df.to_csv(path_or_buf=os.path.join(RESULTS_DIR,'data_frame.csv'))

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
with open(os.path.join(RESULTS_DIR,'meanDF-logfile.txt'), 'w+') as writer:
    writer.write('Codeversion: {} \n'.format(git_hash))
    writer.write('Smoothing kernel used: {}\n'.format(str(SMOOTHING_SIZE)))
    writer.write('Time for computation: {}h'.format(str((time.time() - T_START)/3600)))
