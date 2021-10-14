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
# The purpose of this script is to read in the beta estimates for every video
# presented for a set of ROIs and add the behavioural rating (TO-DO add even
# the pupil data).

################################
# FUNCTIONALITY OF THIS SCRIPT #
################################
# In order to do so, we first read in the SPM.mat file created by
# SPM12 when the GLM is estimated. From this SPM.mat file we read out the
# names of the beta NIfTI images and the names of the regressors, that 
# correspond to the beta image.
# At the same time the behavioral data from subjects are read in and added 
# This is done for every subject and with these information a pandas DataFrame
# is created looking like that: 
# Idx   Regressors          BetaNames       Subject_ID  Rating
# 0     Card_Appear_Magic   beta0001.nii    1           4
# 1     Card_Change_Control beta0002.nii    1           1
# .
# .
# .
# X     Constant            beta0180        24          nan

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
# Idx   Regressors          BetaNames       Subject_ID      V1      V2      ...
# 0     Card_Appear_Magic   beta0001.nii    1           4   0.25    0.1
# 1     Card_Change_Control beta0002.nii    1           1   -0.32   0.7
# .
# .
# .
# X     Ball_Surprise3      beta0180        24          3   0.11    -2.1
#
# Finally for all subjects the events.tsv files of the 12 runs are read in and 
# surprise ratings are added to the data_frame
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
PROJ_DIR        = os.path.join(HOME, 'Documents/Master_Thesis/DATA/MRI')
RAWDATA_DIR     = os.path.join(PROJ_DIR, 'rawdata')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')
RESULTS_DIR     = os.path.join(DERIVATIVES_DIR, 'univariate-ROI',
                               data_analyzed,'EveryVideo')
FREESURFER_DIR  = os.path.join(DERIVATIVES_DIR, 'freesurfer')
if SMOOTHING_SIZE > 0:
    GLM_DATA_DIR    = str(SMOOTHING_SIZE)+'mm-smoothed-nativespace' 
    FLA_DIR         = os.path.join(DERIVATIVES_DIR,'spm12',
                               'spm12-fla','WholeBrain',
                               'EveryVideo',GLM_DATA_DIR,
                               data_analyzed)
else:
    GLM_DATA_DIR    = 'nativespace' 
    FLA_DIR         = os.path.join(DERIVATIVES_DIR,'spm12',
                               'spm12-fla','WholeBrain',
                               'EveryVideo',GLM_DATA_DIR,
                               data_analyzed)
if not os.path.isdir(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)
    
SUBJECTS_WHOLE_PATH = glob.glob(os.path.join(FLA_DIR,'sub-*'))
SUBJECTS_WHOLE_PATH.sort()
SUBJECTS = [os.path.basename(sub) for sub in SUBJECTS_WHOLE_PATH]

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
data_dict = {
'regressors': [],
'beta_dirs': [],
'ids': [],
'runs': [],
    }
for roi in ROIS:
    data_dict[roi] = []
ratings = []

# iterate over subjects to read in ALL the data in one huge data frame
for s, sub in enumerate(SUBJECTS[0:3]):

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
    # This complex loop is necessary to get the run number out of the regressor name
    x       = [' '.join(re.findall(r"\((\d+)\)",string)) for string in sub_regressors]
    sub_run = [int(s_filter.split()[0]) for s_filter in x]
    
    data_dict['regressors'].extend(sub_regressors)
    data_dict['beta_dirs'].extend(sub_beta_names)
    data_dict['ids'].extend(sub_id)
    data_dict['runs'].extend(sub_run)
    
    # read in the beta images from current subject and calculate its mean for 
    # each roi
    
    for roi in ROIS:
        # get ROI mask
        roi_dir = os.path.join(FREESURFER_DIR,sub, 'corrected_ROIs')
    
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
            
        betas = []
        for beta_file in sub_beta_names:
            beta_dir = os.path.join(FLA_DIR, sub, beta_file)
            beta_img = nib.load(beta_dir)
            beta_data = beta_img.get_fdata()
            beta_data = beta_data.flatten()
            beta_roi    = beta_data[ROI & ~np.isnan(beta_data)]
            betas.append(beta_roi.mean())
            
        data_dict[roi].extend(betas)
        
    
# turn dictionary into a pandas data frame
data_frame = pd.DataFrame(data=data_dict,columns=data_dict.keys())
# again a complex process to throw out regressors of no interest (like realignment)
regressors_of_interest  = [True if ('Magic' in n) 
                           or ('Control' in n) 
                           or ('Surprise' in n) 
                           else False for n in data_frame.regressors]
data_frame = data_frame.iloc[regressors_of_interest]

data_frame['Objects']   = np.nan                # Objects used in the videos
data_frame['Effect']    = np.nan                # Which magic effect was performed
data_frame['Type']      = np.nan                # is it magic, control or surprise
# data_frame['Rating']    = np.nan                # behavioral surprise rating
pre_post_inblock        = (data_frame.runs-1)//2
data_frame['PrePost']    = pre_post_inblock%2    # Labels

# Check for every entry in Regressors if it contains one of the effect names. 
# If so, assign the effect name
for e in EFFECTS:
    data_frame.Effect = np.where(data_frame.regressors.str.contains(e),e,data_frame.Effect)
    
# Check for every entry in regressors if it contains one of the object names. 
# If so, assign the object name
for o in OBJECTS:
    data_frame.Objects = np.where(data_frame.regressors.str.contains(o),o,data_frame.Objects)

# Check for every entry in regressors if it contains one of the type names. 
# If so, assign the type name
for t in VIDEO_TYPE:
    data_frame.Type = np.where(data_frame.regressors.str.contains(t),t,data_frame.Type)

###########################
# get the behavioral data #
###########################
# The behavioral data is stored in events.tsv files in the rawdata 
# directory. 
# iterate over subjects to read in ALL the data in one huge data frame
for s, sub in enumerate(SUBJECTS[0:3]):
    run_events = glob.glob(os.path.join(RAWDATA_DIR,sub,'func','*events.tsv'))
    run_events.sort()
    for event_file in run_events:
        event_df = pd.read_csv(filepath_or_buffer=event_file,delimiter='\t')
        ratings.extend(event_df.value[event_df.index%2==1])
        
data_frame['Ratings'] = ratings
data_frame.to_csv(path_or_buf=os.path.join(RESULTS_DIR,'ratings_df.csv'))


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
with open(os.path.join(RESULTS_DIR,'data_frame_creation-logfile.txt'), 'w+') as writer:
    writer.write('Codeversion: {} \n'.format(git_hash))
    writer.write('Time for computation: {}h'.format(str((time.time() - T_START)/3600)))
