#!/gpfs01/bartels/user/vplikat/anaconda3/bin/python

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

################################
# Command line input arguments #
################################

# create a parser object which handles the input variables
parser = argparse.ArgumentParser()

# add all the input arguments
parser.add_argument("--sub",        "-s",                               default='sub-01')         # subject
parser.add_argument("--smooth",             nargs='?',  const=0,        default=0, type=int)         # what data should be used
parser.add_argument("--algorythm",  "-a",   nargs='?',  const='LDA',    default='LDA')
parser.add_argument("--scaling",            nargs='?',  const='None',   default='None', type=str)
parser.add_argument("--cutoff",     "-c",   nargs='?',  const=np.inf,   default=np.inf, type=float) # if and with which value (in std) data is cut off 
parser.add_argument("--feature",    "-f",   nargs='?',  const='None',   default='None', type=str)
parser.add_argument("--kernels",    "-k",   nargs='?',  const=12,       default=1, type=int)   # how many processes should be run in parallel
# parse the arguments to a parse-list(???)
ARGS = parser.parse_args()
# assign values 
SUB             = ARGS.sub
N_PROC          = ARGS.kernels
SMOOTHING_SIZE  = ARGS.smooth
DECODER         = ARGS.algorythm
CUTOFF          = ARGS.cutoff
FEAT_TRANS      = ARGS.feature
SCALE           = ARGS.scaling

# make sure if scaling is None that cutoff value is inf and there is no need for a cutoff value level in file hierarchy
if SCALE == 'None':
    CUTOFF = np.inf
    decoder_parameters = 'scale_'+SCALE
else:
    decoder_parameters = os.path.join('scale_'+SCALE,'cutoff_'+str(CUTOFF))
    
    
if DECODER =='LDA':
    my_decoder          = LDA(solver='lsqr', shrinkage='auto')
    decoder_parameters  = os.path.join(
        decoder_parameters,
        'feat_trans_'+FEAT_TRANS)
elif DECODER == 'SVM':
    SVM_C = 1
    decoder_parameters  = os.path.join(
        decoder_parameters,
        'feat_trans_'+FEAT_TRANS,
        'C_'+str(SVM_C))
    if FEAT_TRANS == 'PCA':
        my_decoder          = SVC(kernel='rbf', C=SVM_C)
    else: 
        my_decoder          = SVC(kernel='linear', C=SVM_C)
# make LDA the default in case something completely different was given
else:
    my_decoder          = LDA(solver='lsqr', shrinkage='auto')
    decoder_parameters  = os.path.join(
        decoder_parameters,
        'feat_trans_'+FEAT_TRANS)
        

# variables for path selection and data access
HOME            = str(Path.home())
PROJ_DIR        = os.path.join(HOME, 'Documents/Master_Thesis/DATA/MRI')
DERIVATIVES_DIR = os.path.join(PROJ_DIR, 'derivatives')
if SMOOTHING_SIZE > 0:
    GLM_DATA_DIR    = str(SMOOTHING_SIZE)+'mm-smoothed-nativespace' 
    FLA_DIR         = os.path.join(DERIVATIVES_DIR,'spm12',
                               'spm12-fla','WholeBrain',
                               'EveryVideo',GLM_DATA_DIR,
                               'WholeVideo',SUB)
else:
    GLM_DATA_DIR    = 'nativespace' 
    FLA_DIR         = os.path.join(DERIVATIVES_DIR,'spm12',
                               'spm12-fla','WholeBrain',
                               'EveryVideo',GLM_DATA_DIR,
                               'WholeVideo',SUB)
FREESURFER_DIR  = os.path.join(DERIVATIVES_DIR, 'freesurfer')
RAWDATA_DIR     = os.path.join(PROJ_DIR, 'rawdata')
ROI_DIR         = os.path.join(FREESURFER_DIR,SUB,'corrected_ROIs')
SPM_MAT_DIR     = os.path.join(FLA_DIR, 'SPM.mat')
ANALYSIS        = 'ROI-analysis'
RESULTS_DIR     = os.path.join(DERIVATIVES_DIR, 'decoding','object_decoding', 'train75_test25_notmagic', ANALYSIS, 
                               DECODER, GLM_DATA_DIR, decoder_parameters, SUB)
if not os.path.isdir(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

#define ROIs
ROIS = [
        'ventricle',
        'V1', 'V2', 'V3', 'hV4', 
        'V3A', 'V3B', 
        'LO1', 'LO2', 
        'VO1', 'VO2', 
        'TO1', 'TO2', 
        'FEF', 'SPL1',
        'PHC1', 'PHC2',
        'IPS1' ,'IPS2', 'IPS3', 'IPS4','IPS5'
      ]

LABEL_NAMES = [
    'Ball',
    'Card',
    'Stick'
]

# Get effect names to train and test on them using cross validation
EFFECT_NAMES = [
    'Appear',
    'Change',
    'Vanish'
]

# empty lists that will be filled with the results to plot after calculation
decode_accuracy = []
decode_p_value  = []

# Optional arguments
rng_seed = 0
n_permutations = 1000
    
print ('Analysing subject: {}'.format(SUB))
print ('Getting ROIs from:	 {}'.format(FREESURFER_DIR))
print ('Saving data at:	 {}'.format(RESULTS_DIR))
print ('rng_seed:	 {}'.format(rng_seed))
print ('Number of permutations:	 {}'.format(n_permutations))


########################################
# reading in the necessary information #
########################################

# From the previously created SPM.mat file we read in the information we need
# The filenames of our beta images
SPM_BETADICT = readmat.load(SPM_MAT_DIR, isStruct=True)['SPM']['Vbeta']
BETA_DIRS = [f['fname'] for f in SPM_BETADICT]
# The corresponding Regressor names - are in form of 'Sn(<run-number>) <Regressor-Name>*bf(1)'
SPM_REGRESSORS = readmat.load(SPM_MAT_DIR,isStruct=True)['SPM']['xX']['name']

# store beta filenames and regressornames in a dictionary
data_dict = {
    'Regressors': SPM_REGRESSORS,
    'BetaNames': BETA_DIRS
}

# convert dictionary into a pandas DataFrame for further analysis
label_df = pd.DataFrame(data_dict, columns=data_dict.keys())

# This complex loop is necessary to get the run number out of the regressor name
x       = [' '.join(re.findall(r"\((\d+)\)",string)) for string in label_df.Regressors]
runs    = [int(s_filter.split()[0]) for s_filter in x]

# add further data to DataFrame
label_df['Runs']    = runs   # In which run
# again a complex list comprehension to only keep magic videos
regressors_of_interest  = [any(i in n for i in LABEL_NAMES) & ('Magic' not in n) for n in SPM_REGRESSORS]
# throw out all rows of regressors of no interest
label_df = label_df.iloc[regressors_of_interest]

label_df['Labels']  = np.nan # Labels
label_df['Chunks']  = label_df.Runs%4 # Chunks needed for cross validation
# Check for every entry in Regressors if it contains one of the label names. If so, assign the label name
for l in LABEL_NAMES:
    label_df.Labels = np.where(label_df.Regressors.str.contains(l),l,label_df.Labels)
    #label_df.Labels[label_df.Regressors.str.contains(l)] = l
    

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
    print ('Doing ROI {}'.format(roi))
    output_dir = os.path.join(RESULTS_DIR,roi + '.hdf5')   # where to store the results

    maskdir_list = glob.glob(os.path.join(ROI_DIR,'*' + roi + '*.nii'))
    masklist = []
    for mask in maskdir_list:
        mask_nii = nib.load(mask)
        mask_img = mask_nii.get_fdata()
        mask_nii.uncache()
        mask_img = np.asarray(mask_img)
        mask_img = mask_img.flatten()
        masklist.append(mask_img)

    ROI = np.sum(masklist,axis=0)
    ROI = ROI>0
    
    # get data you want to use in the ROI analysis
    ROI_data = []
    ROI_data.append([beta[ROI & ~np.isnan(beta)] for beta in betas])
    ROI_data = np.array(ROI_data [0])               # the list comprehension wraps the matrix in an additional, unnecessary array
    # convert data into z values and cut off data
    
    if SCALE == 'z':
        # z-transform data within features
        ROI_data = ROI_data - ROI_data.mean(axis=0)
        ROI_data = ROI_data / ROI_data.std(axis=0)
    elif SCALE == 'min0max1':
        ROI_data = ROI_data - ROI_data.min(axis=0)
        ROI_data = ROI_data / ROI_data.max(axis=0)
    elif SCALE == 'mean':
        ROI_data = ROI_data - ROI_data.mean(axis=0)
        
    if FEAT_TRANS == 'PCA':
        n_components    = min ([ROI_data.shape[1], ROI_data.shape[0]//2])   # take maximally half of the data points as used dimensions
        my_PCA          = PCA(n_components=n_components)
        my_PCA.fit(ROI_data)
        ROI_data        = my_PCA.transform(ROI_data)
        #ROI_data        = ROI_data[:,my_PCA.explained_variance_ratio_>0.03]
    
    # set outliers to CUTOFF value
    ROI_data[ROI_data>CUTOFF]   = CUTOFF
    ROI_data[ROI_data<-CUTOFF]  = -CUTOFF
    
    # the actual decoding
    targets = np.asarray(label_df.Labels)
    chunks = np.asarray(label_df.Chunks)
    if n_permutations > 0:
        res = permutation_test_score(
            estimator=my_decoder,
            X=ROI_data,
            y=targets,
            scoring="accuracy",
            groups=chunks,
            cv=PredefinedSplit(chunks),
            n_permutations=n_permutations,
            random_state=rng_seed,
            n_jobs=N_PROC,
            verbose=3)
        accuracy = res[0]
        null_distribution = res[1]
        p_value = res[2]

    decode_accuracy.append(accuracy-1/3)
    decode_p_value.append(p_value)

    with h5py.File(output_dir, 'w') as f:
        f.create_dataset('accuracy', data=accuracy)
        if n_permutations > 0:
            f.create_dataset('null_distribution', data=null_distribution)
            f.create_dataset('p_value', data=p_value)
            
    
    # let the programm 'sleep' for some time so cpu_usage is calculated correctly
    time.sleep(5)
    
del label_df
del betas

decode_accuracy = np.array(decode_accuracy)
decode_p_value = np.array(decode_p_value)

x = np.arange(len(decode_accuracy))
ps = decode_p_value<0.05

fig = plt.figure()
plt.bar(x, decode_accuracy)
plt.plot(x[ps],decode_accuracy[ps]+0.1,'*')
plt.xticks(np.arange(len(ROIS)),ROIS,rotation=45)
d = '_'
d = d.join(decoder_parameters.split(os.sep))
fig.savefig(os.path.join(RESULTS_DIR,SUB + '_' + DECODER + '_' + GLM_DATA_DIR + '_' + d +'.png'))

rep = git.Repo(search_parent_directories=True)
git_hash = rep.head.object.hexsha

# create a log file, that saves some information about the run script
with open(os.path.join(RESULTS_DIR,'logfile.txt'), 'w+') as writer:
    writer.write('Codeversion: {} \n'.format(git_hash))
    writer.write('Number of permutations: {}\n'.format(n_permutations))
    writer.writer('Number of kernels used: {}\n'.format(str(N_PROC)))
    writer.write('Time for computation: {}h'.format(str((time.time() - T_START)/3600)))

# print time the whole processe took
print ((time.time() - T_START)/3600)