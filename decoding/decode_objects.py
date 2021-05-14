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
import pandas as pd  # to create data frames
import numpy as np   # most important numerical calculations
# read in mat files
import readmat
# needed to extract the run number out of the parentesis of the string in the SPM.mat file
import re
# library for neuroimaging
import nibabel as nib
# machine learning algorithms and stuff
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
parser.add_argument("--sub", "-s",default='sub-01')         # subject
parser.add_argument("--kernels", "-k",default=1,type=int)   # how many processes should be run in parallel
parser.add_argument("--smooth", default=0,type=int)         # what data should be used
parser.add_argument("--algorythm", "-a", default='LDA')
# parse the arguments to a parse-list(???)
ARGS = parser.parse_args()
# assign values 
SUB = ARGS.sub
N_PROC = ARGS.kernels
SMOOTHING_SIZE  = ARGS.smooth
DECODER = ARGS.algorythm
if DECODER =='LDA':
    my_decoder          = LDA(solver='lsqr', shrinkage='auto')
    decoder_parameters  = ''
elif DECODER == 'SVM':
    SVM_C = 0.0001
    decoder_parameters  = 'C_{}'.format(SVM_C)
    my_decoder          = SVC(kernel='linear', C=SVM_C)
        

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
RESULTS_DIR     = os.path.join(DERIVATIVES_DIR, 'decoding', ANALYSIS, DECODER, decoder_parameters, GLM_DATA_DIR, SUB)
if not os.path.isdir(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

# define ROIs
ROIS = [
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

# empty lists that will be filled with the results to plot after calculation
decode_accuracy = []
decode_p_value = []

# Optional arguments
rng_seed = 0
n_permutations = 1000
    
print ('SUB: {}'.format(SUB))
print ('surfer_dir:	 {}'.format(FREESURFER_DIR))
print ('RESULTS_DIR:	 {}'.format(RESULTS_DIR))
print ('rng_seed:	 {}'.format(rng_seed))
print ('n_perm:	 {}'.format(n_permutations))


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
label_df['Runs']    = runs                # In which run
label_df['Chunks']  = label_df.Runs%2     # The chunks (needed for cross validation)
label_df['Labels']  = np.nan              # Labels
# Check for every entry in Regressors if it contains one of the label names. If so, assign the label name
for l in LABEL_NAMES:
    label_df.Labels[label_df.Regressors.str.contains(l)] = l

# again a complex process to throw out regressors of no interest (like realignment)
regressors_of_interest = [True if any(i in n for i in LABEL_NAMES) else False for n in SPM_REGRESSORS]
# throw out all rows of regressors of no interest
label_df = label_df.iloc[regressors_of_interest]

# inner loop - iterating over mask (=ROIs)
for r, roi in tqdm(enumerate(ROIS)):
    output_dir = os.path.join(RESULTS_DIR,roi)   # where to store the results
    #if not os.path.isdir(output_dir):
    #    os.mkdir(output_dir)

    # call combineROIs with the selected ROI and ROI directory
    #ROI = combineROIs(roi, ROI_dir)
    maskdir_list = glob.glob(os.path.join(ROI_DIR,'*' + roi + '*.nii'))
    masklist = []
    for mask in maskdir_list:
        mask_nii = nib.load(mask)
        mask_img = mask_nii.get_fdata()
        mask_img = np.asarray(mask_img)
        mask_img = mask_img.flatten()
        masklist.append(mask_img)

    ROI = np.sum(masklist,axis=0)
    ROI = ROI>0
    # read in all beta image files, convert them into a one-dimensional numpy array and select the entries of ROI
    betas = []                                              # empty list to store data arrays in
    for b, beta in enumerate(label_df.BetaNames):
        beta_nii = nib.load(os.path.join(FLA_DIR,beta)) # read in beta NIfTI image
        beta_data = beta_nii.get_fdata()                    # get data from NIfTI image
        beta_data = beta_data.flatten()                     # convert into one-dimensional array
        beta_data = beta_data[ROI]                          # select entries from ROI
        beta_data = beta_data[~np.isnan(beta_data)]         # throw out all NaN values. This can happen, when the mask selects Voxels that are on the skull or skin
        betas.append(beta_data)                             # append array on betas list

    # convert list into numpy array 
    betas = np.array(betas)

    # the actual decoding
    targets = np.asarray(label_df.Labels)
    chunks = np.asarray(label_df.Chunks)
    if n_permutations > 0:
        res = permutation_test_score(
            my_decoder,
            betas,
            targets,
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
    time.sleep(30)
    
del label_df
del betas

decode_accuracy = np.array(decode_accuracy)
decode_p_value = np.array(decode_p_value)

x = np.arange(len(decode_accuracy))
ps = decode_p_value<0.05

fig = plt.figure()
plt.bar(x, decode_accuracy)
plt.plot(x[ps],decode_accuracy[ps]+0.1,'*')
plt.xticks(np.arange(len(decode_accuracy)),ROIS,rotation=45)
fig.savefig(os.path.join(RESULTS_DIR,SUB+'.png'))

# print time the whole processe took
print ((time.time() - T_START)/3600)