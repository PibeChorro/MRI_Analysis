#!/gpfs01/bartels/user/vplikat/anaconda3/bin/python

# Script to load in ROIs and test if they have overlapping voxels. If so, try to assign them to one ROI

#############
# LIBRARIES #
#############

# libraries to interact with the operating system 
import os
import sys
from pathlib import Path
import glob
import numpy as np   # most important numerical calculations
# library for neuroimaging
import nibabel as nib
from sklearn.neighbors import KNeighborsClassifier

# define NIfTI image shape
NIfTI_shape = (96, 96, 62)

sub = sys.argv[1]

# define ROIs
all_ROIs = ['V1', 'V2', 'V3', 'hV4', 
        'V3A', 'V3B', 
        'LO1', 'LO2', 
        'VO1', 'VO2', 
        'TO1', 'TO2', 
        'FEF', 'SPL1',
        'PHC1', 'PHC2',
        'IPS1' ,'IPS2', 'IPS3', 'IPS4','IPS5']

# declare all path names
home = str(Path.home())
surfer_dir = os.path.join (home, 'Documents/Master_Thesis/DATA/MRI/derivatives/freesurfer/')    # the freesurfer directory              # 

ROI_dir = os.path.join(surfer_dir,sub,'ROIs')
corrected_ROI_dir = os.path.join(surfer_dir,sub,'corrected_ROIs')
if not os.path.exists(corrected_ROI_dir):
    os.makedirs(corrected_ROI_dir)

# iterate over all ROIs except the last one, those will be ROI0
for roi0, ROI0 in enumerate (all_ROIs[:-1]):
    # iterate over all following ROIs, those will be ROI1
    for roi1, ROI1 in enumerate (all_ROIs[roi0+1:]):

        # get a list containing all NIfTI files belonging to the two ROIs
        ROI0_list = glob.glob(os.path.join(ROI_dir,'*'+ROI0+'.nii'))
        ROI1_list = glob.glob(os.path.join(ROI_dir,'*'+ROI1+'.nii'))
        
        # check if the current ROIs already have been changed. If so, select the changed one
        if os.path.exists(os.path.join (corrected_ROI_dir,ROI0 + '.nii')):
            ROI0_list = glob.glob(os.path.join(corrected_ROI_dir,ROI0 +'.nii'))
        if os.path.exists(os.path.join (corrected_ROI_dir,ROI1 + '.nii')):
            ROI1_list = glob.glob(os.path.join(corrected_ROI_dir,ROI1 +'.nii'))

        # create an empty list for both ROIs. For each ROI we have at least 2 nifti files (for each hemisphere)
        # Some are combined (like V1d and V1v resulting in 4 nifti files). We take all files and add them up
        ROI0_matrix = np.zeros(shape=NIfTI_shape)
        ROI1_matrix = np.zeros(shape=NIfTI_shape)

        # read in the nifti files for the first ROI 
        for part in ROI0_list:
            #print ('Reading in {}\n'.format(part))
            R0_mask_nii = nib.load(part)
            mask_img = R0_mask_nii.get_fdata()
            mask_img = np.asarray(mask_img)
            ROI0_matrix = np.add(ROI0_matrix,mask_img)

        # read in the nifti files for the first ROI 
        for part in ROI1_list:
            #print ('Reading in {}\n'.format(part))
            R1_mask_nii = nib.load(part)
            mask_img = R1_mask_nii.get_fdata()
            mask_img = np.asarray(mask_img)
            ROI1_matrix = np.add(ROI1_matrix,mask_img)


        # convert into boolian matrix
        ROI0_matrix = ROI0_matrix>0
        ROI1_matrix = ROI1_matrix>0
        print ('Number of Voxels in {}: {}'.format(ROI0, ROI0_matrix.sum()))
        print ('Number of Voxels in {}: {}'.format(ROI1, ROI1_matrix.sum()))
        
        # look for overlapping voxels
        overlap = np.bitwise_and(ROI0_matrix,ROI1_matrix)
        if overlap.sum() == 0:
            print('{} and {} do not have any overlapping voxels.\n Continue with the next ROI'.format(ROI0,ROI1))
            # Create new clean ROI mask
            # first convert to a form needed
            ROI_mask0 = nib.Nifti1Image(ROI0_matrix.astype(float), R0_mask_nii.affine)
            ROI_mask1 = nib.Nifti1Image(ROI1_matrix.astype(float), R1_mask_nii.affine)
            # second, write to new location
            nib.save(ROI_mask0, os.path.join(corrected_ROI_dir,ROI0 + '.nii'))
            nib.save(ROI_mask1, os.path.join(corrected_ROI_dir,ROI1 + '.nii'))
            continue
            
        else:
            print ('{} and {} do have {} overlapping voxels.\n These will now be assigned'.format(ROI0,ROI1,overlap.sum()))
        

        # set voxels tha overlap to False in both ROIs
        indices_of_interest = np.where(overlap)
        ROI0_matrix[indices_of_interest] = False
        ROI1_matrix[indices_of_interest] = False

        # Combine both ROIs in one matrix and create a label matrix containing the labeling values for the ROIs
        combined_ROIs = np.add(ROI0_matrix,ROI1_matrix)
        train_data = np.where(combined_ROIs)
        label = np.zeros_like(combined_ROIs,dtype=int)
        label[np.where(ROI0_matrix)] = 1
        label[np.where(ROI1_matrix)] = 2

        # cut off all False and 0 values
        label = label[label>0]

        # np.where returns a tuple and we need an array
        train_data_arr = np.asarray(train_data)

        # iniciate the classifier with k=10 and weight the distance
        my_knn = KNeighborsClassifier(n_neighbors=10, weights="distance")
        # fit the classifier on the unambiguous data points (transposing the array is necessary)
        my_knn.fit(train_data_arr.T,label)

        # get the ambiguous datapoints
        test_data = np.where(overlap)
        # again turn into numpy array
        test_data_arr = np.asarray(test_data)
        # let the classifier predict the right labels for the ambiguous datapoints (transposing is necessary)
        result=my_knn.predict(test_data_arr.T)
        # labels are '1' and '2' by logical comparison we get false values for ROI0 and true values for ROI1
        result=result>1

        # set the ambiguous values right
        ROI0_matrix[test_data]= ~result
        ROI1_matrix[test_data]= result
        
        overlap = np.bitwise_and(ROI0_matrix,ROI1_matrix)
        if overlap.sum() == 0:
            # Create new clean ROI mask
            # first convert to a form needed
            ROI_mask0 = nib.Nifti1Image(ROI0_matrix.astype(float), R0_mask_nii.affine)
            ROI_mask1 = nib.Nifti1Image(ROI1_matrix.astype(float), R1_mask_nii.affine)
            # second, write to new location
            nib.save(ROI_mask0, os.path.join(corrected_ROI_dir,ROI0 + '.nii'))
            nib.save(ROI_mask1, os.path.join(corrected_ROI_dir,ROI1 + '.nii'))
            print('Assigning was successful. {} and {} do not have any overlaying voxels anymore.'.format(ROI0,ROI1))
            continue
            
        else:
            print ('Something went wrong. {} and {} still have {} overlapping voxels.\n Abborting here'.format(ROI0,ROI1,overlap.sum))
            raise Exception ('Abborted!')
