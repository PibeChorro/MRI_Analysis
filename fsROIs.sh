# Script to use freesurfer for ROI analysis

# 0. Step: Set enviriounment variables
# 1. Step: Convert structural images from NIFTIs to mgz using mri_convert
# 2. Step: Inflate brains using recon-all
# 3. Step: create ROIs using Benson atlas
# 3a. Step: use neuropythy to create atlas-images. Resulting in benson and wang atlas images.
# 3b. Step: from the atlases create label images. The atlas assign a number 1-25 to each individual region
# 3c. Step: merge ROI labels of V1-V3
# 4. Step: ???
# 5. Step: create ROI NiFTIs 

# 0. Step
# MAKE SURE FREESURFER IS SET UP
export FREESURFER_HOME=/Users/vpl/Documents/freesurfer 
source $FREESURFER_HOME/SetUpFreeSurfer.sh 

# The directory where FREESURFER should store all the processed images and steps in between -- empty if run for the first time !!!
export SUBJECTS_DIR=$HOME/Documents/Master_Thesis/DATA/BIDS_playground/derivative_data/surfer

# The directory where the unprocessed NiFTI images are -- important for the anatomical scan
export raw_dir=$HOME/Documents/Master_Thesis/DATA/BIDS_playground/raw_data

# The directory where the corregistered functional data are stored -- needed for the mean EPI image
export corregistered_dir=$HOME/Documents/Master_Thesis/DATA/BIDS_playground/derivative_data/coregistered

# Those are the ROI names you just need to know them -- V1v=1, V1d=2 in the label images (... I guess ???)
# (https://hub.docker.com/r/nben/occipital_atlas) "older version" of neuropythy)
export roiname_array=("V1v" "V1d" "V2v" "V2d" "V3v" "V3d" "hV4" "VO1" "VO2" "PHC1" "PHC2" \
"TO2" "TO1" "LO2" "LO1" "V3B" "V3A" "IPS0" "IPS1" "IPS2" "IPS3" "IPS4" \
"IPS5" "SPL1" "FEF")

# ROIs we want to merge
export rois_to_merge=("V1" "V2" "V3")

# all ROIs combined -- this strange syntax is needed
allROIs=( "${roiname_array[@]}" "${rois_to_merge[@]}" )
export allROIs

# change into raw_data direcotry to get all subjects and create a folder for each in the FREESURFER SUBJECTS_DIR
cd $raw_dir

for subject in sMag*
do
	mkdir $SUBJECTS_DIR/$subject
	mkdir $SUBJECTS_DIR/$subject/mri
done

cd $SUBJECTS_DIR

for subject in sMag*
do
	# 1. Step
	#mri_convert $raw_dir/$subject/anat/*.nii $SUBJECTS_DIR/$subject/mri/001.mgz
	
	# 2. Step
	#recon-all -autorecon-all -subjid $subject
	
	# 3a. Step
	#python -m neuropythy atlas --verbos $subject

	# 3b. Step
	# iterate 25 times -- number of ROIs
	for i in {0..24}
	do
		# command: mri_cor2label
		#1 the wang atlas created by neuropythy atlas (two wang atlases are created ending on mplbl.mgz and fplbl.mgz. The second does not work)
		#2 the number the ROI is assigned to in the atlas image
		#3 output file name
		#4 subject with information about the hemisphere and if it is inflated or not (???)
		
		# LEFT HEMISPHERE
	 	mri_cor2label --i $SUBJECTS_DIR/$subject/surf/lh.wang15_mplbl.mgz\
	 				--id $(($i+1)) \
	 				--l lh.wang15atlas.${roiname_array[$i]}.label \
	 				--surf $subject lh inflated
	 				
	 	# RIGHT HEMISPHERE
		mri_cor2label --i $SUBJECTS_DIR/$subject/surf/rh.wang15_mplbl.mgz \
					--id $(($i+1)) \
					--l rh.wang15atlas.${roiname_array[$i]}.label \
					--surf $subject rh inflated
	done
	
	# 3c. Step
	# merge V1, V2 and V3
	for i in {0..2}
	do
		# command: mri_mergelabels
		# -i an imput image to merge with other input images
		# -o output image
		
		# LEFT HEMISPHERE
		mri_mergelabels -i ${SUBJECTS_DIR}/${subject}/label/lh.wang15atlas.${rois_to_merge[${i}]}v.label \
						-i ${SUBJECTS_DIR}/${subject}/label/lh.wang15atlas.${rois_to_merge[${i}]}d.label \
						-o ${SUBJECTS_DIR}/${subject}/label/lh.wang15atlas.${rois_to_merge[${i}]}.label 
						
		# RIGHT HEMISPHERE
		mri_mergelabels -i ${SUBJECTS_DIR}/${subject}/label/rh.wang15atlas.${rois_to_merge[${i}]}v.label \
						-i ${SUBJECTS_DIR}/${subject}/label/rh.wang15atlas.${rois_to_merge[${i}]}d.label \
						-o ${SUBJECTS_DIR}/${subject}/label/rh.wang15atlas.${rois_to_merge[${i}]}.label
	done
	
	
	# 4. Step
	# ASK PABLO WHAT THIS MEANS -- is this really a step after label creation???
	tkregister2 --mov $corregistered_dir/${subject}/func/mean*.nii --s ${subject} --regheader --reg $SUBJECTS_DIR/$subject/register_${subject}.dat
	
	# 5. Step
	# create a folder to store ROIs in
	mkdir $SUBJECTS_DIR/$subject/ROIs
	
	# get the meanEPI image
	meanNiFTI=$corregistered_dir/$subject/func/mean*.nii
	
	# convert labels into ROIs
	for roi in "${allROIs[@]}"
	do
		
		# command: mri_label2vol 
		#1: the label image as input
		#2: the mean EPI NiFTI 
		#3: the registration file from Step 4
		#4: ???
		#5: ???
		#6: subject
		#7: hemisphere
		#8: output directory
		
		# LEFT HEMISPHERE
		mri_label2vol --label $subject/label/lh.wang15atlas.$roi.label \
		--temp $meanNiFTI \
		--reg $SUBJECTS_DIR/$subject/register_${subject}.dat \
		--fillthresh 0 \
		--proj frac 0 1 0.1 \
		--subject $subject \
		--hemi lh \
		--o $SUBJECTS_DIR/$subject/ROIs/lh.wang15atlas.$roi.nii
		
		# RIGHT HEMISPHERE
		mri_label2vol --label $subject/label/rh.wang15atlas.$roi.label \
		--temp $meanNiFTI \
		--reg $SUBJECTS_DIR/$subject/register_${subject}.dat \
		--fillthresh 0 \
		--proj frac 0 1 0.1 \
		--subject $subject \
		--hemi rh \
		--o $SUBJECTS_DIR/$subject/ROIs/rh.wang15atlas.$roi.nii
	done
done

cd $HOME