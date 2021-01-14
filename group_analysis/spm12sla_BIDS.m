%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for the second level analysis of the "Magic" Experiment.
% Idea for study: Dr. Pablo R Grassi in cooperation with Prof. Dr. Andreas
% Bartels in the 'Vision and Cognition Lab' at the University of Tübingen.
% Development, implementation and execution: Vincent Plikat and Dr. Pablo R
% Grassi in cooperation with Prof. Dr. Adreas Bartels
% Started August 2019 as the Master Thesis of Vincent Plikat
% Data: functional MRI, pupil dilation, gaze position and behavioural
% rating (1-5). In this analysis only the fMRI data will be analyzed.
% Experimental design: The Experiment was divided into 3 blocks. Each block
% consisted of 4 experimental runs. In each run subjects viewed 24 videos
% (each video is considered a trial) of three different categories: Magic
% Control and Surprise. 
% The videos in each block were associated with one object (Balls, Cards
% and Sticks)
% There were 6 Magic videos with 6 coresponding control videos and 3
% surpise videos per block. 3 Different magic effects were used: Object
% appearing, vanishing and color change. The magic and surprise videos were 
% presented twice, the control videos only once. (2*6 magic + 6 control + 2*3
% surpsise = 24 videos). Half of all video presentations were flipped along
% the y-axis and the order of videos was pseudorandomized in a way that the
% same video was never presented in two consecutive trials. 
% The videos were exactly 14 seconds long and followed by a 2 second
% answering phase in which the subject was asked to rate how surprising the
% video's content was. 
% After the second run in each block the underlying methods behind each
% magic trick was presented (during this period we did not measure fMR
% data). Only when every trick was understood the following two runs were
% presented
% --> a 2(pre post revelation)*3(objects)*3(magic effect)*3(video type)
% design.
% In this analysis we analyze the contrasts of the whole video presentation
% first and the contrast of the special moment afterwards.
% In the magic videos the special moment was the moment the
% magical effect happened, in the control videos it was the corresponding 
% moment where one would expect a magic effect to happen and in the
% surprise videos it was the moment the surprising action happened. Those
% moments were selected by Vincent Plikat.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clear and close EVERYTHING
clear all;
close all;

%% Define important details of your file structure and location
% Set root directory
fprintf(['Please select your project folder.'...
    '(ideally it should contain a folder named "raw_data")'])
root_dir    = uigetdir(homedir, 'Select Project Folder');
if root_dir == 0
    error('No folder was selected --> I terminate the script')
end

% Set source_data directory. This is only needed for slice time correction
source_dir  = fullfile(root_dir,'source_data');
if ~isfolder(source_dir)
    fprintf(['It appears you do not have a "source_data" folder.\n'...
        'Please select the folder that contains your unprocessed niftis.'])
    source_dir  = uigetdir(root_dir, 'Select DICOM folder');
    if source_dir == 0
        error('No folder was selected --> I terminate the script')
    end
end

% Set raw_data directory.
raw_dir     = fullfile(root_dir, 'raw_data');
if ~isfolder(raw_dir)
    fprintf(['It appears you do not have a "raw_data" folder.\n'...
        'Please select the folder that contains your unprocessed niftis.'])
    raw_dir  = uigetdir(root_dir, 'Select unprocessed nifti folder');
    if raw_dir == 0
        error('No folder was selected --> I terminate the script')
    end
end

% Set derivative directory
derived_dir = fullfile (root_dir, 'derivative_data');
if ~isfolder(derived_dir)
    fprintf(['It appears you do not have a "derivative_data" folder.\n'...
        'Please select the folder that contains your preprocessed niftis.'])
    derived_dir  = uigetdir(root_dir, 'Select preprocessed nifti folder');
    if derived_dir == 0
        error('No folder was selected --> I terminate the script')
    end
end

%% create a folder that contains the results of the second level analysis

sec_level_dir = fullfile(derived_dir,'GroupAnalysis');
if ~isfolder(derived_dir)
    mkdir(sec_level_dir)
end

%% Define what to do
do.SpecifyDesign      = 0;
do.estimate           = 0;
do.DefContrasts       = 1;

%% Settings
settings.matprefix = input (['Please specify the prefix of your participant data.\n' ...
    '(like p for participant or s for subject. It has to be unique so that only subject folders are selected):\n'],'s');

% specify format for folder numeration
formatSpec = '%04i';

folders = dir(fullfile(derived_dir,[settings.matprefix, '*']));
subNames = {folders(:).name}; 

% load in one SPM.mat file to read out the contrasts
load(fullfile(derived_dir,subNames{1}, 'MRI/analysis/first_level_analysis/WholeVideo/SPM.mat'));
nContrasts = length(SPM.xCon);
% Iterate over all contrasts
for C = 1:nContrasts 

    if do.SpecifyDesign
    
        contrasts_dirs = {};

        for s = 1:length(subNames)

            current_dir = fullfile(derived_dir,subNames{s},'MRI/analysis/first_level_analysis/WholeVideo');
            contrasts_dirs{end+1} = fullfile(current_dir,['con_' num2str(C,formatSpec) '.nii']);

        end

        %% specify analysis parameter
        matlabbatch{1}.spm.stats.factorial_design.dir = {fullfile(sec_level_dir,'WholeVideo',SPM.xCon(C).name)};
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = contrasts_dirs';
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;


        jobs = matlabbatch;
        spm('defaults', 'FMRI');
        spm_jobman('run', jobs);


        %% respecify the contrasts to include in the analysis
        contrasts_dirs = {};

        for s = 1:length(subNames)

            current_dir = fullfile(derived_dir,subNames{s},'MRI/analysis/first_level_analysis/SpecialMoment');
            contrasts_dirs{end+1} = fullfile(current_dir,['con_' num2str(C,formatSpec) '.nii']);

        end

        %% respecify analysis parameter
        matlabbatch{1}.spm.stats.factorial_design.dir           = {fullfile(sec_level_dir,'SpecialMoment',SPM.xCon(C).name)};
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans  = contrasts_dirs';

        jobs = matlabbatch;
        spm('defaults', 'FMRI');
        spm_jobman('run', jobs);
        clear jobs; clear matlabbatch;
    end % estimate Model
    
    if do.estimate
        
        matlabbatch{1}.spm.stats.fmri_est.spmmat           = {fullfile(sec_level_dir,'WholeVideo',SPM.xCon(C).name, 'SPM.mat')};
        matlabbatch{1}.spm.stats.fmri_est.write_residuals  = 0;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        
        jobs = matlabbatch;
        spm('defaults', 'FMRI');
        spm_jobman('run', jobs);        
        clear jobs; clear matlabbatch;
        
        matlabbatch{1}.spm.stats.fmri_est.spmmat           = {fullfile(sec_level_dir,'SpecialMoment',SPM.xCon(C).name, 'SPM.mat')};
        matlabbatch{1}.spm.stats.fmri_est.write_residuals  = 0;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        
        jobs = matlabbatch;
        spm('defaults', 'FMRI');
        spm_jobman('run', jobs);        
        clear jobs; clear matlabbatch;
    end
    
    if do.DefContrasts
        
        %% Define contrasts:
        matlabbatch{1}.spm.stats.con.spmmat = {fullfile(sec_level_dir,'WholeVideo',SPM.xCon(C).name, 'SPM.mat')};
        matlabbatch{1}.spm.stats.con.delete = 1;
        
        % Name:
        ContrastName    = 'mean';
        Contrast        = 1;
        
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name       = ContrastName;
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep    = 'none'; %'replsc' if repeat for sessions
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec     = Contrast; % or weights?
        
        spm_jobman('run', matlabbatch);
        
        matlabbatch{1}.spm.stats.con.spmmat = {fullfile(sec_level_dir,'SpecialMoment',SPM.xCon(C).name, 'SPM.mat')};
        spm_jobman('run', matlabbatch);
        
        clear matlabbatch;
    end
end
clear all;
close all;