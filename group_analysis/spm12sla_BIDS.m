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
    '(ideally it should contain a folder named "derivatives")\n\n'])
rootDir    = uigetdir(homedir, 'Select Project Folder');
if rootDir == 0
    error('No folder was selected --> I terminate the script')
end

% Set sourcedata directory. This is only needed for slice time correction
sourceDir  = fullfile(rootDir,'sourcedata');
if ~isfolder(sourceDir)
    fprintf(['It appears you do not have a "sourcedata" folder.\n'...
        'Please select the folder that contains your DICOMS.'])
    sourceDir  = uigetdir(rootDir, 'Select DICOM folder');
    if sourceDir == 0
        error('No folder was selected --> I terminate the script')
    end
end

% Set rawdata directory.
rawDir     = fullfile(rootDir, 'rawdata');
if ~isfolder(rawDir)
    fprintf(['It appears you do not have a "rawdata" folder.\n'...
        'Please select the folder that contains your unprocessed niftis.'])
    rawDir  = uigetdir(rootDir, 'Select unprocessed nifti folder');
    if rawDir == 0
        error('No folder was selected --> I terminate the script')
    end
end

derivesDir = fullfile (rootDir, 'derivatives');
if ~isfolder(derivesDir)
    fprintf(['It appears you do not have a "derivatives" folder.\n'...
        'Please select the folder that contains your preprocessed niftis.'])
    derivesDir  = uigetdir(rootDir, 'Select preprocessed nifti folder');
    if derivesDir == 0
        error('No folder was selected --> I terminate the script')
    end
end

%% Data locations 
softwareName        = 'spm12';  
% get the data from the first level analysis pipeline you want
pipelineName        = 'spm12-fla';
brainMask           = 'WholeBrain';         % whole brain or ROI
conditionsAnalyzed  = 'MagicEffects';
smoothKernelSize	= 9;   % in mm
smoothKernelSpace   = 'mni';
% combine above specifications for a well structured file hierarchy
smoothnessDir       = [num2str(smoothKernelSize) 'mm-smoothed-' smoothKernelSpace 'space'];
flaDir              = fullfile(derivesDir, softwareName, pipelineName, brainMask, conditionsAnalyzed, smoothnessDir);
wholeVideoGLMName   = 'WholeVideo';
specialMomentGLMName= 'SpecialMoment';
% specify the name of the processing pipeline
analysisPipeline    = 'spm12-sla';

%% create a folder that contains the results of the second level analysis

secLevelDir = fullfile(derivesDir, softwareName, analysisPipeline, brainMask, conditionsAnalyzed, smoothnessDir);
if ~isfolder(secLevelDir)
    mkdir(secLevelDir)
end

%% Define what to do
do.SpecifyDesign      = 1;
do.estimate           = 1;
do.DefContrasts       = 1;
% Which model to do
do.wholeVideo       = 1;
do.specialMoment    = 1;

%% Settings
% specify format for folder numeration
formatSpec = '%04i';

folders = dir(fullfile(flaDir, wholeVideoGLMName,'sub-*'));
subNames = {folders(:).name}; 

% load in one SPM.mat file to read out the contrasts
load(fullfile(flaDir, wholeVideoGLMName,subNames{1},'SPM.mat'));
nContrasts = length(SPM.xCon);

%% The actual second level analysis
% Iterate over all contrasts
for C = 1:nContrasts 

    if do.SpecifyDesign
        %% specify general analysis parameter
        matlabbatch{1}.spm.stats.factorial_design.cov                       = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.multi_cov                 = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none        = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im                = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em                = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit            = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no    = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm           = 1;
        
        %% whole video analysis
        if do.wholeVideo
            dataDir         = fullfile(flaDir, wholeVideoGLMName);
            contrastsDirs   = {};

            for s = 1:length(subNames)

                currentDir              = fullfile(dataDir,subNames{s});
                contrastsDirs{end+1}    = fullfile(currentDir,['con_' num2str(C,formatSpec) '.nii']);

            end

            matlabbatch{1}.spm.stats.factorial_design.dir                       = {fullfile(secLevelDir,wholeVideoGLMName,SPM.xCon(C).name)};
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans              = contrastsDirs';

            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
        end

        %% special moment analysis
        if do.specialMoment
            dataDir         = fullfile(flaDir, specialMomentGLMName);
            contrastsDirs   = {};

            for s = 1:length(subNames)

                currentDir = fullfile(dataDir,subNames{s});
                contrastsDirs{end+1} = fullfile(currentDir,['con_' num2str(C,formatSpec) '.nii']);

            end

            matlabbatch{1}.spm.stats.factorial_design.dir           = {fullfile(secLevelDir,specialMomentGLMName,SPM.xCon(C).name)};
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans  = contrastsDirs';

            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
        end
        clear matlabbatch;
    end % estimate Model
    
    if do.estimate
        
        matlabbatch{1}.spm.stats.fmri_est.write_residuals  = 0;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        
        %% whole video analysis
        if do.wholeVideo
            matlabbatch{1}.spm.stats.fmri_est.spmmat           = {fullfile(secLevelDir,wholeVideoGLMName,SPM.xCon(C).name, 'SPM.mat')};

            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch); 
        end
        
        %% special moment analysis
        if do.specialMoment
            matlabbatch{1}.spm.stats.fmri_est.spmmat           = {fullfile(secLevelDir,specialMomentGLMName,SPM.xCon(C).name, 'SPM.mat')};

            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);        
        end
        clear matlabbatch;
    end
    
    if do.DefContrasts
        %% Define contrasts:
        matlabbatch{1}.spm.stats.con.delete = 1;
        
        % Name:
        ContrastName    = 'mean';
        Contrast        = 1;
        
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.name       = ContrastName;
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep    = 'none'; %'replsc' if repeat for sessions
        matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec     = Contrast; % or weights?
        
        %% whole video analysis
        if do.wholeVideo
            matlabbatch{1}.spm.stats.con.spmmat = {fullfile(secLevelDir,wholeVideoGLMName,SPM.xCon(C).name, 'SPM.mat')};
            spm_jobman('run', matlabbatch);
        end
        
        %% special moment analysis
        if do.specialMoment
            matlabbatch{1}.spm.stats.con.spmmat = {fullfile(secLevelDir,specialMomentGLMName,SPM.xCon(C).name, 'SPM.mat')};
            spm_jobman('run', matlabbatch);
        end
        
        clear matlabbatch;
    end
end
clear all;
close all;