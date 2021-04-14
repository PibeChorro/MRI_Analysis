%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for the first level analysis of the "Magic" Experiment.
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
% In this analysis we analyze the data of the whole video presentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; % first clear all variables, so we don't have any intervening variables.
tic; % start script.

%% Define important details of your file structure and location
% Set root directory
fprintf(['Please select your project folder.'...
    '(ideally it should contain a folder named "rawdata")\n\n'])
rootDir    = '/Users/vpl/Documents/Master_Thesis/DATA/MRI'; % uigetdir(homedir, 'Select Project Folder');
if rootDir == 0
    error('No folder was selected --> I terminate the script')
end

% Set sourcedata directory. This is needed to get DICOM header information
% and the names of subjects in the PsychoPhysics data folder
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
softwareName        = 'spm12';              % software used to create preprocessed data

% specify the name of the analysis pipeline
analysisPipeline    = 'spm12-fla';          % how is the folder named that contains first level results
brainMask           = 'WholeBrain';         % whole brain or ROI
conditionsAnalyzed  = 'EveryVideo';      % Every magic effect and every version of an effect is a regressor (Appear1, Appear2, Change1, etc.)
smoothKernelSize	= 9;                    % in mm
smoothKernelSpace   = 'mni';                % mni or native (mni makes more sense, native rather for explorative analysis... maybe)
% combine above specifications for a well structured file hierarchy
smoothnessDir       = [num2str(smoothKernelSize) 'mm-smoothed-' smoothKernelSpace 'space'];                     % Name of smoothed data directory
destDir             = fullfile (derivesDir, softwareName, analysisPipeline, brainMask, conditionsAnalyzed, smoothnessDir);  % where all the results of the TWO GLMs are stored
% Get name, location and number of sourcedata subjects

DICOMprefix         = 'sMag'; % input (['Please specify the prefix of your participant data in your SOURCE DATA.\n' ...
%'(like p for participant or s for subject. It has to be unique so that only subject folders are selected):\n'],'s');
DICOMsubNames       = spm_select('List', sourceDir, 'dir', ['^' DICOMprefix]);  % ^ is needed so that the filter only searches for folders STARTING with DICOMprefix
DICOMsubNames       = cellstr(DICOMsubNames);                                   % cellstring format is needed for spm

% Get name, location and number of preprocessed subjects
pipelineName        = 'spm12-preproc';                                                  % how is the folder named that contains preprocessed data
dataDir             = fullfile(derivesDir, softwareName, pipelineName, smoothnessDir);  % data that is used for the analysis
realignedDir        = fullfile(derivesDir, softwareName, pipelineName, 'realigned');    % needed for realignment files as regressors of no interest
subNames            = spm_select('List', dataDir, 'dir', 'sub-');
subNames            = cellstr(subNames);
%% Multiband Factor
multibandFactor = 2;

%% Get the information about the videos
% read in the .mat file, that contains the information about the magic
% moment
fprintf('Please select the .mat file, that contains the information about the special moments.\n\n')
videoInfoMatfileDir = uigetfile(pwd, 'Select the .mat file');
if videoInfoMatfileDir == 0
        error('No .mat file was selected --> I terminate the script')
end
load(videoInfoMatfileDir);

% information about the videos. Important is only the framerate and thus
% the time every frame was presented
fps         = 25;
frameTime   = 1/fps;

%% Define what to do
do.SpecifyDesign    = 1;
do.loadlog          = 1; % load LOG files!
do.estimate         = 1;
% Which model to do
do.wholeVideo       = 1;
do.specialMoment    = 1;

% should the movement be used as regressors of no interest
fla.realignmentParametersFlag  = 1;

for s =  1:length(subNames)
    %% Define where to look for functional MRI data and the logs that contain information about stimulus on/offsets
    smoothedDataDir     = fullfile(dataDir,         subNames{s},'func');
    realignedDataDir    = fullfile(realignedDir,    subNames{s},'func');
    psyphysicDataDir    = fullfile(derivesDir, 'PsychoPhysic',DICOMsubNames{s});
    % Further information - number of runs and where a DICOM file can be
    % found
    runs                = cellstr(spm_select('List', smoothedDataDir, '.nii')); 
    numRuns             = length(runs); 
    sourcedataRuns      = spm_select('FPList',fullfile(sourceDir,DICOMsubNames{s},'func'), 'dir','^run*');
    sourcedataRuns      = cellstr(sourcedataRuns);
    %% load a dicom header that contains information needed for analysis
    dicomFiles  = spm_select('FPList', sourcedataRuns{1}, '**.IMA');    % ** is needed because one * stops at the first period. ** filters for the whole filename
    dicomFiles  = cellstr(dicomFiles);
    hdr         = spm_dicom_headers(dicomFiles{1});
    
    %% Model specification of FLA
    if do.SpecifyDesign == 1 % Model specification
        
        %% DEFINE MODEL PARAMETERS GENERAL
        matlabbatch{1}.spm.stats.fmri_spec.timing.units     = 'secs';                       % 'secs' or 'scans' unit
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT        = hdr{1}.RepetitionTime/1000;   % TR in seconds!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=spm;70bf5c7f.0904
        % suggest to set the microtime resolution to the number of slices,
        % if slice time correction was implemented and microtime onset to
        % the half. 
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t    = hdr{1}.Private_0019_100a/multibandFactor;             % Number of slices - before it was 16
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0   = round(hdr{1}.Private_0019_100a/(multibandFactor*2));	% Number of slices - before it was 8 --> T0 = reference slice is the one in the middle. 18 or 8
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        matlabbatch{1}.spm.stats.fmri_spec.fact             = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];                % derivatives
        matlabbatch{1}.spm.stats.fmri_spec.volt             = 1;                    % CHECK
        matlabbatch{1}.spm.stats.fmri_spec.global           = '';                   % Global scaling? if so: 'Scaling' if not empty: ''
        matlabbatch{1}.spm.stats.fmri_spec.mthresh          = 0.8;                  % Masking threshold
        matlabbatch{1}.spm.stats.fmri_spec.mask             = {''};                 % Mask?
        matlabbatch{1}.spm.stats.fmri_spec.cvi              = 'AR(1)';
        
        %% Get the realigment files
        AlignmentFiles = spm_select('FPList',realignedDataDir,'^rp_.*.txt');
        
        % Load .mat files that contain condition information, such as
        % condition names, onsets, offsets, etc.
        if do.loadlog == 1 % get the mat files
            logMatfiles = spm_select('FPList', psyphysicDataDir, 'log.mat'); % SELECT MAT FILES
        end
        
        
        for r = 1:numRuns % For each run
            
            % Find functional files
            dirfiles    = spm_select('ExtFPList',smoothedDataDir, runs(r), Inf);
            
            % check if any files were selected. If not stop the script
            if strcmp(dirfiles,'')
                warning('No files selected!');
                return;
            else
                nfiles     = length(dirfiles);
                fprintf('We got %s files! \n', num2str(nfiles));
            end
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans        = cellstr(dirfiles); %For every run scans are added to the matlabbatch
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi        = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg    = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf          = 128;        % High-pass filter
            
            % delete the variable 'dirfiles'
            clear dirfiles
            
            %......Include realignment parameters for each run
            if fla.realignmentParametersFlag == 1
                
                fprintf('Adding realignment parameters to design matrix! \n');
                
                % Define names for 'Motion Regressors' (aka. Realignment
                % Parameters)
                raParamNames = [
                    'rp 1';...
                    'rp 2';...
                    'rp 3';...
                    'rp 4';...
                    'rp 5';...
                    'rp 6'
                    ];
                numRaPara = length(raParamNames);
                
                % Load alignment parameters
                raParamValues    = load(AlignmentFiles(r,:)); % get realignment parameters
                
                % Add realignment parameters as regressors of no interest.
                for P = 1:length (raParamNames)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(P).name = raParamNames(P,:);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(P).val  = raParamValues(:,P);
                end
            end
        end
        
        if do.wholeVideo
            %% Name, create and set directory for this analysis
            analysisName    = 'WholeVideo';
            betaLoc         = fullfile(destDir, analysisName, subNames{s});

            if ~isfolder(betaLoc)
                mkdir(betaLoc);
            end

            matlabbatch{1}.spm.stats.fmri_spec.dir = {betaLoc};  % The directory, the SPM.mat and all betas are written

            for r = 1:numRuns % For each run

                if do.loadlog == 1

                    % IMPORTANT: this throws an error if there is a space in
                    % your path.
                    load(strtrim(logMatfiles(r,:)));

                    % specify the trial onsets and durations.
                    VideoNames          = log.data.Condition;
                    fla.conditionNames  = convertStringsToChars(VideoNames);
                    fla.nconditions     = length(fla.conditionNames);
                    fla.trialOnset     = {};
                    fla.trialDuration  = {};
                    for reg = 1:length (VideoNames)
                        % get the video name so you can extract the Timepoint
                        % of magic moment from the 'do' struct
                        % The videoname may contain an '_F' for the flip
                        % condition. Therefore we need to get the normal video
                        % name
                        currentVideo = VideoNames(reg);
                        fla.trialOnset{end+1}      = log.data.VideoStart(contains(log.data.Condition,currentVideo))-log.Keys.trig(end);

                        % maybe remove the duration
                        fla.trialDuration{end+1}   = log.data.Rating_stimOn(contains(log.data.Condition,currentVideo))...
                            -fla.trialOnset{end};
                    end
                    % We just add 2 seconds for the response in the last entry of trialDuration and the end of each Video as trialOnset
                    fla.conditionNames{end+1}       = 'Response';
                    fla.trialOnset{end+1}           = log.data.Rating_stimOn - log.Keys.trig(end);
                    fla.trialDuration{end+1}        = ones(1, length (fla.trialOnset{end}))*2;
                else % don't get the matfile log files.

                    %PLACE HOLDER

                end

                % Condition names:
                for cc = 1:length(fla.trialOnset) % for CurrentCondition
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).name        = fla.conditionNames{cc};
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).onset       = fla.trialOnset{cc};
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).duration    = fla.trialDuration{cc};

                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).pmod        = struct('name', {}, 'param', {}, 'poly', {});
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).tmod        = 0;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).orth        = 1;
                end
            end
        
        %% Specify design in SPM.mat
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
        end
        
        if do.specialMoment
            %% Name, create and set directory for this analysis
            analysisName    = 'SpecialMoment';
            betaLoc         = fullfile(destDir, analysisName, subNames{s});

            if ~isfolder(betaLoc)
                mkdir(betaLoc);
            end

            matlabbatch{1}.spm.stats.fmri_spec.dir              = {betaLoc};   % The directory, the SPM.mat and all betas are written

            for r = 1:numRuns % For each run

                if do.loadlog == 1

                    % IMPORTANT: this throws an error if there is a space in
                    % your path. 
                    load(strtrim(logMatfiles(r,:)));

                    % specify the trial onsets and durations.
                    VideoNames          = log.data.Condition;
                    fla.conditionNames  = convertStringsToChars(VideoNames);
                    fla.nconditions     = length(fla.conditionNames);
                    fla.trialOnset      = {}; 
                    fla.trialDuration   = {};
                    for reg = 1:length (VideoNames)
                        % get the video name so you can extract the Timepoint
                        % of magic moment from the 'do' struct
                        % The videoname may contain an '_F' for the flip
                        % condition. Therefore we need to get the normal video
                        % name
                        currentVideo = fla.conditionNames{reg};
                        if contains(currentVideo,'_F')
                            normalVideoName=currentVideo(1:end-2);
                        else
                            normalVideoName=currentVideo;
                        end
                        SpecialMoment               = do.all_frames_of_effect(contains(do.ListOfVideos,normalVideoName));
                        SpecialMomentOnset          = SpecialMoment{1}*frameTime;
                        fla.trialOnset{end+1}      = log.data.VideoStart(contains(log.data.Condition,currentVideo))-log.Keys.trig(end)+...
                            SpecialMomentOnset;

                        fla.trialDuration{end+1}   = 0;
                    end
                    % We just add 2 seconds for the response in the last entry of trialDuration and the end of each Video as trialOnset
                    fla.conditionNames{end+1}       = 'Response';
                    fla.trialOnset{end+1}         = log.data.Rating_stimOn - log.Keys.trig(end);
                    fla.trialDuration{end+1}      = ones(1, length (fla.trialOnset{end}))*2;
                else % don't get the matfile log files.

                  %PLACE HOLDER

                end % End load mat files

                % Condition names:
                for cc = 1:length(fla.trialOnset) % for CurrentCondition
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).name     = fla.conditionNames{cc};
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).onset    = fla.trialOnset{cc};
                    % For this analysis we do  set the duration to 0, since we
                    % are only interested in the effect of the special onset
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).duration = fla.trialDuration{cc};

                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).pmod        = struct('name', {}, 'param', {}, 'poly', {});
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).tmod        = 0;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).orth        = 1;
                end
            end
        
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
        end
        clear matlabbatch;
        
    end % model specification
    
    %% Estimate design:
    if do.estimate == 1
        % General settings
        matlabbatch{1}.spm.stats.fmri_est.write_residuals  = 0;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        
        if do.wholeVideo
            % Set directory for whole Video analysis
            analysisName    = 'WholeVideo';
            betaLoc         = fullfile(destDir, analysisName, subNames{s});

            matlabbatch{1}.spm.stats.fmri_est.spmmat           = {fullfile(betaLoc, 'SPM.mat')};
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
        end
        
        if do.specialMoment
            % Set directory for special moment analysis
            analysisName    = 'SpecialMoment';
            betaLoc         = fullfile(destDir, analysisName, subNames{s});

            matlabbatch{1}.spm.stats.fmri_est.spmmat           = {fullfile(betaLoc, 'SPM.mat')};
        
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
        end
        clear matlabbatch;
    end % end estimate
end

clear all