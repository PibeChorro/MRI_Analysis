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
conditionsAnalyzed  = 'MagicEffects';         % Magic, Control and Surprise videos are one regressor each (Response as well, but that is less important)
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
do.DefContrasts     = 1;
% Which model to do
do.wholeVideo       = 1;
do.specialMoment    = 1;


%% Get experimental design parameters
fla.conditionNames  = {
    'Appear*Magic'  ; ...
    'Vanish*Magic'  ; ...
    'Change*Magic'  ; ...
    'Appear*Control'; ...
    'Vanish*Control'; ...
    'Change*Control'; ...
    'Surprise'      ; ...
    'Response'        ...
    };
effectVersion = {
    '1'; '2'...
    };

% how many conditions we have (in total and per "type")
fla.numConditions   = length(fla.conditionNames);
numMag              = 3;
numCon              = 3;
numSur              = 1;
numRaPara           = 6; % default number of realignment parameter, in model estimation overwritten but still 6
% how many blocks we had
numBlocks           = 3;
% should the movement be used as regressors of no interest
fla.realignmentParametersFlag  = 1;

for s = 1:length(subNames)
    %% Define where to look for functional MRI data and the logs that contain information about stimulus on/offsets
    smoothedDataDir     = fullfile(dataDir,         subNames{s}, 'func');
    realignedDataDir    = fullfile(realignedDir,    subNames{s}, 'func');
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
                    % We iterate over our conditions specified above
                    fla.trialOnset={}; fla.trialDuration={};
                    for reg = 1:length (fla.conditionNames) - 2 % The last two conditions are the surprise and response. We set them manually at the end
                        % from the matlab log file select the condition names.
                        % Look which index contains the current condition.
                        % Those indices are used to select the trial onset from
                        % the matlab log file. The very first EXPERIMENTAL
                        % trigger is substracted from Video start, to set
                        % the measurement onset = 0
                        tmpOnset    = []; % An empty array to contain specialmoment onsets all videos belonging to one condition
                        tmpDuration = [];
                        for magEff = 1:length(effectVersion)
                            indexCondition  = regexp(log.data.Condition, regexptranslate('wildcard',fla.conditionNames(reg)));  % gives a cell array containing index of wildcard start
                            indexCondition  = find(not(cellfun('isempty',indexCondition)));                                     % gives an array containing the non empty entries in indexCondition

                            indexEffectVersion          = contains(log.data.Condition(indexCondition),effectVersion(magEff));           % gives us the indeces of indexCondition containing the current magic effect version
                            indicesToUse                = indexCondition(indexEffectVersion);
                            tmpOnset = [tmpOnset; log.data.VideoStart(indicesToUse) - log.Keys.trig(end)];
                            tmpDuration = [tmpDuration; log.data.Rating_stimOn(indicesToUse) - log.data.VideoStart(indicesToUse)];
                        end
                        fla.trialOnset{end+1}       = tmpOnset;
                        fla.trialDuration{end+1}    = tmpDuration;
                    end
                    %% Surprise
                    % Do the same as in the for loop above, but for the
                    % three surprise videos
                    reg = 3;
                    tmpOnset = [];
                    tmpDuration = [];
                    for surVer = 1:3
                        indexCondition  = regexp(log.data.Condition, regexptranslate('wildcard',fla.conditionNames(end-1)));    % gives a cell array containing index of wildcard start
                        indexCondition  = find(not(cellfun('isempty',indexCondition)));                                         % gives an array containing the non empty entries in indexCondition

                        indexEffectVersion          = contains(log.data.Condition(indexCondition),effectVersion(magEff));           % gives us the indeces of indexCondition containing the current magic effect version
                        indicesToUse                = indexCondition(indexEffectVersion);
                        
                        % The very first EXPERIMENTAL trigger is 
                        % substracted from Video start, to set the
                        % measurement onset = 0
                        tmpOnset    = [tmpOnset; log.data.VideoStart(indicesToUse) - log.Keys.trig(end)];
                        tmpDuration = [tmpDuration; log.data.Rating_stimOn(indicesToUse) - log.data.VideoStart(indicesToUse)];
                    end

                    fla.trialOnset{end+1}       = tmpOnset;
                    fla.trialDuration{end+1}    = tmpDuration;
                    
                    %% Response
                    % We just add 2 seconds for the response in the last entry of trialDuration and the end of each Video as trialOnset
                    fla.trialOnset{end+1}         = log.data.Rating_stimOn - log.Keys.trig(end);
                    fla.trialDuration{end+1}      = ones(1, length (fla.trialOnset{end}))*2;
                else % don't get the matfile log files.

                    %PLACE HOLDER

                end

                % Condition names:
                for cc = 1:fla.numConditions % for CurrentCondition
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

            matlabbatch{1}.spm.stats.fmri_spec.dir = {betaLoc};   % The directory, the SPM.mat and all betas are written

            for r = 1:numRuns % For each run

                if do.loadlog == 1

                    % IMPORTANT: this throws an error if there is a space in
                    % your path. 
                    load(strtrim(logMatfiles(r,:)));

                    % specify the trial onsets and durations.
                    % We iterate over Magic and Control conditions 
                    % specified above
                    fla.trialOnset={}; fla.trialDuration={};
                    %% Magic and Control
                    for reg = 1:length (fla.conditionNames) - 2 % the last two conditions are Surprise and Rating
                        tmp = []; % An empty array to contain specialmoment onsets all videos belonging to one condition
                        for magEff = 1:length(effectVersion)
                            % Get the indices for for every video belonging
                            % to the current condition
                            indexCondition  = regexp(log.data.Condition, regexptranslate('wildcard',fla.conditionNames(reg)));  % gives a cell array containing index of wildcard start
                            indexCondition  = find(not(cellfun('isempty',indexCondition)));                                     % gives an array containing the non empty entries in indexCondition

                            indexEffectVersion          = contains(log.data.Condition(indexCondition),effectVersion(magEff));           % gives us the indeces of indexCondition containing the current magic effect version
                            indicesToUse                = indexCondition(indexEffectVersion);
                            
                            VideoNames = log.data.Condition(indicesToUse);
                            % The videoname may contain an '_F' for the flip
                            % condition. Therefore we need to get the normal video
                            % name
                            if contains(VideoNames{1},'_F')
                                VideoName   =VideoNames{1};
                                VideoName   =VideoName(1:end-2);
                            else
                                VideoName   =VideoNames{1};
                            end
                            % Select the frame of the special moment from the
                            % 'do' struct
                            SpecialMoment               = do.all_frames_of_effect(contains(do.ListOfVideos,VideoName));
                            SpecialMomentOnset          = SpecialMoment{1}*frameTime; % multiply with the frameduration

                            % The very first EXPERIMENTAL trigger is 
                            % substracted from Video start, to set the
                            % measurement onset = 0
                            tmp = [tmp; log.data.VideoStart(indicesToUse) - log.Keys.trig(end) + SpecialMomentOnset];
                        end
                       
                        fla.trialOnset{end+1}       = tmp;
                        fla.trialDuration{end+1}    = zeros(1,length(tmp));
                    end
                    
                    %% Surprise
                    % Do the same as in the for loop above, but for the
                    % three surprise videos
                    reg = 3;
                    tmp = [];
                    
                    for surVer = 1:3
                        VideoNames = log.data.Condition(contains(log.data.Condition,fla.conditionNames(end-1)) & contains(log.data.Condition, num2str(surVer)));
                        % The videoname may contain an '_F' for the flip
                        % condition. Therefore we need to get the normal video
                        % name
                        if contains(VideoNames{1},'_F')
                            VideoName   =VideoNames{1};
                            VideoName   =VideoName(1:end-2);
                        else
                            VideoName   =VideoNames{1};
                        end
                        % Select the frame of the special moment from the
                        % 'do' struct
                        SpecialMoment               = do.all_frames_of_effect(contains(do.ListOfVideos,VideoName));
                        SpecialMomentOnset          = SpecialMoment{1}*frameTime; % multiply with the frameduration
                        
                        % The very first EXPERIMENTAL trigger is 
                        % substracted from Video start, to set the
                        % measurement onset = 0
                        tmp = [tmp; log.data.VideoStart(contains(log.data.Condition,VideoName))-log.Keys.trig(end) + SpecialMomentOnset];
                    end

                    fla.trialOnset{end+1}       = tmp;
                    fla.trialDuration{end+1}    = zeros(1,length(tmp));
                    %% Response
                    % We just add 2 seconds for the response in the last entry of trialDuration and the end of each Video as trialOnset
                    fla.trialOnset{end+1}       = log.data.Rating_stimOn - log.Keys.trig(end);
                    fla.trialDuration{end+1}    = ones(1, length (fla.trialOnset{end}))*2;
                else % don't get the matfile log files.

                  %PLACE HOLDER

                end % End load mat files

                % Condition names:
                for cc = 1:fla.numConditions % for CurrentCondition
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
    
    
    %% Define contrast
    if do.DefContrasts == 1
        % General settings
        % Contrast Names:
        ContrastNames = ...
            {'Magic > Control';                     ... %1
            'Control > Magic';                      ... %2
            'Magic > Control Before';               ... %3
            'Control > Magic Before';               ... %4
            'Magic Before > Magic After';           ... %5
            'Magic After > Magic Before';           ... %6
            'Magic > Surprise Before';              ... %7
            'Surpise > Magic Before';               ... %8
            'Surprise > Control';                   ... %9
            'Control > Surprise';                   ... %10
            'Magic > Control After';                ... %11
            'Control > Magic After';                ... %12
            'Magic > Surprise After';               ... %13
            'Surpise > Magic After';                ... %14
            'MagPre-ConPre vs MagPost-ConPost';     ... %15
            'MagPost-ConPost vs MagPre-ConPre';     ... %16
            'Appear Before > Appear After';         ... %17
            'Vanish Before > Vanish After';         ... %18
            'Change Before > Change After';         ... %19
            'Appear After > Appear Before';         ... %20
            'Vanish After > Vanish Before';         ... %21
            'Change After > Change Before';         ... %22
            'Appear > Control';                     ... %23
            'Control > Appear';                     ... %24
            'Vanish > Control';                     ... %25
            'Control > Vanish';                     ... %26
            'Change > Control';                     ... %27
            'Control > Change';                     ... %28
            'Appear > Control Before';              ... %29
            'Control > Appear Before';              ... %30
            'Vanish > Control Before';              ... %31
            'Control > Vanish Before';              ... %32
            'Change > Control Before';              ... %33
            'Control > Change Before';              ... %34
            'Appear > Surprise Before';             ... %35
            'Vanish > Surprise Before';             ... %36
            'Change > Surprise Before';             ... %37
            'Appear > Control After';               ... %38
            'Vanish > Control After';               ... %39
            'Change > Control After';               ... %40
            'Appear > Surprise After';              ... %41
            'Vanish > Surprise After';              ... %42
            'Change > Surprise After';              ... %43
            % Interaction effects
            'AppPre-ConPre vs AppPost-ConPost';     ... %44
            'VanPre-ConPre vs Vanpost-ConPost';     ... %45
            'ChaPre-ConPre vs ChaPost-ConPost';     ... %46
            'AppPost-ConPost vs AppPre-ConPre';     ... %47
            'Vanpost-ConPost vs Vanpost-ConPost';   ... %48
            'ChaPost-ConPost vs ChaPre-ConPre';     ... %49
            % Contrats to outrule the timeconfound by comparing run 1vs2
            % and run 2vs3 - the same time difference, but the first is pre
            % vs pre and the other is pre vs post revelation
            'Magic PreVsPre (run 1vs2)';            ... %50
            'Magic PreVsPost (run 2vs3)';           ... %51
            'Video vs Response';                    ... %52
            'Response vs Video'                     ... %53
            };
        
        % Contrast values
        
        %       PreRevelation       Magic videos            Control             Surprise    Response    Realignment      PostRevelation  Magic videos           Control             Surprise           Response   Realignment
        C1 = repmat([repmat([ones(1,numMag)         ones(1,numCon)*(-1) zeros(1,numSur)     0           zeros(1,numRaPara)],1,2) repmat([ones(1,numMag)         ones(1,numCon)*(-1) zeros(1,numSur)     0          zeros(1,numRaPara)],1,2)],1,numBlocks);
        C2 = repmat([repmat([ones(1,numMag)*(-1)    ones(1,numCon)      zeros(1,numSur)     0           zeros(1,numRaPara)],1,2) repmat([ones(1,numMag)*(-1)    ones(1,numCon)      zeros(1,numSur)     0          zeros(1,numRaPara)],1,2)],1,numBlocks);
        C3 = repmat([repmat([ones(1,numMag)         ones(1,numCon)*(-1) zeros(1,numSur)     0           zeros(1,numRaPara)],1,2) repmat([zeros(1,numMag)        zeros(1,numCon)     zeros(1,numSur)     0          zeros(1,numRaPara)],1,2)],1,numBlocks);
        C4 = repmat([repmat([ones(1,numMag)*(-1)    ones(1,numCon)      zeros(1,numSur)     0           zeros(1,numRaPara)],1,2) repmat([zeros(1,numMag)        zeros(1,numCon)     zeros(1,numSur)     0          zeros(1,numRaPara)],1,2)],1,numBlocks);
        C5 = repmat([repmat([ones(1,numMag)         zeros(1,numCon)     zeros(1,numSur)     0           zeros(1,numRaPara)],1,2) repmat([ones(1,numMag)*(-1)    zeros(1,numCon)     zeros(1,numSur)     0          zeros(1,numRaPara)],1,2)],1,numBlocks);
        C6 = repmat([repmat([ones(1,numMag)*(-1)    zeros(1,numCon)     zeros(1,numSur)     0           zeros(1,numRaPara)],1,2) repmat([ones(1,numMag)         zeros(1,numCon)     zeros(1,numSur)     0          zeros(1,numRaPara)],1,2)],1,numBlocks);
        C7 = repmat([repmat([ones(1,numMag)         zeros(1,numCon)     ones(1,numSur)*(-3) 0           zeros(1,numRaPara)],1,2) repmat([zeros(1,numMag)        zeros(1,numCon)     zeros(1,numSur)     0          zeros(1,numRaPara)],1,2)],1,numBlocks);
        C8 = repmat([repmat([ones(1,numMag)*(-1)    zeros(1,numCon)     ones(1,numSur)*3    0           zeros(1,numRaPara)],1,2) repmat([zeros(1,numMag)        zeros(1,numCon)     zeros(1,numSur)     0          zeros(1,numRaPara)],1,2)],1,numBlocks);
        C9 = repmat([repmat([zeros(1,numMag)        ones(1,numCon)*(-1) ones(1,numSur)*3    0           zeros(1,numRaPara)],1,2) repmat([zeros(1,numMag)        ones(1,numCon)*(-1) ones(1,numSur)*3    0          zeros(1,numRaPara)],1,2)],1,numBlocks);
        C10= repmat([repmat([zeros(1,numMag)        ones(1,numCon)      ones(1,numSur)*(-3) 0           zeros(1,numRaPara)],1,2) repmat([zeros(1,numMag)        ones(1,numCon)      ones(1,numSur)*(-3) 0          zeros(1,numRaPara)],1,2)],1,numBlocks);
        C11= repmat([repmat([zeros(1,numMag)        zeros(1,numCon)     zeros(1,numSur)     0           zeros(1,numRaPara)],1,2) repmat([ones(1,numMag)         ones(1,numCon)*(-1) zeros(1,numSur)     0          zeros(1,numRaPara)],1,2)],1,numBlocks);
        C12= repmat([repmat([zeros(1,numMag)        zeros(1,numCon)     zeros(1,numSur)     0           zeros(1,numRaPara)],1,2) repmat([ones(1,numMag)*(-1)    ones(1,numCon)      zeros(1,numSur)     0          zeros(1,numRaPara)],1,2)],1,numBlocks);
        C13= repmat([repmat([zeros(1,numMag)        zeros(1,numCon)     zeros(1,numSur)     0           zeros(1,numRaPara)],1,2) repmat([ones(1,numMag)         zeros(1,numCon)     ones(1,numSur)*(-3) 0          zeros(1,numRaPara)],1,2)],1,numBlocks);
        C14= repmat([repmat([zeros(1,numMag)        zeros(1,numCon)     zeros(1,numSur)     0           zeros(1,numRaPara)],1,2) repmat([ones(1,numMag)*(-1)    zeros(1,numCon)     ones(1,numSur)*3    0          zeros(1,numRaPara)],1,2)],1,numBlocks);
        %       Interaction Effects
        C15= repmat([repmat([ones(1,numMag)         ones(1,numCon)*(-1) zeros(1,numSur)     0           zeros(1,numRaPara)],1,2) repmat([ones(1,numMag)*(-1) ones(1,numCon)      zeros(1,numSur)     0          zeros(1,numRaPara)],1,2)],1,numBlocks);
        C16= repmat([repmat([ones(1,numMag)*(-1)    ones(1,numCon)      zeros(1,numSur)     0           zeros(1,numRaPara)],1,2) repmat([ones(1,numMag)      ones(1,numCon)*(-1) zeros(1,numSur)     0          zeros(1,numRaPara)],1,2)],1,numBlocks);
        %                   AppearM VanishM ChangeM AppearC VanishC ChangeC Surprise            Response    Realignment     PostRevelation  AppearM VanishM ChangeM AppearC VanishC ChangeC Surprise            Response    Realignment
        C17= repmat([repmat([1      0       0       0       0       0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       -1      0       0       0       0       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C18= repmat([repmat([0      1       0       0       0       0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       -1      0       0       0       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C19= repmat([repmat([0      0       1       0       0       0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       0       -1      0       0       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C20= repmat([repmat([-1     0       0       0       0       0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       1       0       0       0       0       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C21= repmat([repmat([0      -1      0       0       0       0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       1       0       0       0       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C22= repmat([repmat([0      0       -1      0       0       0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       0       1       0       0       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C23= repmat([repmat([1      0       0       -1      0       0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       1       0       0       -1      0       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C24= repmat([repmat([-1     0       0       1       0       0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       -1      0       0       1       0       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C25= repmat([repmat([0      1       0       0       -1      0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       1       0       0       -1      0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C26= repmat([repmat([0      -1      0       0       1       0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       -1      0       0       1       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C27= repmat([repmat([0      0       1       0       0       -1      zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       0       1       0       0       -1      zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C28= repmat([repmat([0      0       -1      0       0       1       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       0       -1      0       0       1       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C29= repmat([repmat([1      0       0       -1      0       0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       0       0       0       0       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C30= repmat([repmat([-1     0       0       1       0       0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       0       0       0       0       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C31= repmat([repmat([0      1       0       0       -1      0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       0       0       0       0       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C32= repmat([repmat([0      -1      0       0       1       0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       0       0       0       0       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C33= repmat([repmat([0      0       1       0       0       -1      zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       0       0       0       0       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C34= repmat([repmat([0      0       -1      0       0       1       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       0       0       0       0       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C35= repmat([repmat([1      0       0       0       0       0       ones(1,numSur)*(-1) 0   zeros(1,numRaPara)],1,2) repmat([       0       0       0       0       0       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C36= repmat([repmat([0      1       0       0       0       0       ones(1,numSur)*(-1) 0   zeros(1,numRaPara)],1,2) repmat([       0       0       0       0       0       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C37= repmat([repmat([0      0       1       0       0       0       ones(1,numSur)*(-1) 0   zeros(1,numRaPara)],1,2) repmat([       0       0       0       0       0       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C38= repmat([repmat([0      0       0       0       0       0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       1       0       0       -1      0       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C39= repmat([repmat([0      0       0       0       0       0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       1       0       0       -1      0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C40= repmat([repmat([0      0       0       0       0       0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       0       1       0       0       -1      zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C41= repmat([repmat([0      0       0       0       0       0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       1       0       0       0       0       0       ones(1,numSur)*(-1)  0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C42= repmat([repmat([0      0       0       0       0       0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       1       0       0       0       0       ones(1,numSur)*(-1)  0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C43= repmat([repmat([0      0       0       0       0       0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       0       1       0       0       0       ones(1,numSur)*(-1)  0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        %       Interaction Effects
        C44= repmat([repmat([1      0       0       -1      0       0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       -1      0       0       1       0       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C45= repmat([repmat([0      1       0       0       -1      0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       -1      0       0       1       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C46= repmat([repmat([0      0       1       0       0       -1      zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       0       -1      0       0       1       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C47= repmat([repmat([-1     0       0       1       0       0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       1       0       0       -1      0       0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C48= repmat([repmat([0      -1      0       0       1       0       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       1       0       0       -1      0       zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        C49= repmat([repmat([0      0       -1      0       0       1       zeros(1,numSur)     0   zeros(1,numRaPara)],1,2) repmat([       0       0       1       0       0       -1      zeros(1,numSur)      0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        %   FirstRun Magic           Control         Surprise        Response   Realignment SecondRun   Magic                Control         Surprise        Response   Realignment Postrevelation  Magic           Control         Surprise        Response    Realignment
        C50= repmat([ones(1,numMag)  zeros(1,numCon) zeros(1,numSur) 0          zeros(1,numRaPara)      ones(1,numMag)*(-1)  zeros(1,numCon) zeros(1,numSur) 0          zeros(1,numRaPara)  repmat([zeros(1,numMag) zeros(1,numCon) zeros(1,numSur) 0           zeros(1,numRaPara)],1,2)],1,numBlocks);
        %   FirstRun Magic           Control         Surprise        Response   Realignment SecondRun   Magic            Control         Surprise        Response   Realignment   ThirdRun  Magic                 Control         Surprise        Response  Realignment FourthRun   Magic           Control         Surprise        Response Realignment
        C51= repmat([zeros(1,numMag) zeros(1,numCon) zeros(1,numSur) 0          zeros(1,numRaPara)      ones(1,numMag)   zeros(1,numCon) zeros(1,numSur) 0          zeros(1,numRaPara)      ones(1,numMag)*(-1)   zeros(1,numCon) zeros(1,numSur) 0         zeros(1,numRaPara)      zeros(1,numMag) zeros(1,numCon) zeros(1,numSur) 0 zeros(1,numRaPara)],1,numBlocks);
        %                  Magic                No Magic            Surprise            Response    Realigment
        C52= repmat(repmat([ones(1,numMag)      ones(1,numCon)      ones(1,numSur)      -7         zeros(1,numRaPara)],1,4),1,numBlocks);
        C53= repmat(repmat([ones(1,numMag)*(-1) ones(1,numCon)*(-1) ones(1,numSur)*(-1) 7          zeros(1,numRaPara)],1,4),1,numBlocks);
        
        % Combine all Contrasts in one Matrix
        Contrasts = [C1; C2; C3; C4; C5; C6; C7; C8; C9; C10; C11; C12; C13; C14;... % All Magic and Control videos combined
            C15; C16; C17; C18; C19; C20; C21; C22; C23; C24; C25; C26; C27; C28; C29; C30; C31; C32; C33; C34; C35; C36; C37; C38; C39; C40; C41; C42; C43; C44; C45; C46; C47; C48; C49; ... % Magic effects seperated
            C50; C51; C52; C53]; % Control contrasts
        
        % safety net: check if sum of contrasts is 0
        if any(sum(Contrasts,2))
            error('One of the contrasts sum is not zero')
        end
        
        
        %% Define contrasts:
        
        matlabbatch{1}.spm.stats.con.delete = 1; % Delete existing Contrasts
        
        % Write in spm-struct:
        for C = 1:size(Contrasts, 1)
            matlabbatch{1}.spm.stats.con.consess{C}.tcon.name       = ContrastNames{C,:};
            matlabbatch{1}.spm.stats.con.consess{C}.tcon.sessrep    = 'none'; %'replsc' if repeat for sessions
            matlabbatch{1}.spm.stats.con.consess{C}.tcon.convec     = Contrasts(C,:); % or weights?
        end
        
        if do.wholeVideo
            % Set directory for whole Video analysis
            analysisName    = 'WholeVideo';
            betaLoc         = fullfile(destDir, analysisName, subNames{s});
            matlabbatch{1}.spm.stats.con.spmmat = {fullfile(betaLoc, 'SPM.mat')};
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
        end
        
        %% Contrasts for analysis of the special moment. 
        
        if do.specialMoment
            % Set directory for special moment analysis
            analysisName    = 'SpecialMoment';
            betaLoc         = fullfile(destDir, analysisName, subNames{s});
            matlabbatch{1}.spm.stats.con.spmmat = {fullfile(betaLoc, 'SPM.mat')};
            spm('defaults', 'FMRI');
            spm_jobman('run', matlabbatch);
        end
        
        clear matlabbatch;
        
    end % define contrast
    
    toc;
    
end

clear all