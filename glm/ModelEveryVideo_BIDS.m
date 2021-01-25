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
% In this analysis we analyze the data dependent on the 'special' moment in
% each video. In the magic videos the special moment was the moment the
% magical effect happened, in the control videos it was the corresponding 
% moment where one would expect a magic effect to happen and in the
% surprise videos it was the moment the surprising action happened. Those
% moments were selected by Vincent Plikat.
% Additionally we model each video by to later correlate beta maps with
% pupil dilation and behavioural ratings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear matlabbatch; % first clear unfinished batches.
tic; % start script.

%% Define important details of your file structure and location
% Set root directory
fprintf(['Please select your project folder.'...
    '(ideally it should contain a folder named "raw_data"'...
    'and a folder named "derivative_data")\n\n'])
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

derived_dir = fullfile (root_dir, 'derivative_data');
if ~isfolder(derived_dir)
    fprintf(['It appears you do not have a "derivative_data" folder.\n'...
        'Please select the folder that contains your preprocessed niftis.'])
    derived_dir  = uigetdir(root_dir, 'Select preprocessed nifti folder');
    if derived_dir == 0
        error('No folder was selected --> I terminate the script')
    end
end

%% Define what to do
do.SpecifyDesign      = 1;
do.loadlog            = 1; % load LOG files!
do.estimate           = 1;

%% Settings
settings.matprefix = input (['Please specify the prefix of your participant data.\n' ...
    '(like p for participant or s for subject. It has to be unique so that only subject folders are selected):\n'],'s');
settings.preprocessing = ['^s6wauf' '.*nii']; % ['^s2au' '.*nii']: realigned, slice time corrected, 6mm-smoothed data

data_dir = fullfile(derived_dir,'6smoothed');
folders = dir(fullfile(data_dir,[settings.matprefix, '*']));

subNames = {folders(:).name}; 
% create a "first_level_analysis" folder
spm_mkdir(fullfile(derived_dir), 'spm12flaModelEveryVideo', subNames);

fla.realignmentParameters_flag  = 1;

for s = 1:length(subNames)
    %% Define where to store and the results and where to look for functional and anatomical data
    beta_loc            = fullfile(derived_dir,'spm12flaModelEveryVideo',subNames{s});
    smoothed_data_dir   = fullfile(derived_dir,'6smoothed',subNames{s},'func'); % TODO: make more elegant
    realigned_data_dir  = fullfile(derived_dir,'realigned',subNames{s},'func'); % TODO: make more elegant
    psyphysic_data_dir  = fullfile(derived_dir,'PsychoPhysic',subNames{s});
    runs                = dir(fullfile(smoothed_data_dir,'run*'));
    nruns               = length(runs); % Number of Runs
    
    %% Model specification of FLA
    if do.SpecifyDesign == 1 % model specification
        
        %% DEFINE MODEL PARAMETERS GENERAL
        matlabbatch{1}.spm.stats.fmri_spec.dir              = {beta_loc};           % The directory, the SPM.mat and all betas are written
        matlabbatch{1}.spm.stats.fmri_spec.timing.units     = 'secs';               % 'secs' or 'scans' unit
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT        = 2;                    % TR
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t    = 16;                   % micro-timing 36 or 16
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0   = 8;                    % if 8 --> T0 = reference slice is the one in the middle. 18 or 8
        matlabbatch{1}.spm.stats.fmri_spec.fact             = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];                % derivatives
        matlabbatch{1}.spm.stats.fmri_spec.volt             = 1;                    % CHECK
        matlabbatch{1}.spm.stats.fmri_spec.global           = '';                   % Global scaling? if so: 'Scaling' if not empty: ''
        matlabbatch{1}.spm.stats.fmri_spec.mthresh          = 0.8;                  % Masking threshold
        matlabbatch{1}.spm.stats.fmri_spec.mask             = {''};                 % Mask?
        matlabbatch{1}.spm.stats.fmri_spec.cvi              = 'AR(1)';
        
        nDummies                                            = 5;
        dummiesTime                                         = nDummies*matlabbatch{1}.spm.stats.fmri_spec.timing.RT;
        
        for r = 1:nruns % For each run
            
            % Find functional files
            Path2Run    = fullfile(smoothed_data_dir,runs(r).name); 
            dirfiles    = spm_select('FPList',Path2Run,settings.preprocessing); 
            
            if strcmp(dirfiles,'')
                warning('No files selected!');
                return;
            else
                nfiles     = length(dirfiles);
                fprintf('We got %s files! \n', num2str(nfiles));
            end
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = cellstr(dirfiles); %JB For every run scans are added to the matlabbatch
            
            clear dirfiles
            
            %......Include realignment parameters for each run
            if fla.realignmentParameters_flag == 1 % Start realignment
                
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
                
                % Load alignment parameters
                AlignmentFile    = spm_select('List',fullfile(realigned_data_dir,runs(r).name),'^rp_.*.txt'); % get realignment parameters
                raParamValues    = load(string(fullfile(realigned_data_dir,runs(r).name,AlignmentFile))); % load values
                
                % Add realignment parameters as regressors of no interest.
                for P = 1:length (raParamNames)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(P).name = raParamNames(P,:);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(P).val  = raParamValues(:,P);
                end
            end % End realignment
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % % % Define parameters such as conditions, button presses... % % % % % %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if do.loadlog == 1 % get the mat files
            log_matfiles        = spm_select('FPList', psyphysic_data_dir, 'log.mat'); % SELECT MAT FILES
        end
        
        for r = 1:nruns % For each run
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi     = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf       = 128;        % High-pass filter
            
            if do.loadlog == 1
                
                load(strrep(log_matfiles(r,:),' ',''));
                
                VideoNames          = log.data.Condition;
                fla.conditionNames  = convertStringsToChars(VideoNames);
                fla.nconditions     = length(fla.conditionNames);
                fla.trial_onset     = {}; 
                fla.trial_duration  = {};
                for reg = 1:length (VideoNames)
                    % get the video name so you can extract the Timepoint
                    % of magic moment from the 'do' struct
                    % The videoname may contain an '_F' for the flip
                    % condition. Therefore we need to get the normal video
                    % name
                    currentVideo = VideoNames(reg);
                    fla.trial_onset{end+1}      = log.data.VideoStart(contains(log.data.Condition,currentVideo))-log.data.VideoStart(1)+dummiesTime;
                    
                    % maybe remove the duration
                    fla.trial_duration{end+1}   = log.data.Rating_stimOn(contains(log.data.Condition,currentVideo))...
                        -fla.trial_onset{end};
                end
            else % don't get the matfile log files.
                
              %PLACE HOLDER
              
            end % End load mat files
            
            % Condition names:
            for cc = 1:fla.nconditions % for CurrentCondition
                matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).name     = fla.conditionNames{cc};
                matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).onset    = fla.trial_onset{cc};
                matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).duration = fla.trial_duration{cc};
                
                % Equal for all conditions (Except trialDuration for button presses)
                matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).pmod        = struct('name', {}, 'param', {}, 'poly', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).tmod        = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(cc).orth        = 1;
                                                
            end
        end
      
    %% Specify design in SPM.mat
    jobs = matlabbatch;
    spm('defaults', 'FMRI');
    spm_jobman('run', jobs);
    clear jobs; clear matlabbatch;         
    end
          
    %% Estimate design:
    % PRG
    if do.estimate == 1
        matlabbatch{1}.spm.stats.fmri_est.spmmat           = {fullfile(beta_loc, 'SPM.mat')};
        matlabbatch{1}.spm.stats.fmri_est.write_residuals  = 0;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        
        jobs = matlabbatch;
        spm('defaults', 'FMRI');
        spm_jobman('run', jobs);
        clear jobs; clear matlabbatch;
    end
    toc;
end 
