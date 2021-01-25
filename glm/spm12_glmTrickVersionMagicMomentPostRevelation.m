%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for the beta estimates of the "Magic" Experiment.
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
% each magic video post revelation. In the magic videos the special moment
% was the moment the magical effect happened, Those moments were selected 
% by Vincent Plikat.
% The resulting betas will be used for the decoding analysis. 
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

%% The selected frames (meaned) Pablo and I selected to be the frames where the magic effect happens -- taken from my python script
% read in the .mat file, that contains the information about the magic
% moment

fprintf('Please select the .mat file, that contains the information about the special moments.\n\n')
matfile_dir = uigetfile(pwd, 'Select the .mat file');

load(matfile_dir);

fps                         = 25;   % framerate of our Video
frame_time                  = 1/fps;% Time every frame is presented
preTime                     = 1;    % Time in seconds prior to magic effect we include in our glm 
postTime                    = 4;    % Time in seconds after the magic effect we include in our glm 

%% Define what to do
do.SpecifyDesign      = 1;
do.loadlog            = 1; % load LOG files!
do.estimate           = 1;

%% Settings
settings.matprefix = input (['Please specify the prefix of your participant data.\n' ...
    '(like p for participant or s for subject. It has to be unique so that only subject folders are selected):\n'],'s');
settings.preprocessing = ['^auf' '.*nii']; % ['^s2au' '.*nii']: realigned, slice time corrected, 6mm-smoothed data
settings.num_presentations_per_video = 2;

data_dir = fullfile(derived_dir,'slice_time_corrected');
folders = dir(fullfile(data_dir,[settings.matprefix, '*']));
subNames = {folders(:).name}; 
% create a "first_level_analysis" folder
spm_mkdir(fullfile(derived_dir), 'spm12flaModelTrickVersionMagicMomentPostRevelation', subNames);

fla.realignmentParameters_flag  = 1;

%% Get experimental design parameters
fla.conditionNames  = {
    'Appear1_Magic'     ; 'Appear2_Magic'; ...
    'Change1_Magic'     ; 'Change2_Magic'; ...
    'Vanish1_Magic'     ; 'Vanish2_Magic';
    };
fla.nconditions     = length(fla.conditionNames);
fla.realignmentParameters_flag = 1;


for s = 1:length(subNames)
    %% Define where to store and the results and where to look for functional and anatomical data
    beta_loc            = fullfile(derived_dir,'spm12flaModelTrickVersionMagicMomentPostRevelation',subNames{s});
    realigned_data_dir  = fullfile(derived_dir,'realigned',subNames{s},'func'); % TODO: make more elegant
    psyphysic_data_dir  = fullfile(derived_dir,'PsychoPhysic',subNames{s});
    runs                = dir(fullfile(smoothed_data_dir,'run*'));
    nruns               = length(runs); % Number of Runs
    
    %% Model especification of FLA
    if do.SpecifyDesign
        
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
        counter                                             = 0;                    % A counter that is used when we are in a run, we really want to use
        
        for r = 1:nruns
            Path2Run                = fullfile(realignment_dir,runs(r).name); 
            Path2RealignmentFiles   = fullfile(realignment_dir, runs(r).name);
            dirfiles                = spm_select('FPList',Path2Run,settings.preprocessing); 
            
            % We are only interested in the last two runs of every block.
            % Four runs per block 
            if mod (r,4) == 3 || mod (r,4) == 0
                counter     = counter +1;       % needed to fill the matlab batch
                dirfiles    = spm_select('FPList',Path2Run,settings.preprocessing); 

                if strcmp(dirfiles,'')
                    warning('No files selected!');
                    return;
                else
                    nfiles     = length(dirfiles);
                    fprintf('We got %s files! \n', num2str(nfiles));
                end

                matlabbatch{1}.spm.stats.fmri_spec.sess(counter).scans = cellstr(dirfiles); 

                clear dirfiles

                %......Include realignment parameters for each run
                if fla.realignmentParameters_flag == 1

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
                    AlignmentFile    = spm_select('List',Path2RealignmentFiles,'^rp_.*.txt'); % get realignment parameters
                    raParamValues    = load(string(fullfile(Path2RealignmentFiles,AlignmentFile))); % load values

                    % Add realignment parameters as regressors of no interest.
                    for P = 1:length (raParamNames)
                        matlabbatch{1}.spm.stats.fmri_spec.sess(counter).regress(P).name = raParamNames(P,:);
                        matlabbatch{1}.spm.stats.fmri_spec.sess(counter).regress(P).val  = raParamValues(:,P);
                    end
                end

            end

        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % % % Define parameters such as conditions, button presses... % % % % % %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if do.loadlog == 1 % get the mat files
                log_matfiles        = spm_select('FPList', psyphysic_data_dir, 'log.mat'); % SELECT MAT FILES
            end
        

            matlabbatch{1}.spm.stats.fmri_spec.sess(counter).multi     = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(counter).multi_reg = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(counter).hpf       = 128;        % High-pass filter
            
            fla.trial_duration={}; fla.trial_onset={};
            % ask if we are in the right run (pre revelation -> 1 & 2
            if do.loadlog == 1

                load(strrep(log_matfiles(r,:),' ',''));

                trial_onset={}; trial_duration={};
                for reg = 1:length (fla.conditionNames)
                    % get the video name so you can extract the Timepoint
                    % of magic moment from the 'do' struct
                    VideoNames                  = log.data.Condition(contains(log.data.Condition,fla.conditionNames(reg)));
                    % The videoname may contain an '_F' for the flip
                    % condition. Therefore we need to get the normal video
                    % name
                    if contains(VideoNames{1},'_F')
                        VideoName=VideoNames{1};
                        VideoName=VideoName(1:end-2);
                    else
                        VideoName=VideoNames{1};
                    end
                    SpecialMoment               = do.all_frames_of_effect(contains(do.ListOfVideos,VideoName));
                    SpecialMomentOnset          = SpecialMoment{1}*frame_time;
                    trial_onset{end+1}      = log.data.VideoStart(contains(log.data.Condition,fla.conditionNames(reg)))-log.data.VideoStart(1)+ ...
                        nDummies*matlabbatch{1}.spm.stats.fmri_spec.timing.RT + SpecialMomentOnset - preTime;
                    trial_duration{end+1}   = repmat(preTime+postTime, length(trial_onset{1}),1); 
                end

                fla.trial_onset     = num2cell(vertcat(trial_onset{:})); 
                fla.trial_duration  = num2cell(vertcat(trial_duration{:}));

            else % don't get the matfile log files.

              %PLACE HOLDER

            end


            % Condition names:
            for cc = 1:fla.nconditions % for CurrentCondition
                % since for every condition name we have two
                % presentations, we iterate twice for every condition
                for p = 1:settings.num_presentations_per_video
                    matlabbatch{1}.spm.stats.fmri_spec.sess(counter).cond((cc-1)*settings.num_presentations_per_video+p).name     = fla.conditionNames{cc};
                    matlabbatch{1}.spm.stats.fmri_spec.sess(counter).cond((cc-1)*settings.num_presentations_per_video+p).onset    = fla.trial_onset{(cc-1)*settings.num_presentations_per_video+p};
                    matlabbatch{1}.spm.stats.fmri_spec.sess(counter).cond((cc-1)*settings.num_presentations_per_video+p).duration = fla.trial_duration{(cc-1)*settings.num_presentations_per_video+p};

                    % Equal for all conditions (Except trialDuration for button presses)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(counter).cond((cc-1)*settings.num_presentations_per_video+p).pmod        = struct('name', {}, 'param', {}, 'poly', {});
                    matlabbatch{1}.spm.stats.fmri_spec.sess(counter).cond((cc-1)*settings.num_presentations_per_video+p).tmod        = 0;
                    matlabbatch{1}.spm.stats.fmri_spec.sess(counter).cond((cc-1)*settings.num_presentations_per_video+p).orth        = 1;
                end
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
    
end