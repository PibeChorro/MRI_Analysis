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

% read in the .mat file, that contains the information about the magic
% moment

fprintf('Please select the .mat file, that contains the information about the special moments.\n\n')
matfile_dir = uigetfile(pwd, 'Select the .mat file');

load(matfile_dir);

%% Define what to do
do.SpecifyDesign      = 1;
do.loadlog            = 1; % load LOG files!
do.estimate           = 1;
do.DefContrasts       = 1;

%% Settings
settings.matprefix = input (['Please specify the prefix of your participant data.\n' ...
    '(like p for participant or s for subject. It has to be unique so that only subject folders are selected):\n'],'s');
settings.preprocessing = ['^s6wauf' '.*nii']; % ['^s2au' '.*nii']: realigned, slice time corrected, 6mm-smoothed data

data_dir = fullfile(derived_dir,'6smoothed');
folders = dir(fullfile(data_dir,[settings.matprefix, '*']));

subNames = {folders(:).name}; 
% create a "first_level_analysis" folder
spm_mkdir(fullfile(derived_dir), 'spm12flaSpecialMoment', subNames);


%% Get experimental design parameters
fla.conditionNames  = {
    'Appear1_Magic'     ; 'Appear2_Magic'; ...
    'Vanish1_Magic'     ; 'Vanish2_Magic';...
    'Change1_Magic'     ; 'Change2_Magic';...
    'Appear1_Control'   ; 'Appear2_Control'; ...
    'Vanish1_Control'   ; 'Vanish2_Control';...
    'Change1_Control'   ; 'Change2_Control';...
    'Surprise1'         ; 'Surprise2';...
    'Surprise3'
    };
fla.nconditions                 = length(fla.conditionNames);
fla.realignmentParameters_flag  = 1;

%% Information about the videos
fps                         = 25;   % framerate of our Video
frame_time                  = 1/fps;% Time every frame is presented

for s = 1:length(subNames)
    %% Define where to store and the results and where to look for functional and anatomical data
    beta_loc            = fullfile(derived_dir,'spm12flaSpecialMoment',subNames{s});
    smoothed_data_dir   = fullfile(derived_dir,'6smoothed',subNames{s},'func'); 
    realigned_data_dir  = fullfile(derived_dir,'realigned',subNames{s},'func'); 
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
                
                % consider extracting dummy values
                %tmp = tmp(n_dummies+1:n_images_per_run(ses)+n_dummies, :); % Remove trailing entries from realignment txt
                
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
                
                fla.trial_onset={}; fla.trial_duration={};
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
                    fla.trial_onset{end+1}      = log.data.VideoStart(contains(log.data.Condition,fla.conditionNames(reg)))-log.data.VideoStart(1)+...
                        SpecialMomentOnset+dummiesTime;
                    
                    % maybe remove the duration
                    fla.trial_duration{end+1}   = log.data.Rating_stimOn(contains(log.data.Condition,fla.conditionNames(reg)))...
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
    
    
    %% Define contrast
    % PRG
    if do.DefContrasts == 1
        
        %% Define contrasts:
        matlabbatch{1}.spm.stats.con.spmmat = {fullfile(beta_loc, 'SPM.mat')};
        matlabbatch{1}.spm.stats.con.delete = 1;
        
        % Names:
        ContrastNames = ...
            {'Magic > NoMagic Before'; ... %1
            'NoMagic > Magic Before'; ...  %2 
            'Magic Before > Magic After';...%3
            'Appear Before > Appear After';...%3a
            'Vanish Before > Vanish After';...%3b
            'Change Before > Change After';...%3c
            'Magic After > Magic Before';...%4
            'Appear After > Appear Before';...%4a
            'Vanish After > Vanish Before';...%4b
            'Change After > Change Before';...%4c
            'Magic > Surprise Before';...  %5
            'Surpise > Magic Before';...  %6
            'Surprise > NoMagic'; ... % 7
            'NoMagic > Surprise'; ... %8
            'Appear > Control Before'; ... 9
            'Control > Appear Before'; ... 9a
            'Vanish > Control Before'; ... %10
            'Control > Vanish Before'; ... %10a
            'Change > Control Before'; ... %11
            'Control > Change Before'; ... %11a
            'Appear > Surprise Before'; ... %12
            'Vanish > Surprise Before'; ... %13
            'Change > Surprise Before'; ... %14
            'Magic > NoMagic After'; ... %15
            'NoMagic > Magic After'; ... %16 
            'Magic > Surprise After'; ... %17
            'Surpise > Magic After'; ...  %18
            'Appear > Control After'; ... %19
            'Vanish > Control After'; ... %20
            'Change > Control After'; ... %21
            'Appear > Surprise After'; ... %22
            'Vanish > Surprise After'; ... %23
            'Change > Surprise After';... %24
            % Interaction effects
            'MagPre-ConPre vs MagPost-ConPost';... %25
            'AppPre-ConPre vs AppPost-ConPost';... %26
            'VanPre-ConPre vs Vanpost-ConPost';... % 27
            'ChaPre-ConPre vs ChaPost-ConPost';... %28
            'MagPost-ConPost vs MagPre-ConPre';... %29
            'AppPost-ConPost vs AppPre-ConPre';... %30
            'Vanpost-ConPost vs Vanpost-ConPost';... %31
            'ChaPost-ConPost vs ChaPre-ConPre';... %32
            % Contrats to outrule the timeconfound by comparing run 1vs2
            % and run 2vs3 - the same time difference, but the first is pre
            % vs pre and the other is pre vs post
            'Magic PreVsPre (run 1vs2)';... %33
            'Magic PreVsPost (run 2vs3)'... %34
            };
            
        
        %       PreRevelation   Magic videos    NoMagic         Surprise        Realignment                         PostRevelation  Magic videos    NoMagic             Surprise       Realignment
        C1 = repmat([repmat(    [ones(1,6)      ones(1,6)*(-1)  zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [zeros(1,6)     zeros(1,6)      zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        C2 = repmat([repmat(    [ones(1,6)*(-1) ones(1,6)       zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [zeros(1,6)     zeros(1,6)      zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        C3 = repmat([repmat(    [ones(1,6)      zeros(1,6)      zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [ones(1,6)*(-1) zeros(1,6)      zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        %                       AppearM VanishM ChangeM AppearC VanishC ChangeC Surprise        Realignment                         PostRevelation  AppearM VanishM ChangeM AppearC     VanishC ChangeC     Surprise        Realignment
        C3a= repmat([repmat(    [1 1    0 0     0 0     0 0     0 0     0 0     zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [-1 -1  0 0     0 0     0 0         0 0     0 0         zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        C3b= repmat([repmat(    [0 0    1 1     0 0     0 0     0 0     0 0     zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [0 0    -1 -1   0 0     0 0         0 0     0 0         zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        C3c= repmat([repmat(    [0 0    0 0     1 1     0 0     0 0     0 0     zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [0 0    0 0     -1 -1   0 0         0 0     0 0         zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        %       PreRevelation   Magic videos    NoMagic         Surprise        Realignment                         PostRevelation  Magic videos    NoMagic             Surprise       Realignment
        C4 = repmat([repmat(    [ones(1,6)*(-1) zeros(1,6)      zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [ones(1,6)      zeros(1,6)      zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        %                       AppearM VanishM ChangeM AppearC VanishC ChangeC Surprise        Realignment                         PostRevelation  AppearM VanishM ChangeM AppearC     VanishC ChangeC     Surprise        Realignment
        C4a= repmat([repmat(    [-1 -1  0 0     0 0     0 0     0 0     0 0     zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [1 1    0 0     0 0     0 0         0 0     0 0         zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        C4b= repmat([repmat(    [0 0    -1 -1   0 0     0 0     0 0     0 0     zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [0 0    1 1     0 0     0 0         0 0     0 0         zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        C4c= repmat([repmat(    [0 0    0 0     -1 -1   0 0     0 0     0 0     zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [0 0    0 0     1 1     0 0         0 0     0 0         zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        %       PreRevelation   Magic videos    NoMagic         Surprise        Realignment                         PostRevelation  Magic videos    NoMagic             Surprise       Realignment        
        C5 = repmat([repmat(    [ones(1,6)      zeros(1,6)      ones(1,3)*(-2)  zeros(1,length(raParamNames))],1,2)     repmat(     [zeros(1,6)     zeros(1,6)      zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        C6 = repmat([repmat(    [ones(1,6)*(-1) zeros(1,6)      ones(1,3)*2     zeros(1,length(raParamNames))],1,2)     repmat(     [zeros(1,6)     zeros(1,6)      zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        C7 = repmat([repmat(    [zeros(1,6)     ones(1,6)*(-1)  ones(1,3)*2     zeros(1,length(raParamNames))],1,2)     repmat(     [zeros(1,6)     ones(1,6)*(-1)  ones(1,3)*2     zeros(1,length(raParamNames))],1,2)],1,3);
        C8 = repmat([repmat(    [zeros(1,6)     ones(1,6)       ones(1,3)*(-2)  zeros(1,length(raParamNames))],1,2)     repmat(     [zeros(1,6)     ones(1,6)       ones(1,3)*(-2)  zeros(1,length(raParamNames))],1,2)],1,3);
        %                       AppearM VanishM ChangeM AppearC VanishC ChangeC Surprise        Realignment                         PostRevelation  AppearM VanishM ChangeM AppearC     VanishC ChangeC     Surprise        Realignment
        C9 = repmat([repmat(    [1 1    0 0     0 0     -1 -1   0 0     0 0     zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [0 0    0 0     0 0     0 0         0 0     0 0         zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        C9a= repmat([repmat(    [-1 -1  0 0     0 0     1 1     0 0     0 0     zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [0 0    0 0     0 0     0 0         0 0     0 0         zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        C10= repmat([repmat(    [0 0    1 1     0 0     0 0     -1 -1   0 0     zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [0 0    0 0     0 0     0 0         0 0     0 0         zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        C10a=repmat([repmat(    [0 0    -1 -1   0 0     0 0     1 1     0 0     zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [0 0    0 0     0 0     0 0         0 0     0 0         zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        C11= repmat([repmat(    [0 0    0 0     1 1     0 0     0 0     -1 -1   zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [0 0    0 0     0 0     0 0         0 0     0 0         zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        C11a=repmat([repmat(    [0 0    0 0     -1 -1   0 0     0 0     1 1     zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [0 0    0 0     0 0     0 0         0 0     0 0         zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        C12= repmat([repmat(    [3 3    0 0     0 0     0 0     0 0     0 0     ones(1,3)*(-1)  zeros(1,length(raParamNames))],1,2)     repmat(     [0 0    0 0     0 0     0 0         0 0     0 0         ones(1,3)*(-1)  zeros(1,length(raParamNames))],1,2)],1,3);
        C13= repmat([repmat(    [0 0    3 3     0 0     0 0     0 0     0 0     ones(1,3)*(-1)  zeros(1,length(raParamNames))],1,2)     repmat(     [0 0    0 0     0 0     0 0         0 0     0 0         ones(1,3)*(-1)  zeros(1,length(raParamNames))],1,2)],1,3);
        C14= repmat([repmat(    [0 0    0 0     3 3     0 0     0 0     0 0     ones(1,3)*(-1)  zeros(1,length(raParamNames))],1,2)     repmat(     [0 0    0 0     0 0     0 0         0 0     0 0         ones(1,3)*(-1)  zeros(1,length(raParamNames))],1,2)],1,3);
        %       PreRevelation       Magic videos    NoMagic     Surprise    Realignment                             PostRevelation  Magic videos    NoMagic          Surprise         Realignment
        C15= repmat([   repmat(     [zeros(1,6)     zeros(1,6)  zeros(1,3)  zeros(1,length(raParamNames))],1,2)     repmat(         [ones(1,6)      ones(1,6)*(-1)   zeros(1,3)       zeros(1,length(raParamNames))],1,2)],1,3);
        C16= repmat([   repmat(     [zeros(1,6)     zeros(1,6)  zeros(1,3)  zeros(1,length(raParamNames))],1,2)     repmat(         [ones(1,6)*(-1) ones(1,6)        zeros(1,3)       zeros(1,length(raParamNames))],1,2)],1,3);
        C17= repmat([   repmat(     [zeros(1,6)     zeros(1,6)  zeros(1,3)  zeros(1,length(raParamNames))],1,2)     repmat(         [ones(1,6)      zeros(1,6)       ones(1,3)*(-2)   zeros(1,length(raParamNames))],1,2)],1,3);
        C18= repmat([   repmat(     [zeros(1,6)     zeros(1,6)  zeros(1,3)  zeros(1,length(raParamNames))],1,2)     repmat(         [ones(1,6)*(-1) zeros(1,6)       ones(1,3)*2      zeros(1,length(raParamNames))],1,2)],1,3);
        %       PreRevelation   AppearM VanishM ChangeM AppearC VanishC     ChangeC     Surprise    Realignment                         PostRevelation  AppearM VanishM ChangeM AppearC     VanishC     ChangeC     Surprise        Realignment                     
        C19= repmat([   repmat( [0 0    0 0     0 0     0 0     0 0         0 0         zeros(1,3)  zeros(1,length(raParamNames))],1,2)    repmat(      [1 1    0 0     0 0     -1 -1   0 0         0 0         zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        C20= repmat([   repmat( [0 0    0 0     0 0     0 0     0 0         0 0         zeros(1,3)  zeros(1,length(raParamNames))],1,2)    repmat(      [0 0    1 1     0 0     0 0         -1 -1   0 0         zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        C21= repmat([   repmat( [0 0    0 0     0 0     0 0     0 0         0 0         zeros(1,3)  zeros(1,length(raParamNames))],1,2)    repmat(      [0 0    0 0     1 1     0 0         0 0         -1 -1   zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        C22= repmat([   repmat( [0 0    0 0     0 0     0 0     0 0         0 0         zeros(1,3)  zeros(1,length(raParamNames))],1,2)    repmat(      [3 3    0 0     0 0     0 0         0 0         0 0         ones(1,3)*(-2)  zeros(1,length(raParamNames))],1,2)],1,3);
        C23= repmat([   repmat( [0 0    0 0     0 0     0 0     0 0         0 0         zeros(1,3)  zeros(1,length(raParamNames))],1,2)    repmat(      [0 0    3 3     0 0     0 0         0 0         0 0         ones(1,3)*(-2)  zeros(1,length(raParamNames))],1,2)],1,3);
        C24= repmat([   repmat( [0 0    0 0     0 0     0 0     0 0         0 0         zeros(1,3)  zeros(1,length(raParamNames))],1,2)    repmat(      [0 0    0 0     3 3     0 0         0 0         0 0         ones(1,3)*(-2)  zeros(1,length(raParamNames))],1,2)],1,3);
        %       PreRevelation   Magic videos    NoMagic         Surprise        Realignment                         PostRevelation  Magic videos    NoMagic         Surprise        Realignment
        C25= repmat([repmat(    [ones(1,6)      ones(1,6)*(-1)  zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [ones(1,6)*(-1) ones(1,6)       zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        %       PreRevelation   AppearM VanishM ChangeM AppearC VanishC ChangeC Surprise        Realignment                         PostRevelation  AppearM VanishM ChangeM AppearC     VanishC ChangeC     Surprise        Realignment
        C26= repmat([repmat(    [1 1    0 0     0 0     -1 -1   0 0     0 0     zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [-1 -1  0 0     0 0     1 1         0 0     0 0         zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        C27= repmat([repmat(    [0 0    1 1     0 0     0 0     -1 -1   0 0     zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [0 0    -1 -1   0 0     0 0         1 1     0 0         zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        C28= repmat([repmat(    [0 0    0 0     1 1     0 0     0 0     -1 -1   zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [0 0    0 0     -1 -1   0 0         0 0     1 1         zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        %       PreRevelation   Magic videos    NoMagic         Surprise        Realignment                         PostRevelation  Magic videos    NoMagic         Surprise        Realignment
        C29= repmat([repmat(    [ones(1,6)*(-1) ones(1,6)       zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [ones(1,6)      ones(1,6)*(-1)  zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        %                       AppearM VanishM ChangeM AppearC VanishC ChangeC Surprise        Realignment                         PostRevelation  AppearM VanishM ChangeM AppearC     VanishC ChangeC     Surprise        Realignment
        C30= repmat([repmat(    [-1 -1  0 0     0 0     1 1     0 0     0 0     zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [1 1    0 0     0 0     -1 -1       0 0     0 0         zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        C31= repmat([repmat(    [0 0    -1 -1   0 0     0 0     1 1     0 0     zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [0 0    1 1     0 0     0 0         -1 -1   0 0         zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        C32= repmat([repmat(    [0 0    0 0     -1 -1   0 0     0 0     1 1     zeros(1,3)      zeros(1,length(raParamNames))],1,2)     repmat(     [0 0    0 0     1 1     0 0         0 0     -1 -1       zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        %       PreRevelation   Magic       NoMagic     Surprise        Realignment            SecondRun   Magic            NoMagic     Surprise        Realignment                      Postrevelation     Magic           NoMagic  Surprise        Realignment
        C33= repmat([           ones(1,6)   zeros(1,6)  zeros(1,3)      zeros(1,length(raParamNames))      ones(1,6)*(-1)   zeros(1,6)  zeros(1,3)      zeros(1,length(raParamNames))           repmat([    zeros(1,6)   zeros(1,6)  zeros(1,3)      zeros(1,length(raParamNames))],1,2)],1,3);
        %       FirstRun        Magic       NoMagic     Surprise        Realignment            SecondRun   Magic            NoMagic     Surprise        Realignment                      ThirdRun   Magic           NoMagic     Surprise    Realignment                     FourthRun   Magic           NoMagic  Surprise        Realignment
        C34= repmat([           zeros(1,6)  zeros(1,6)  zeros(1,3)      zeros(1,length(raParamNames))      ones(1,6)        zeros(1,6)  zeros(1,3)      zeros(1,length(raParamNames))               ones(1,6)*(-1)  zeros(1,6)  zeros(1,3)  zeros(1,length(raParamNames))               zeros(1,6)  zeros(1,6)  zeros(1,3)      zeros(1,length(raParamNames))],1,3);
        Contrasts = [C1; C2; C3; C3a; C3b; C3c; C4; C4a; C4b; C4c; C5; C6; C7; C8; C9; C9a; C10; C10a; C11; C11a; C12; C13; C14; C15; C16; C17; C18; C19; C20;...
            C21; C22; C23; C24; C25; C26; C27; C28; C29; C30; C31; C32; C33; C34];
        
                
        % safety net: check if sum of contrasts is 0
        if any(sum(Contrasts,2))
            error('One of the contrasts sum is not zero')
        end
        
        % Write in spm-struct:
        for C = 1:size(Contrasts, 1)
            matlabbatch{1}.spm.stats.con.consess{C}.tcon.name       = ContrastNames{C,:};
            matlabbatch{1}.spm.stats.con.consess{C}.tcon.sessrep    = 'none'; %'replsc' if repeat for sessions
            matlabbatch{1}.spm.stats.con.consess{C}.tcon.convec     = Contrasts(C,:); % or weights?
        end
                
        spm_jobman('run', matlabbatch);
        clear matlabbatch;
        
    end % define contrast
    
    toc;
     
end 
