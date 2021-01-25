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

clear matlabbatch; % first clear unfinished batches.
tic; % start script.

%% Define important details of your file structure and location
% Set root directory
fprintf(['Please select your project folder.'...
    '(ideally it should contain a folder named "raw_data")\n\n'])
root_dir    = uigetdir(homedir, 'Select Project Folder');
if root_dir == 0
    error('No folder was selected --> I terminate the script')
end

% Set source_data directory. This is only needed for slice time correction
source_dir  = fullfile(root_dir,'source_data');
if ~isfolder(source_dir)
    fprintf(['It appears you do not have a "source_data" folder.\n'...
        'Please select the folder that contains your DICOMS.'])
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
do.DefContrasts       = 1;

%% Settings
settings.matprefix = input (['Please specify the prefix of your participant data.\n' ...
    '(like p for participant or s for subject. It has to be unique so that only subject folders are selected):\n'],'s');
settings.preprocessing = ['^s6wauf' '.*nii']; % realigned, slice time corrected, normalized, 6mm-smoothed data

data_dir = fullfile(derived_dir,'6smoothed');
folders = dir(fullfile(data_dir,[settings.matprefix, '*']));
subNames = {folders(:).name}; 
% create a "first_level_analysis" folder
spm_mkdir(fullfile(derived_dir), 'spm12flaWholeVideo', subNames);


%% Get experimental design parameters
fla.conditionNames  = {
    'Appear1_Magic'     ; 'Appear2_Magic'; ...
    'Vanish1_Magic'     ; 'Vanish2_Magic';...
    'Change1_Magic'     ; 'Change2_Magic';...
    'Appear1_Control'   ; 'Appear2_Control'; ...
    'Vanish1_Control'   ; 'Vanish2_Control';...
    'Change1_Control'   ; 'Change2_Control';...
    'Surprise'};
%     'Surprise1'         ; 'Surprise2';...
%     'Surprise3'
%     };
fla.nconditions                 = length(fla.conditionNames);
numMag = 6;
numCon = 6;
numSur = 1;
fla.realignmentParameters_flag  = 1;

for s = 1:length(subNames)
    %% Define where to store and the results and where to look for functional and anatomical data
    beta_loc            = fullfile(derived_dir,'spm12flaWholeVideo',subNames{s});
    smoothed_data_dir   = fullfile(derived_dir,'6smoothed',subNames{s},'func'); 
    realigned_data_dir  = fullfile(derived_dir,'realigned',subNames{s},'func'); 
    psyphysic_data_dir  = fullfile(derived_dir,'PsychoPhysic',subNames{s});
    runs                = dir(fullfile(smoothed_data_dir,'run*'));
    nruns               = length(runs); % Number of Runs
    
    %% Model specification of FLA
    if do.SpecifyDesign == 1 % Model specification
        
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
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = cellstr(dirfiles); %For every run scans are added to the matlabbatch
            
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
                numRaPara = length(raParamNames);
                
                % Load alignment parameters
                AlignmentFile    = spm_select('List',fullfile(realigned_data_dir,runs(r).name),'^rp_.*.txt'); % get realignment parameters
                raParamValues    = load(string(fullfile(realigned_data_dir,runs(r).name,AlignmentFile))); % load values
                
                % Add realignment parameters as regressors of no interest.
                for P = 1:length (raParamNames)
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(P).name = raParamNames(P,:);
                    matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress(P).val  = raParamValues(:,P);
                end
            end
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
                    fla.trial_onset{end+1}      = log.data.VideoStart(contains(log.data.Condition,fla.conditionNames(reg))) - log.data.VideoStart(1) + ...
                        nDummies * matlabbatch{1}.spm.stats.fmri_spec.timing.RT;
                    fla.trial_duration{end+1}   = log.data.Rating_stimOn(contains(log.data.Condition,fla.conditionNames(reg)))...
                        -log.data.VideoStart(contains(log.data.Condition,fla.conditionNames(reg)));
                end
            else % don't get the matfile log files.
                
              %PLACE HOLDER
              
            end
            
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
    end % model specification
          
    %% Estimate design:
    if do.estimate == 1
        matlabbatch{1}.spm.stats.fmri_est.spmmat           = {fullfile(beta_loc, 'SPM.mat')}; 
        matlabbatch{1}.spm.stats.fmri_est.write_residuals  = 0;
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        
        jobs = matlabbatch;
        spm('defaults', 'FMRI');
        spm_jobman('run', jobs);
        clear jobs; clear matlabbatch;
    end % end estimate
    
    
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
            'Magic After > Magic Before';...%4
            'Magic > Surprise Before';...  %5
            'Surpise > Magic Before';...  %6
            'Surprise > NoMagic'; ... % 7
            'NoMagic > Surprise'; ... %8
            'Magic > NoMagic After'; ... %9
            'NoMagic > Magic After'; ... %10 
            'Magic > Surprise After'; ... %11
            'Surpise > Magic After'; ...  %12
            'MagPre-ConPre vs MagPost-ConPost';... %13
            'MagPost-ConPost vs MagPre-ConPre';... %14
            'Appear Before > Appear After';...%15
            'Vanish Before > Vanish After';...%16
            'Change Before > Change After';...%17
            'Appear After > Appear Before';...%18
            'Vanish After > Vanish Before';...%19
            'Change After > Change Before';...%20
            'Appear > Control Before'; ... %21
            'Control > Appear Before'; ... %22
            'Vanish > Control Before'; ... %23
            'Control > Vanish Before'; ... %24
            'Change > Control Before'; ... %25
            'Control > Change Before'; ... %26
            'Appear > Surprise Before'; ... %27
            'Vanish > Surprise Before'; ... %28
            'Change > Surprise Before'; ... %29
            'Appear > Control After'; ... %30
            'Vanish > Control After'; ... %31
            'Change > Control After'; ... %32
            'Appear > Surprise After'; ... %33
            'Vanish > Surprise After'; ... %34
            'Change > Surprise After';... %35
            % Interaction effects
            'AppPre-ConPre vs AppPost-ConPost';... %36
            'VanPre-ConPre vs Vanpost-ConPost';... %37
            'ChaPre-ConPre vs ChaPost-ConPost';... %38
            'AppPost-ConPost vs AppPre-ConPre';... %39
            'Vanpost-ConPost vs Vanpost-ConPost';... %40
            'ChaPost-ConPost vs ChaPre-ConPre';... %41
            % Contrats to outrule the timeconfound by comparing run 1vs2
            % and run 2vs3 - the same time difference, but the first is pre
            % vs pre and the other is pre vs post
            'Magic PreVsPre (run 1vs2)';... %42
            'Magic PreVsPost (run 2vs3)'... %43
            };
%       PreRevelation       Magic videos            NoMagic             Surprise           Realignment      PostRevelation  Magic videos         NoMagic             Surprise           Realignment
        C1 = repmat([repmat([ones(1,numMag)         ones(1,numCon)*(-1) zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([zeros(1,numMag)     zeros(1,numCon)     zeros(1,numSur)     zeros(1,numRaPara)],1,2)],1,3);
        C2 = repmat([repmat([ones(1,numMag)*(-1)    ones(1,numCon)      zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([zeros(1,numMag)     zeros(1,numCon)     zeros(1,numSur)     zeros(1,numRaPara)],1,2)],1,3);
        C3 = repmat([repmat([ones(1,numMag)         zeros(1,numCon)     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([ones(1,numMag)*(-1) zeros(1,numCon)     zeros(1,numSur)     zeros(1,numRaPara)],1,2)],1,3);
        C4 = repmat([repmat([ones(1,numMag)*(-1)    zeros(1,numCon)     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([ones(1,numMag)      zeros(1,numCon)     zeros(1,numSur)     zeros(1,numRaPara)],1,2)],1,3);
        C5 = repmat([repmat([ones(1,numMag)         zeros(1,numCon)     ones(1,numSur)*(-6) zeros(1,numRaPara)],1,2) repmat([zeros(1,numMag)     zeros(1,numCon)     zeros(1,numSur)     zeros(1,numRaPara)],1,2)],1,3);
        C6 = repmat([repmat([ones(1,numMag)*(-1)    zeros(1,numCon)     ones(1,numSur)*6    zeros(1,numRaPara)],1,2) repmat([zeros(1,numMag)     zeros(1,numCon)     zeros(1,numSur)     zeros(1,numRaPara)],1,2)],1,3);
        C7 = repmat([repmat([zeros(1,numMag)        ones(1,numCon)*(-1) ones(1,numSur)*6    zeros(1,numRaPara)],1,2) repmat([zeros(1,numMag)     ones(1,numCon)*(-1) ones(1,numSur)*6    zeros(1,numRaPara)],1,2)],1,3);
        C8 = repmat([repmat([zeros(1,numMag)        ones(1,numCon)      ones(1,numSur)*(-6) zeros(1,numRaPara)],1,2) repmat([zeros(1,numMag)     ones(1,numCon)      ones(1,numSur)*(-6) zeros(1,numRaPara)],1,2)],1,3);
        C9 = repmat([repmat([zeros(1,numMag)        zeros(1,numCon)     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([ones(1,numMag)      ones(1,numCon)*(-1) zeros(1,numSur)     zeros(1,numRaPara)],1,2)],1,3);
        C10= repmat([repmat([zeros(1,numMag)        zeros(1,numCon)     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([ones(1,numMag)*(-1) ones(1,numCon)      zeros(1,numSur)     zeros(1,numRaPara)],1,2)],1,3);
        C11= repmat([repmat([zeros(1,numMag)        zeros(1,numCon)     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([ones(1,numMag)      zeros(1,numCon)     ones(1,numSur)*(-6) zeros(1,numRaPara)],1,2)],1,3);
        C12= repmat([repmat([zeros(1,numMag)        zeros(1,numCon)     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([ones(1,numMag)*(-1) zeros(1,numCon)     ones(1,numSur)*6    zeros(1,numRaPara)],1,2)],1,3);
%       Interaction Effects
        C13= repmat([repmat([ones(1,numMag)         ones(1,numCon)*(-1) zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([ones(1,numMag)*(-1) ones(1,numCon)      zeros(1,numSur)     zeros(1,numRaPara)],1,2)],1,3);
        C14= repmat([repmat([ones(1,numMag)*(-1)    ones(1,numCon)      zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([ones(1,numMag)      ones(1,numCon)*(-1) zeros(1,numSur)     zeros(1,numRaPara)],1,2)],1,3);
        %                   AppearM VanishM ChangeM AppearC VanishC ChangeC Surprise            Realignment     PostRevelation  AppearM VanishM ChangeM AppearC VanishC ChangeC Surprise             Realignment
        C15= repmat([repmat([1 1    0 0     0 0     0 0     0 0     0 0     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([-1 -1  0 0     0 0     0 0     0 0     0 0     zeros(1,numSur)      zeros(1,numRaPara)],1,2)],1,3);
        C16= repmat([repmat([0 0    1 1     0 0     0 0     0 0     0 0     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([0 0    -1 -1   0 0     0 0     0 0     0 0     zeros(1,numSur)      zeros(1,numRaPara)],1,2)],1,3);
        C17= repmat([repmat([0 0    0 0     1 1     0 0     0 0     0 0     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([0 0    0 0     -1 -1   0 0     0 0     0 0     zeros(1,numSur)      zeros(1,numRaPara)],1,2)],1,3);
        C18= repmat([repmat([-1 -1  0 0     0 0     0 0     0 0     0 0     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([1 1    0 0     0 0     0 0     0 0     0 0     zeros(1,numSur)      zeros(1,numRaPara)],1,2)],1,3);
        C19= repmat([repmat([0 0    -1 -1   0 0     0 0     0 0     0 0     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([0 0    1 1     0 0     0 0     0 0     0 0     zeros(1,numSur)      zeros(1,numRaPara)],1,2)],1,3);
        C20= repmat([repmat([0 0    0 0     -1 -1   0 0     0 0     0 0     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([0 0    0 0     1 1     0 0     0 0     0 0     zeros(1,numSur)      zeros(1,numRaPara)],1,2)],1,3);
        C21= repmat([repmat([1 1    0 0     0 0     -1 -1   0 0     0 0     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([0 0    0 0     0 0     0 0     0 0     0 0     zeros(1,numSur)      zeros(1,numRaPara)],1,2)],1,3);
        C22= repmat([repmat([-1 -1  0 0     0 0     1 1     0 0     0 0     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([0 0    0 0     0 0     0 0     0 0     0 0     zeros(1,numSur)      zeros(1,numRaPara)],1,2)],1,3);
        C23= repmat([repmat([0 0    1 1     0 0     0 0     -1 -1   0 0     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([0 0    0 0     0 0     0 0     0 0     0 0     zeros(1,numSur)      zeros(1,numRaPara)],1,2)],1,3);
        C24= repmat([repmat([0 0    -1 -1   0 0     0 0     1 1     0 0     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([0 0    0 0     0 0     0 0     0 0     0 0     zeros(1,numSur)      zeros(1,numRaPara)],1,2)],1,3);
        C25= repmat([repmat([0 0    0 0     1 1     0 0     0 0     -1 -1   zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([0 0    0 0     0 0     0 0     0 0     0 0     zeros(1,numSur)      zeros(1,numRaPara)],1,2)],1,3);
        C26= repmat([repmat([0 0    0 0     -1 -1   0 0     0 0     1 1     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([0 0    0 0     0 0     0 0     0 0     0 0     zeros(1,numSur)      zeros(1,numRaPara)],1,2)],1,3);
        C27= repmat([repmat([1 1    0 0     0 0     0 0     0 0     0 0     ones(1,numSur)*(-2) zeros(1,numRaPara)],1,2) repmat([0 0    0 0     0 0     0 0     0 0     0 0     zeros(1,numSur)      zeros(1,numRaPara)],1,2)],1,3);
        C28= repmat([repmat([0 0    1 1     0 0     0 0     0 0     0 0     ones(1,numSur)*(-2) zeros(1,numRaPara)],1,2) repmat([0 0    0 0     0 0     0 0     0 0     0 0     zeros(1,numSur)      zeros(1,numRaPara)],1,2)],1,3);
        C29= repmat([repmat([0 0    0 0     1 1     0 0     0 0     0 0     ones(1,numSur)*(-2) zeros(1,numRaPara)],1,2) repmat([0 0    0 0     0 0     0 0     0 0     0 0     zeros(1,numSur)      zeros(1,numRaPara)],1,2)],1,3);
        C30= repmat([repmat([0 0    0 0     0 0     0 0     0 0     0 0     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([1 1    0 0     0 0     -1 -1   0 0     0 0     zeros(1,numSur)      zeros(1,numRaPara)],1,2)],1,3);
        C31= repmat([repmat([0 0    0 0     0 0     0 0     0 0     0 0     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([0 0    1 1     0 0     0 0     -1 -1   0 0     zeros(1,numSur)      zeros(1,numRaPara)],1,2)],1,3);
        C32= repmat([repmat([0 0    0 0     0 0     0 0     0 0     0 0     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([0 0    0 0     1 1     0 0     0 0     -1 -1   zeros(1,numSur)      zeros(1,numRaPara)],1,2)],1,3);
        C33= repmat([repmat([0 0    0 0     0 0     0 0     0 0     0 0     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([1 1    0 0     0 0     0 0     0 0     0 0     ones(1,numSur)*(-2)  zeros(1,numRaPara)],1,2)],1,3);
        C34= repmat([repmat([0 0    0 0     0 0     0 0     0 0     0 0     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([0 0    1 1     0 0     0 0     0 0     0 0     ones(1,numSur)*(-2)  zeros(1,numRaPara)],1,2)],1,3);
        C35= repmat([repmat([0 0    0 0     0 0     0 0     0 0     0 0     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([0 0    0 0     1 1     0 0     0 0     0 0     ones(1,numSur)*(-2)  zeros(1,numRaPara)],1,2)],1,3);
%       Interaction Effects
        C36= repmat([repmat([1 1    0 0     0 0     -1 -1   0 0     0 0     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([-1 -1  0 0     0 0     1 1     0 0     0 0     zeros(1,numSur) zeros(1,numRaPara)],1,2)],1,3);
        C37= repmat([repmat([0 0    1 1     0 0     0 0     -1 -1   0 0     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([0 0    -1 -1   0 0     0 0     1 1     0 0     zeros(1,numSur) zeros(1,numRaPara)],1,2)],1,3);
        C38= repmat([repmat([0 0    0 0     1 1     0 0     0 0     -1 -1   zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([0 0    0 0     -1 -1   0 0     0 0     1 1     zeros(1,numSur) zeros(1,numRaPara)],1,2)],1,3);
        C39= repmat([repmat([-1 -1  0 0     0 0     1 1     0 0     0 0     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([1 1    0 0     0 0     -1 -1   0 0     0 0     zeros(1,numSur) zeros(1,numRaPara)],1,2)],1,3);
        C40= repmat([repmat([0 0    -1 -1   0 0     0 0     1 1     0 0     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([0 0    1 1     0 0     0 0     -1 -1   0 0     zeros(1,numSur) zeros(1,numRaPara)],1,2)],1,3);
        C41= repmat([repmat([0 0    0 0     -1 -1   0 0     0 0     1 1     zeros(1,numSur)     zeros(1,numRaPara)],1,2) repmat([0 0    0 0     1 1     0 0     0 0     -1 -1   zeros(1,numSur) zeros(1,numRaPara)],1,2)],1,3);
        %   FirstRun Magic           NoMagic         Surprise        Realignment SecondRun   Magic                NoMagic         Surprise        Realignment Postrevelation  Magic           NoMagic         Surprise        Realignment
        C42= repmat([ones(1,numMag)  zeros(1,numCon) zeros(1,numSur) zeros(1,numRaPara)      ones(1,numMag)*(-1)  zeros(1,numCon) zeros(1,numSur) zeros(1,numRaPara)  repmat([zeros(1,numMag) zeros(1,numCon) zeros(1,numSur) zeros(1,numRaPara)],1,2)],1,3);
        %   FirstRun Magic           NoMagic         Surprise        Realignment SecondRun   Magic            NoMagic         Surprise        Realignment   ThirdRun  Magic                 NoMagic         Surprise        Realignment FourthRun   Magic           NoMagic         Surprise        Realignment
        C43= repmat([zeros(1,numMag) zeros(1,numCon) zeros(1,numSur) zeros(1,numRaPara)      ones(1,numMag)   zeros(1,numCon) zeros(1,numSur) zeros(1,numRaPara)      ones(1,numMag)*(-1)   zeros(1,numCon) zeros(1,numSur) zeros(1,numRaPara)      zeros(1,numMag) zeros(1,numCon) zeros(1,numSur) zeros(1,numRaPara)],1,3);
        
        Contrasts = [C1; C2; C3; C5; C6; C7; C8; C9; C10; C11; C12; C13; C14; C15; C16; C17; C18; C19; C20;...
            C21; C22; C23; C24; C25; C26; C27; C28; C29; C30; C31; C32; C33; C34; C35; C36; C37; C38; C39; C40; C41; C42; C43];
        
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

