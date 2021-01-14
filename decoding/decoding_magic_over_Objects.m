%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for decoding analysis of the "Magic" Experiment.
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
% In this analysis we train a decoder on 2/3 of the betas created by a glm that
% estimated those betas based on functional MRI data during the
% presentation of magic videos pre revelation and test the prediction based
% on the remaining 1/3 of the betas. In this script we use cross validation
% to improve the performance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = decoding_magic_over_Objects()

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


%% Settings
settings.matprefix = input (['Please specify the prefix of your participant data.\n' ...
    '(like p for participant or s for subject. It has to be unique so that only subject folders are selected):\n'],'s');

folders = dir(fullfile(derived_dir,[settings.matprefix, '*']));
subNames = {folders(:).name}; 
% create a "first_level_analysis" folder
spm_mkdir(derived_dir, subNames, 'MRI/analysis/decodingAnalysis/SpecialMoment');

% Set the label names to the regressor names which you want to use for
% decoding, e.g. 'button left' and 'button right'
labelnames = {
    'Appear1_Magic'     ; 'Appear2_Magic'; ...
    'Vanish1_Magic'     ; 'Vanish2_Magic';...
    'Change1_Magic'     ; 'Change2_Magic';
    };

% This means that the 2 conditions defined above are sorted depending on
% the dedicated label:
labels = [1; 1; 2; 2; 3; 3];

% Set the name for the output directory where data will be saved
output_name = 'OverObjects_MagicMomentSeachlight';

for s = 1:length(subNames)
    % Set the filepath where your SPM.mat and all related betas are
    beta_loc = fullfile(derived_dir, subNames{s}, 'MRI/analysis/glm/TrickVersionMagicMoment'); %Fill in beta folder name

    %% Settings for the TDT
    % Get TDT defaults
    cfg = decoding_defaults;
    cfg.results.overwrite = 1;
    % Set the analysis that should be performed (default is 'searchlight')
    cfg.analysis = 'searchlight';
    % cfg.testmode = 1; % use if you just want a quick test, calculates 1 SL

    % Set the output directory where data will be saved
    cfg.results.dir = fullfile(derived_dir, subNames{s}, 'MRI/analysis/decodingAnalysis/SpecialMoment', output_name);

    % Set the filename of your brain mask (or your ROI masks as cell matrix) 
    cfg.files.mask = fullfile(beta_loc, 'mask.nii');

    regressor_names = design_from_spm(beta_loc);
    cfg = decoding_describe_data(cfg,labelnames,labels,regressor_names,beta_loc);
    % Set additional parameters manually if you want (see decoding.m or
    % decoding_defaults.m). Below some example parameters that you might want 
    % to use:

    % cfg.searchlight.unit = 'mm';
    % cfg.searchlight.radius = 12; % this will yield a searchlight radius of 12mm.
    % cfg.searchlight.spherical = 1;
    % cfg.verbose = 2; % you want all information to be printed on screen
    % cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; 
    cfg.results.output = {'accuracy_minus_chance','binomial_probability'};


    % Decide whether you want to see the searchlight/ROI/... during decoding
    cfg.plot_selected_voxels = 100; % 0: no plotting, 1: every step, 2: every second step, 100: every hundredth step...

    %% Nothing needs to be changed below for a standard leave-one-run out cross
    %% validation analysis.


    %% Select the classifier
    % Standard is LibSVM:
    cfg.decoding.software = 'libsvm';
    cfg.decoding.method   = 'classification_kernel';

    % Another classifier is liblinear:
    %cfg.decoding.software = 'liblinear';
    %cfg.decoding.method   = 'classification';
    %cfg.decoding.train.classification.model_parameters = '-s 2 -c 1 -q';
    %cfg.decoding.test.classification.model_parameters  = '-q';

    %% Scale the data:
    cfg.scale.method     = 'z';  % 'z' 'min0max1', 'mean', 'none', etc.
    cfg.scale.estimation = 'all'; %'all', 'all_iter', 'across', 'separate', or 'none'.
    cfg.scale.cutoff     = [-2 2]; %before: [-inf inf];

    %% Feature transformation
    % cfg.feature_transformation.method = 'PCA';
    % cfg.feature_transformation.estimation = 'across';
    % cfg.feature_transformation.critical_value = 0.001;
    % cfg.feature_transformation.scale.method = 'none'; % mean, none, z, min0max1, etc.
    % cfg.feature_transformation.scale.cutoff = [-inf inf];

    %% Feature selection
    % cfg.feature_selection.method = 'embedded';    
    % cfg.feature_selection.embedded = 'RFE';
    % cfg.feature_selection.n_vox = [50:10:100];
    % cfg.feature_selection.nested_n_vox = [50:10:100];
    % cfg.feature_selection.direction = 'backward';

    %% Parameter selection
    %cfg.parameter_selection.method = 'grid'; % only method implemented for now in TDT  
    %cfg.parameter_selection.parameters = {'-c'};
    %cfg.parameter_selection.parameter_range = {[0.0001 0.001 0.01 0.1 1 10]};
    %cfg.parameter_selection.format.name = 'string_number';
    %cfg.parameter_selection.format.separator = ' ';
    %cfg.parameter_selection.optimization_criterion = 'max';

    %% This creates our crosss validation, design matrix
    % we have six runs. Each run has 12 presentations. We train on the two runs
    % and test on the other
    cfg.design = make_design_cv_overObjects(cfg);


    %% Run decoding
    results = decoding(cfg);
    save([cfg.results.dir 'cfg.mat'],'cfg');
end