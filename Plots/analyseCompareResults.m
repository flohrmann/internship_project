%% load all data & define variables
% save path
comparison_results_folder = 'C:\Users\flohrmann\Documents\Analysis\';
% save plots yes/no 1/0 
safe = 0;
% colour scheme
color_map = containers.Map({'a', 'a_simple', 'b', 'b_simple', 'ADHD', 'nonADHD'}, {
    [0.9, 0.5, 0  ]  % orange
    [1.0, 0  , 0  ]  % red
    [0  , 0  , 1  ]  % blue
    [0  , 0.5, 0.1]  % green
    [0.5, 0.7, 0.5]  % orange ish
    [0.7, 0.2, 0.2]  % green ish
    
    });

folders = {
    'C:\Users\flohrmann\Documents\Results\1_20240714_140048'; % alex
    'C:\Users\flohrmann\Documents\Results\2_20240726_191755'; % mara adhd
    'C:\Users\flohrmann\Documents\Results\3_20240805_105213'; % tilo adhd
    'C:\Users\flohrmann\Documents\Results\4_20240811_131601'; % anu adhd
    'C:\Users\flohrmann\Documents\Results\5_20240813_114700'; % kieran adhd
    };

% Define the classification of the folders (ADHD or non-ADHD)
group_labels = {
    'nonADHD';  % alex
    'ADHD';     % mara adhd
    'ADHD';     % tilo adhd
    'ADHD';     % anu adhd
    'ADHD';     % kieran adhd
    };

% Initialize a struct to store data for each person
data_struct = struct();

% Loop through each folder and load the data
for i = 1:length(folders)
    folder = folders{i};
    %folder = strcat(folders{i}, '\results');
    group = group_labels{i};
    
    % Load the data (assuming data is stored in a .mat file named 'trial_results.mat')
    data_file   = dir(fullfile(folder, '\results\trial_results*.mat')); % original data
    eye_rt_file = dir(strcat(folder, '\analysis\eye_rt.mat')); % calculated eye rt from analyseResults.m
    if ~isempty(data_file)
        load(fullfile(data_file.folder, data_file.name), 'trial_results');
        load(fullfile(eye_rt_file.folder, eye_rt_file.name), 'eye_rt');
        % Store the data in the struct
        data_struct(i).id                = i;
        data_struct(i).group             = group;
        data_struct(i).Condition         = trial_results.Condition;
        data_struct(i).rt                = trial_results.rt;
        data_struct(i).accuracy          = trial_results.correct;
        data_struct(i).StimulusOnsetTime = double(trial_results.StimulusOnsetTime/ 1e6);
        data_struct(i).trialEndTime      = double(trial_results.trialEndTime/ 1e6);
        data_struct(i).rt_right          = eye_rt.RightEyeRT;
        data_struct(i).rt_left           = eye_rt.LeftEyeRT;
    else
        warning(['Data file not found in folder: ', folder]);
    end
end

%% normalized RT per group and condition
% Define the simple conditions to use for normalization
condition_A = 'a_simple';
condition_B = 'b_simple';

% get the faster of the two eye RTs; ignore 0s where target was not looked at
data_struct = processEyeRTs(data_struct);
% normalize RTs
data_struct_norm = normalizeRTsBySimpleConditions(data_struct, condition_A, condition_B);

% set data to normalized data
data = data_struct_norm;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% violin plot rt time across groups and conditions
plotViolinRTGroupCondition(data, color_map, comparison_results_folder, safe)


%% normalized RT across groups
plotNormalizedRTCcomparison(data, color_map, comparison_results_folder, safe);


%% accuracy comparison mean points 
plotMeanAccuracyPerGroupCondition(data_struct, color_map, safe, comparison_results_folder)


%% mean accuracy comparison group and condition + mean points 
%plotBarAccuracyGroupCondition(data, color_map, comparison_results_folder, safe)
% TODO: change bars to violin!

%% rt button press vs rt eyes; scatter plot
plotRTButtonPressVsEye(data, comparison_results_folder, safe)


%% rt button press vs rt eyes; means of groups per condition
plotMeanRTButtonPressVsEyecomparison(data, color_map, comparison_results_folder, safe)










