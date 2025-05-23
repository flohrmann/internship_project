function [rand_trials, trial_results, samp, cutData] = loadData(folder, subfolder_name)
%% load data
load(strcat(folder, '\rand_trials.mat')); % load trial infos; rand_trials

file_pattern = fullfile(folder, 'results', 'trial_results_*');
file_info = dir(file_pattern); % trial results/behavioural data
load(strcat(folder, '\results\', file_info.name)); % load results; trial_results

% if fixation was used interactively and theres snippets use the continuous data
try 
    file_pattern = fullfile(folder, 'results', 'continuous_gaze*');
    file_info = dir(file_pattern); % eyetracking data
    load(strcat(folder, '\results\', file_info.name)); % eyetrackgin data; samp
    samp = continuous_gaze; % "rename" so its the same name as in the old data
catch % if not data is continuous in this place anyways
    file_pattern = fullfile(folder, 'results', 'eyetracking_results*');
    file_info = dir(file_pattern); % eyetracking data
    load(strcat(folder, '\results\', file_info.name)); % eyetrackgin data; samp
end
%% make folder for analysis plots and files
analysis_folder = strcat(folder, subfolder_name);
try
    mkdir(analysis_folder);
catch % folder already exists
end

%% try to load cut data (already cut in trials)
try 
    file_pattern = fullfile(folder, subfolder_name, 'cut_trials*');
    file_info = dir(file_pattern);
    load(strcat(folder, subfolder_name, '\', file_info.name)); % cut eyetracking; cut_data
catch % cut data if not already cut
    cutData = cutEyeTrackingData(strcat(folder, subfolder_name), trial_results, samp); % cut eyetracking; cut_data
end