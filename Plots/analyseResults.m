%function analyseResults(results, num_rows, num_columns)
%% for testing hardcoded
num_rows = 9;
num_columns = 12;
folder = 'C:\Users\flohrmann\Documents\Results\1_20240714_140048';

% path for saving
analysis_folder = strcat(folder, '\analysis');

load(strcat(folder, '\rand_trials.mat')); % load trial infos; rand_trials

% get results file
file_pattern = fullfile(folder, '\results', 'trial_results_*');
file_info = dir(file_pattern);
load(strcat(folder, '\results\', file_info.name)); % load results; trial_results

% get eye tracking data
file_pattern = fullfile(folder, 'results', 'eyetracking_results*');
file_info = dir(file_pattern);
load(strcat(folder, '\results\', file_info.name)); % eyetrackgin data; samp
% get cut eye tracking data
try % try cut data
    file_pattern = fullfile(analysis_folder, 'cut_trials_*');
    file_info = dir(file_pattern);
    load(strcat(analysis_folder, '\', file_info.name)); % cut eyetracking; cut_data
catch % cut data if not already
    cutData = cutEyeTrackingData(analysis_folder, trial_results, samp); % cut eyetracking; cut_data
end

%% reaction time numbers and plots
rt_per_condition = rtTimes(cutData, analysis_folder);
save(fullfile(analysis_folder, 'rt_per_condition.mat'), 'rt_per_condition');
plotRT(trial_results, analysis_folder)
plotRTDistributionPerCondition(trial_results, analysis_folder)

%% rt/accuracy
plotAccuracyPerCondition(trial_results, analysis_folder)
plotRTvsAccuracy(trial_results, analysis_folder)

%% plot verteilung stims
plotConditionSpreadAndStimPosition(rand_trials, num_rows, num_columns, analysis_folder)
%plotTargetPositionHeatmap(rand_trials, num_rows, num_columns)

%% plot eye gaze and stimulation per trial
num_plots = 10; % how many trials you want plotted, starts with first
plotStimAndEye(analysis_folder, cutData, num_plots)

%% pupil size
plotPupilDiameterOverTime(cutData, samp, trial_results, analysis_folder)
plotPupilDiameterAverageOverTrials(cutData, analysis_folder)
plotPupilSizeAroundEvents(cutData)
%end