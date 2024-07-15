%function analyseResults(results, num_rows, num_columns)
function analyseResults(n_rows, n_columns, folder, rand_trials, trial_results, samp)
% gets results from experiment including data from eyetracker
% calls all the analysis/plotting functions
analysis_folder = strcat(folder, '\analysis');

%% for testing hardcoded
% % n_rows = 9;
% % n_columns = 12;
% %folder = 'C:\Users\flohrmann\Documents\Results\1_20240714_140048';
% load(strcat(folder, '\rand_trials.mat')); % load trial infos; rand_trials
% % get results file
% file_pattern = fullfile(folder, '\results', 'trial_results_*');
% file_info = dir(file_pattern);
% load(strcat(folder, '\results\', file_info.name)); % load results; trial_results
% % get eye tracking data
% file_pattern = fullfile(folder, 'results', 'eyetracking_results*');
% file_info = dir(file_pattern);
% load(strcat(folder, '\results\', file_info.name)); % eyetrackgin data; samp
% % get cut eye tracking data
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
plotConditionSpreadAndStimPosition(rand_trials, n_rows, n_columns, analysis_folder)
%plotTargetPositionHeatmap(rand_trials, num_rows, num_columns)

%% plot eye gaze and stimulation per trial
show = false;
num_plots = size(cutData, 1); % how many trials you want plotted, starts with first
eye_rt = plotStimAndEye(analysis_folder, cutData, num_plots, show);

%% did they look at fixation?
screenXpixels = 3240;
screenYpixels = 2160;
fixationThreshold = 150; % threshold for how close the gaze needs to be to the fixation cross
lookedAtFixation = checkFixation(trial_results, samp, screenXpixels, screenYpixels, fixationThreshold);
fixationSummary = any(lookedAtFixation, 2);
fixationSummary = double(fixationSummary); % Convert logical array to double for 1 and 0 output
numOnes = sum(fixationSummary(:)); 
percentageOnes = numOnes / numel(lookedAtFixation) * 100;

%% pupil size
plotPupilDiameterOverTime(cutData, samp, trial_results, analysis_folder)
plotPupilDiameterAverageOverTrials(cutData, analysis_folder)
plotPupilSizeAroundEvents(cutData)
%end
%% distance to target at onset vs rt (gaze) to target

%% rt vs pupil size scatterplot per condition
