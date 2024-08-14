%function analyseResults(results, num_rows, num_columns)
%function analyseResults(n_rows, n_columns, folder, rand_trials, trial_results, samp)
% gets results from experiment including data from eyetracker
% calls all the analysis/plotting functions
% colours: a        orange
%          a_simple red
%          b        blue
%          b_simple green
% analysis_folder = strcat(folder_name, '\analysis');
% mkdir(analysis_folder);
color_map = containers.Map({'a', 'a_simple', 'b', 'b_simple'}, {
    [0.9, 0.5, 0],  % orange
    [1.0, 0  , 0],  % red
    [0  , 0  , 1],  % blue
    [0  , 0.5, 0]   % green
});
%% for testing hardcoded
n_rows = 9;
n_columns = 12;
screenXpixels = 3240;
screenYpixels = 2160;


%folder = 'C:\Users\flohrmann\Documents\Results\1_20240714_140048'; %alex
%folder = 'C:\Users\flohrmann\Documents\Results\2_20240726_191755'; % mara
%folder = 'C:\Users\flohrmann\Documents\Results\3_20240805_105213'; % tilo 
%folder = 'C:\Users\flohrmann\Documents\Results\4_20240811_131601'; % anu 
folder = 'C:\Users\flohrmann\Documents\Results\5_20240813_114700'; % kieran

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

%% 
analysis_folder = strcat(folder, '\analysis');
try
    mkdir(analysis_folder);
catch % folder already exists
end

try % try get cut data
    file_pattern = fullfile(analysis_folder, 'cut_trials_*');
    file_info = dir(file_pattern);
    load(strcat(analysis_folder, '\', file_info.name)); % cut eyetracking; cut_data
catch % cut data if not already cut
    cutData = cutEyeTrackingData(analysis_folder, trial_results, samp); % cut eyetracking; cut_data
end

%% plot eye gaze and stimulation per trial
try % load data (takes forever to calc/plot, dont wanna do this twice)
    load(strcat(analysis_folder, '\eye_rt.mat')); % eye_rt
    %show = true; % show plot (very slow)?
    %num_plots = size(cutData, 1); % how many trials you want plotted, starts with first
    %eye_rt = plotStimAndEye(analysis_folder, cutData, num_plots, show);
catch % calculate/plot if first time
     %load(strcat(analysis_folder, '\eye_rt.mat')); % eye_rt
     show = false; % show plot (very slow)?
     num_plots = size(cutData, 1); % how many trials you want plotted, starts with first
     eye_rt = plotStimAndEye(analysis_folder, cutData, num_plots, show);
end

%% plot  button press - gaze rt (speed of pressing button once stim found)
plotButtonPressMinusGazeRT(trial_results, eye_rt, analysis_folder);

%% diff if stim in inner vs outer circle of screen
plotRTbyEccentricity(trial_results, eye_rt, screenXpixels, screenYpixels, n_rows, n_columns, analysis_folder);

%% distance to target at onset vs rt (gaze) to target (did they even look at the stim)
plotDistanceVsRT(trial_results, samp, eye_rt, screenXpixels, screenYpixels, analysis_folder, color_map)

%% rt vs pupil size scatterplot per condition
plotRTvsPupilSizePerCondition(trial_results, samp, analysis_folder, color_map)

%% diff button press rt vs gaze rt
plotButtonPressVsGazeRT(trial_results, eye_rt, analysis_folder)

%% reaction time numbers and plots
%% todo check for differences between rt plots/if the same/debug
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

%% did they look at fixation?
fixationThreshold = 150; % threshold for how close the gaze needs to be to the fixation cross
lookedAtFixation = checkFixation(trial_results, samp, screenXpixels, screenYpixels, fixationThreshold);
fixationSummary = any(lookedAtFixation, 2);
fixationSummary = double(fixationSummary); % Convert logical array to double for 1 and 0 output
numOnes = sum(fixationSummary(:)); 
percentageOnes = numOnes / size(lookedAtFixation,1) * 100;

%% rt trial with looking at fixation vs without
plotRTwithAndWithoutFixation(trial_results, eye_rt, lookedAtFixation, analysis_folder);

%% pupil size
plotPupilDiameterOverTime(cutData, samp, trial_results, analysis_folder)
plotPupilDiameterAverageOverTrials(cutData, analysis_folder)

%% todo:
%plotPupilSizeAroundEvents(cutData)


%% todo: diff rt over under 300 ms topdown/bottom up




%% todo: questionaire answers (look for cutoff)





