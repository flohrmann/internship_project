%function analyseResults(results, num_rows, num_columns)
%% for testing hardcoded
num_rows = 9;
num_columns = 12;

%% plot verteilung stims
%old
%load('C:\Users\flohrmann\Documents\Results\fani_test\rand_trials.mat'); % rand_trials
%new
load('C:\Users\flohrmann\Documents\Results\fani_1_20240701_135200\rand_trials.mat') % rand_trials
plotConditionSpreadAndStimPosition(rand_trials, num_rows, num_columns)
%plotTargetPositionHeatmap(rand_trials, num_rows, num_columns)

%% plot eye gaze and trials
plotStimAndEye(trial_results)

%% get latency exetracker/monitor
load('C:\Users\flohrmann\Documents\Results\fani_test\results\trial_results_20240701_1123.mat'); %trial_results
plotRT(trial_results)
plotRTDistributionPerCondition(trial_results)
plotAccuracyPerCondition(trial_results)
plotRTvsAccuracy(trial_results)

%end