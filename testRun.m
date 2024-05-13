%% Load custom path for testing experiment
data_struct = load('C:\Users\idm\Desktop\Semester4\Internship\Matlab\ExperimentRepo\test\rand_trials.mat');
data = data_struct.rand_trials;
folder_name = 'C:\Users\idm\Desktop\Semester4\Internship\Matlab\ExperimentRepo\test';
results = startPsychToolboxEyeTracking(data, folder_name, 9, 12, 10, 2);

%% Fix fixation data
data_struct = load('C:\Users\idm\Desktop\Semester4\Internship\Matlab\ExperimentRepo\test\results\trial_results_with_fix.mat');
trial_results = data_struct.trial_results;
folder_name = 'C:\Users\idm\Desktop\Semester4\Internship\Matlab\ExperimentRepo\test\results';

clean_data = aggregateFixationData(trial_results);
trial_results.sampFixAll = table2array(clean_data);
%results_file_name = [folder_name, '\trial_results_fixed.mat'];
%save(results_file_name, 'trial_results');
plot_this(trial_results)
