function mainAnalysis()

% colour scheme
color_map = containers.Map({'a',  'b', 'a_simple', 'b_simple', 'ADHD', 'nonADHD'}, {
    [0.9, 0.5, 0  ]  % orange
    [0  , 0.5, 0.1]  % green
    [1.0, 0  , 0  ]  % red
    [0  , 0  , 1  ]  % blue
    [0.5, 0.7, 0.5]  % green ish
    [0.7, 0.2, 0.2]  % orange ish
    });

% slightly transparent colours
alpha = 0.2;
color_map_trans = containers.Map({'a', 'a_simple', 'b', 'b_simple'}, {
    [0.9 0.5 0 alpha]  % orange
    [1.0 0   0 alpha]  % red
    [0   0.5 0 alpha]  % green
    [0   0   1 alpha]  % blue
    });
% colors for non condition/group plots
color_map_other = containers.Map({'1',  '2', '3', '4', '5', '6'}, {
    [0.878, 0.349, 0.682]  % pink
    [0.259, 0.651, 0.533]  % green turquis
    [0.506, 0.212, 0.6]    % purple
    [0.259, 0.427, 0.651]  % blueish
    [0.506, 0.112, 0.5]    % 
    [0.559, 0.327, 0.55]   % 
    });
conditions = {'a', 'b','a_simple', 'b_simple'};
condition_labels = {'a', 'b','a simple', 'b simple'};
groups = {'ADHD', 'nonADHD'};

n_rows = 9;
n_columns = 12;
screenXpixels = 3240;
screenYpixels = 2160;
sr = 1/60; % sampling rate of eyetracker
safe = 1; % safe plots
do_plots = 1;

fullscreen = 1; % make plots fullscreen

% folders and IDs for each participant
results_path = 'C:\Users\idm\Desktop\Semester4\Internship\Results_2025_03\';
subfolders = { ...
    '1_20240714_140048', ... % only 160 trials
    '2_20240726_191755', ... 
    '3_20240805_105213', ... 
    '4_20240811_131601', ... 
    '5_20240813_114700', ... 
    '6_20240821_191408', ... 
    '7_20240821_194651', ...
    '7_20240823_162058', ... 
    '9_20240829_101613', ... 
    '10_20240929_151434' ... 
    '11_20241006_153704' ... 
    '12_20241006_150321' ... 
    '13_20241025_195047' ... 
    '14_20241028_084922' ... 
    '15_20241114_200512' ... 
    '16_20241128_113348' ... 
    '17_20250126_202511' ... 
    '18_20250222_153229' ... % chance performance
    '19_20250224_120844' ... 
    '20_20250224_171837' ... % took medication 
    '21_20250224_174504' ... 
    };

ids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21];
comparison_results_folder = 'C:\Users\idm\Desktop\Semester4\Internship\analysis_group_2025_03\'; % safe some plots here for easier comparison
analysis_subfolder = '\analysis'; %'\analysis_final';


%% individual results - needs to run first
% takes long to plot each individual trials gaze path
 fixation_threshold = 200; % threshold for how close the gaze needs to be to the fixation cross/target to count as it being fixated
 fix_cluster_threshold = 50; % max distance between consecutive points (in pixels) to count as fixation cluster
 analyseResults(color_map, color_map_trans, conditions, condition_labels, ...
                n_rows, n_columns, screenXpixels, screenYpixels, sr, safe, do_plots, fullscreen, ...
                results_path, subfolders, ids, comparison_results_folder, analysis_subfolder,...
                fixation_threshold, fix_cluster_threshold)

%% group results
% save path
%analysis_subfolder = '\analysis'; % subfolders in participants where fixation independet results were saved earlier
fixation_duration = 100; % 50,100,200 
subfolder_fixation = strcat('\analysis_fixation_', num2str(fixation_duration),'ms');
try
    mkdir(comparison_results_folder);
catch % folders already exists
end

analyseCompareResults(color_map, color_map_other, conditions, condition_labels, groups, ...
                      screenXpixels, screenYpixels, sr, safe, ...
                      results_path, subfolders(2:end), ids(2:end), analysis_subfolder, comparison_results_folder, ...
                      subfolder_fixation, fixation_duration);




%% infos:
% subject gaming was accidentally safed as false instead of in hours
% 
%