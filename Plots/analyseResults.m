%function analyseResults(results, num_rows, num_columns)
%function analyseResults(n_rows, n_columns, folder, rand_trials, trial_results, samp)
% gets results from experiment including data from eyetracker
% calls all the analysis/plotting functions
% analysis_folder = strcat(folder_name, '\analysis');
% mkdir(analysis_folder);
color_map = containers.Map({'a', 'a_simple', 'b', 'b_simple'}, {
    [0.9, 0.5, 0],  % orange
    [1.0, 0  , 0],  % red
    [0  , 0  , 1],  % blue
    [0  , 0.5, 0]   % green
    });
% slightly transparent colours
alpha = 0.2;
color_map_trans = containers.Map({'a', 'a_simple', 'b', 'b_simple'}, {
    [0.9 0.5 0 alpha],  % orange
    [1.0 0   0 alpha],  % red
    [0   0   1 alpha],  % blue
    [0   0.5 0 alpha]   % green
    });

conditions = {'a', 'a_simple', 'b', 'b_simple'};


%% for testing hardcoded
n_rows = 9;
n_columns = 12;
screenXpixels = 3240;
screenYpixels = 2160;



% Define folders and IDs for each participant
subfolders = { ...
    '1_20240714_140048', ... % alex
    '2_20240726_191755', ... % mara
    '3_20240805_105213', ... % tilo
    '4_20240811_131601', ... % anu
    '5_20240813_114700', ... % kieran
    '6_20240821_191408', ... % sura
    '7_20240821_194651', ... % hamit
    '7_20240823_162058', ... % jannik
    '9_20240829_101613', ... % farn
    '10_20240929_151434' ... % julia
    };

ids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]; 


for i = 1:length(subfolders)
    name = subfolders{i}; 
    folder = strcat('C:\Users\flohrmann\Documents\Results\', name);
    id = ids(i);         
    % Display the folder and ID being processed
    disp(['Processing folder: ', folder, ' with ID: ', num2str(id)]);
    
    
    load(strcat(folder, '\rand_trials.mat')); % load trial infos; rand_trials
    % get results file
    file_pattern = fullfile(folder, '\results', 'trial_results_*');
    file_info = dir(file_pattern);
    load(strcat(folder, '\results\', file_info.name)); % load results; trial_results
    % get eye tracking data
    file_pattern = fullfile(folder, 'results', 'eyetracking_results*');
    file_info = dir(file_pattern);
    load(strcat(folder, '\results\', file_info.name)); % eyetrackgin data; samp
    
    % make folder
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
    
    %% plot spatial spread of stims over trials and possible positions
    plotConditionSpreadAndStimPosition(rand_trials, n_rows, n_columns, analysis_folder);
    %plotTargetPositionHeatmap(rand_trials, num_rows, num_columns)
    
    
    %% plot eye gaze and stimulation per trial
    try % load data (takes forever to calc/plot, dont wanna do this twice)
        %a = notaFunction(); % fail try block
        load(strcat(analysis_folder, '\eye_rt.mat')); % eye_rt
        %show = true; % show plot (very slow)?
        %num_plots = size(cutData, 1); % how many trials you want plotted, starts with first
        %eye_rt = plotStimAndEye(analysis_folder, cutData, num_plots, show);
    catch % calculate/plot if first time
        %load(strcat(analysis_folder, '\eye_rt.mat')); % eye_rt
        %show = false; % show plot (very slow)?
        show = true; % show plot (very slow)?
        num_plots = size(cutData, 1); % how many trials you want plotted, starts with first
        eye_rt = plotStimAndEye(analysis_folder, cutData, 200, show);
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
    
    
    %% did they look at fixation?
    fixationThreshold = 200; % threshold for how close the gaze needs to be to the fixation cross
    lookedAtFixation = checkFixation(trial_results, samp, screenXpixels, screenYpixels, fixationThreshold);
    fixation_summary = any(lookedAtFixation, 2);
    
    fixationSummary = table();
    fixationSummary.fix_count = sum(fixation_summary(:,1));
    fixationSummary.fix_perc  = fixationSummary.fix_count / size(lookedAtFixation,1) * 100;
    fixationSummary.fix_l     = {double(lookedAtFixation(:, 2))};
    fixationSummary.fix_r     = {double(lookedAtFixation(:, 1))};
    fixationSummary.fix_any   = {double(fixation_summary(:))};
    
    save(fullfile(analysis_folder, 'fixationSummary.mat'), 'fixationSummary');
    
    
    %% rt trial with looking at fixation vs without
    % TODO Missing legend, Labels, Axis
    plotRTwithAndWithoutFixation(id, trial_results, eye_rt, lookedAtFixation, analysis_folder);
    
    
    
    %% --- pupil diameter ---
    % all trials appended to one another
    plotPupilDiameterOverTime(id, cutData, samp, trial_results, analysis_folder);
    
    %% data split into fixation screen, blank screen and stimulation screen
    % interpolated (mistakes in interpolation??)
    %plotPupilDiameterAverageOverTrials(id, cutData, analysis_folder);
    
    % interpolate properly
    plotNormalizedPupilDiameterAverageOverTrials(id, cutData, conditions, analysis_folder);
    % 4 tiles per condition; movemean for smoothing after interpolation; average both eyes (should have same dilation, otherwise artefact)
    plotAvgNormalizedPupilDiameterByCondition(id, cutData, conditions, analysis_folder);
    
    %% only stimulation data (until button press)
    % all conditions in one plot
    avg_results = plotAvgNormalizedPupilDiameterByConditionInOne(id, cutData, conditions, analysis_folder, color_map);
    
    % average number of datapoints during stim presentation per condition for this person
    % safes avg_num_stim_data_points_per_condition.mat
    avg_stim_data_points = calculateAverageStimDataPoints(cutData, conditions, analysis_folder);
    
    % avg is at least 30 datapoints -> 30*16ms (60hz sampling rate) = 480 ms ~ 0.5 s
    % -> only take 30 datapoints before stim found (button pressed) and average this!
    % pad with NaNs at beginning if less datapoints to ensure correctness
    avg_stim_30 = calcAvgStimFound30DataPoints(cutData, conditions, 30);
    
    % 30 points before trial end/ button pressed per condition
    plotStimFound30Points(avg_stim_30, id, analysis_folder, color_map_trans); % all
    % mean and median with SE and interquartile ranges
    plotAvgStimFound30Points(avg_stim_30, conditions, id, analysis_folder, color_map);
    
    % per condition; avg pupil diameter across trials; fit linear regression line to avg data points -> slope
    % slope represents how much the pupil diameter changes per unit time before the stimulus is found
    slope_table = analyzeSlopesByCondition(avg_stim_30, conditions, id, analysis_folder, color_map);
    
    %% only stimulation data (until target found with gaze)
    % check target coords/eye coords (+ tolerance)
    % get 30 datapoints before target is reached with eyes + 10 datapoints after
    % padded with NaNs if not enough datapoints
    tolerance = 100;  % Tolerance in pixels for detecting gaze near the target
    num_before = 30; % Number of datapoints before target found
    num_after = 10;  % Number of datapoints after target found
    
    [diam_around_stim_table, tnf] = findTargetAndExtractData(cutData, screenXpixels, screenYpixels, id, analysis_folder, tolerance, num_before, num_after);
    
    % plot 30 datapoints before target was reached with gaze and 10 after, average per condition
    plotAvgBeforeAfterStimFound(diam_around_stim_table, tnf, conditions, id, analysis_folder, color_map, num_before, num_after)
    
    
    
    %% todo: diff rt over under 300 ms topdown/bottom up
    
    
    close all
    
end

