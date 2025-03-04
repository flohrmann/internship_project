%% Load all the varibles and get all the paths
%function analyseResults(results, num_rows, num_columns)
%function analyseResults(n_rows, n_columns, folder, rand_trials, trial_results, samp)
% gets results from experiment including data from eyetracker
% calls all the analysis/plotting functions
% analysis_folder = strcat(folder_name, '\analysis');
% mkdir(analysis_folder);

color_map = containers.Map({'a', 'a_simple', 'b', 'b_simple'}, {
    [0.9, 0.5, 0]  % orange
    [1.0, 0  , 0]  % red
    [0  , 0.5, 0]  % green
    [0  , 0  , 1]  % blue
    });
% slightly transparent colours
alpha = 0.2;
color_map_trans = containers.Map({'a', 'a_simple', 'b', 'b_simple'}, {
    [0.9 0.5 0 alpha]  % orange
    [1.0 0   0 alpha]  % red
    [0   0.5 0 alpha]  % green
    [0   0   1 alpha]  % blue
    });

conditions = {'a', 'a_simple', 'b', 'b_simple'};
condition_labels = {'a', 'a simple', 'b', 'b simple'};

% for testing hardcoded
n_rows = 9;
n_columns = 12;
screenXpixels = 3240;
screenYpixels = 2160;
sr = 1/60;              % Sampling rate of eyetracker
safe = 1; % safe plots
fullscreen = 1; % make plots fullscreen
do_plots = 3;
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
    '11_20241006_153704' ... % felix
    '12_20241006_150321' ... % florian
    '13_20241025_195047' ... % leandra
    '14_20241028_084922' ... % finn
    '15_20241114_200512' ... % david
    '16_20241128_113348' ... % bennet
    '17_20250126_202511' ... % kerstin
    '18_20250222_153229' ... % stacey
    '19_20250224_120844' ... % maria
    '20_20250224_171837' ... % lea
    '21_20250224_174504' ... % giove
    };

ids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21];
compare_folder = 'C:\Users\flohrmann\Documents\Analysis_2025_final'; % safe some plots here for easier comparison
analysis_subfolder = '\analysis_final';

%%

for subject = 15: length(subfolders)
    name = subfolders{subject};
    folder = strcat('C:\Users\flohrmann\Documents\Results\', name);
    analysis_folder = strcat(folder, analysis_subfolder);
    id = ids(subject);
    disp(['Processing folder: ', folder, ' with ID: ', num2str(id)]);
    
    % load data & cut into trials
    [rand_trials, trial_results, samp, cut_data] = loadData(folder, analysis_subfolder);
    
    %% --- Analysis ---
    % Plot gaze and stimulation per trial (only uses data recorded during stimulation screen)
    % eyeTrial: starts with fixation onset
    % stimulusTrial: starts with Stimulation onset
    try % load data (takes forever to calc/plot, dont wanna do this twice)
        %a = notaFunction(); % fail try block
        load(strcat(analysis_folder, '\eye_rt.mat')); % eye_rt
    catch % calculate/plot if first time
        show = false; % dont show plots (slightly faster)
        num_plots = 10; % size(cut_data, 1); % how many trials you want plotted, starts with first
        eye_rt = plotStimAndEye(analysis_folder, cut_data, num_plots, show, 'plotsanddata');
    end
    
    %Behavioural Data (RT, accuracy, RTV, confusion)
    % Reaction time distribution per condition
    plotRTDistributionPerCondition(trial_results, analysis_folder, condition_labels, color_map);

    % Plot the speed of pressing button once the target was found (button press - gaze rt)
    plotButtonPressMinusGazeRT(trial_results, eye_rt, analysis_folder);
    
    %Plot the RT per trial and violin plot/average/std error of mean RT per condition
    plotRTOverTimeColouredByCondition(trial_results, color_map, analysis_folder, condition_labels)
    plotRT(trial_results, analysis_folder, condition_labels)

    %Accuracy per condition
    plotAccuracyPerCondition(trial_results, color_map, safe, analysis_folder)
    
    % rt/accuracy
    plotRTvsAccuracy(trial_results, analysis_folder)
    
    %Reaction Time Variability
    rt_variability = getRTvariabilityPerParticipant(trial_results, conditions, color_map, analysis_folder);
    
    % diff if stim in inner vs outer circle of screen
    % TODO check visul angle 
    plotRTbyEccentricity(trial_results, eye_rt, screenXpixels, screenYpixels, n_rows, n_columns, analysis_folder);
    
    % distance to target at onset vs rt (gaze) to target (did they even look at the stim)
    plotDistanceVsRT(trial_results, samp, eye_rt, screenXpixels, screenYpixels, analysis_folder, color_map)
    
    % rt vs pupil size scatterplot per condition
    plotRTvsPupilSizePerCondition(trial_results, samp, analysis_folder, color_map)
    
    % did they look at fixation?
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
    
    % rt trial with looking at fixation vs without
    plotRTwithAndWithoutFixation(id, trial_results, eye_rt, lookedAtFixation, analysis_folder, condition_labels);
    
    % Spatial spread of stims over trials and possible positions
    plotConditionSpreadAndStimPosition(rand_trials, n_rows, n_columns, analysis_folder, condition_labels);
    
    
    %% --- Eyetracking Data  ---
    % plotPupilDiameterOverTime(id, cut_data, analysis_folder);
    plot_these = [1, 2, 5, 20, 100]; % just plot some trials for checking
    fixation_durations = [50, 100, 200]; % in ms
    
    for fd_idx=1:(size(fixation_durations,2))
        fd_label = strcat(num2str(fixation_durations(fd_idx)), 'ms');
        analysis_folder = strcat(folder, '\analysis_fixation_', fd_label);
        compare_folder = strcat('C:\Users\flohrmann\Documents\Analysis_2025\PupilDiam_fix_', fd_label); % safe some plots here for easier comparison
        
        sampling_rate = sr *1000; % in seconds
        min_duration = round(fixation_durations(fd_idx)/sampling_rate); % min num of consecutive datapoints to count as a fixation
        
        try
            mkdir(analysis_folder);
            mkdir(compare_folder);
        catch % folders already exists thats fine
        end
        
        % - 1. Detect fixations -
        dist_threshold = 50; % max distance between consecutive points (in pixels) to count as fixation cluster
        fixations = analyseFixation(id, plot_these, cut_data, dist_threshold, min_duration, screenXpixels, screenYpixels, safe, analysis_folder, fd_label);
        
        % - 2. Reverse engineer "saccades" from fixation clusters -
        saccades = reverseEngineerSaccades(cut_data, fixations, screenXpixels, screenYpixels, plot_these, analysis_folder, safe);
        
        % - 3. Get num saccades, distance, angle, cosine from trial start to target center/saccades,... -
        % fill in baselines with not enough data with NaNs
        % find saccades towards target AND find all saccades NOT towards target
        % remove saccades with less than half of the datapoints
        % cut saccades towards target (TS): aligned by end of saccade (TSE) and start of saccade (TSS)
        % cut saccades NOT towards target (NTS): aligned by end of saccade (NTSE) and start of saccade (NTSS)
        baseline_length = 15;   % Length for baseline used for normalizing later; 250 ms
        num_before = 20;        % Number of datapoints before target found
        num_after = 31;         % Number of valid pupil diam datapoints after target found 0.5/0.016=31.3
        min_length = 10;        % number of valid datapoints needed for data about 1/5
        tolerance = 200;        % num pixels away to count as target found
        trial_metrics = calcGazeALLSaccadesDistanceDirection(cut_data, fixations, saccades, tolerance, num_before, num_after, baseline_length, screenXpixels, screenYpixels, min_length);
        save(fullfile(analysis_folder, strcat('\trial_metrics.mat')), 'trial_metrics');
                
        %  - 4. Look at changes in pupil dilation -      
        % Reverse engineer Saccades from Fixations and take as t0:
        % 1. [Start of Saccade that Reached Target]
        % 2. [End of Saccade that Reached Target]
        % 3. [Start of Saccade that didnt Reach Target]
        % 4. [End of Saccade that didnt Reach Target]
        % 5. [FIRST Start Target Saccade of trial]
        % 6. [FIRST End of Target Saccade of trial]
        % 7. [LAST Start of Target Saccades per trial] (uneccessary bc fixation as marker for fixation ensures the saccade was goal oriented)
        time_vector = (-num_before:num_after - 1) * sr; % Time in seconds
        x_label_text = 'Time from t0 (s)';
        analyseChangePupilDilation(id, cut_data, trial_metrics, conditions, condition_labels, analysis_folder, compare_folder, do_plots, time_vector, x_label_text, color_map);
        close all; % close plots
        
        % search efficiency
        % average number of saccades until target found per condition + sem
        % average number of saccades after target found per condition + sem
        % average distance of saccade before saccade towards target and average distance of saccade towards target per condition (corective saccades)
        % average difference between optimal angle and angle of first saccade per condition
        % --- overall sacacdes ---
        % per condition and once for all trials
        plotHistSaccadeMetricsParticipant(trial_metrics, cut_data, conditions,color_map, sr, ...
                                          fullfile(compare_folder, strcat('saccade_stats_hist_all_id_',num2str(id))), ...
                                          fullfile(compare_folder, strcat('saccade_stats_hist_condition_id_',num2str(id))))
        plotSaccadeMetricsViolin(trial_metrics, cut_data, conditions, color_map, condition_labels, fullfile(compare_folder, ...
                                strcat('saccade_num_dist_angle_violin_id_',num2str(id),'.svg')));
         
        % --- first sacacde ---
        % - reaction time/ start of first sacacde -
        first_sacc = plotAverageStartFirstSaccade(id, trial_metrics, trial_results, conditions, condition_labels, color_map, compare_folder, analysis_folder);
        
        % difference of first saccade angle & distance compared to optimal one
        % pretty but useless plot
        %plotDiffFirstSaccadeAngleDistance(trial_metrics, compare_folder, safe)
        
        % cosine similarity
        plotCosineSimilarity(trial_results, trial_metrics, saccades, fixations, screenXpixels, screenYpixels, ...
            conditions, condition_labels, color_map, safe, analysis_folder)
        close all
    end
end


