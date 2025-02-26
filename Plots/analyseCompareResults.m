%% load all data & define variables
% save path
analysis_subfolder = '\analysis_new';
subfolder_fixation = '\analysis_fixation_100ms'; fixation_duration = 100;

comparison_results_folder = 'C:\Users\flohrmann\Documents\Analysis_2025\';
%comparison_results_folder = strcat(comparison_results_folder, num2str(fixation_duration), 'ms\');
try
    mkdir(comparison_results_folder);
catch % folders already exists
end

% save plots yes/no 1/0
safe = 1;

% colour scheme
color_map = containers.Map({'a',  'b', 'a_simple', 'b_simple', 'ADHD', 'nonADHD'}, {
    [0.9, 0.5, 0  ]  % orange
    [0  , 0.5, 0.1]  % green
    [1.0, 0  , 0  ]  % red
    [0  , 0  , 1  ]  % blue
    [0.5, 0.7, 0.5]  % orange ish
    [0.7, 0.2, 0.2]  % green ish
    });

ids = [1, 2, 3, 4, 5, 6, 7, 8, 9 , 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21];
folders = {
    'C:\Users\flohrmann\Documents\Results\1_20240714_140048'; % alex
    'C:\Users\flohrmann\Documents\Results\2_20240726_191755'; % mara
    'C:\Users\flohrmann\Documents\Results\3_20240805_105213'; % tilo
    'C:\Users\flohrmann\Documents\Results\4_20240811_131601'; % anu
    'C:\Users\flohrmann\Documents\Results\5_20240813_114700'; % kieran
    'C:\Users\flohrmann\Documents\Results\6_20240821_191408'; % sura
    'C:\Users\flohrmann\Documents\Results\7_20240821_194651'; % hamit
    'C:\Users\flohrmann\Documents\Results\7_20240823_162058'; % jannik
    'C:\Users\flohrmann\Documents\Results\9_20240829_101613'; % farn
    'C:\Users\flohrmann\Documents\Results\10_20240929_151434'; % julia
    'C:\Users\flohrmann\Documents\Results\11_20241006_153704'; % felix
    'C:\Users\flohrmann\Documents\Results\12_20241006_150321'; % florian
    'C:\Users\flohrmann\Documents\Results\13_20241025_195047'; % leandra
    'C:\Users\flohrmann\Documents\Results\14_20241028_084922'; % finn
    'C:\Users\flohrmann\Documents\Results\15_20241114_200512'; % david
    'C:\Users\flohrmann\Documents\Results\16_20241128_113348'; % bennet
    'C:\Users\flohrmann\Documents\Results\17_20250126_202511'; % kerstin
    'C:\Users\flohrmann\Documents\Results\18_20250222_153229'; % stacey
    'C:\Users\flohrmann\Documents\Results\19_20250224_120844'; % maria
    'C:\Users\flohrmann\Documents\Results\20_20250224_171837'; % lea
    'C:\Users\flohrmann\Documents\Results\21_20250224_174504'; % giove
    };

% Define the classification of the folders (ADHD or non-ADHD)
group_labels = {'nonADHD'; 'ADHD'   ; 'ADHD'   ; 'ADHD'   ; 'ADHD'   ; 
                'ADHD'   ; 'nonADHD'; 'nonADHD'; 'ADHD'   ; 'nonADHD'; 
                'ADHD'   ; 'nonADHD'; 'nonADHD'; 'ADHD'   ; 'ADHD'   ; 
                'nonADHD'; 'nonADHD'; 'nonADHD'; 'nonADHD'; 'ADHD'   ;
                'nonADHD'};

% individual symbol and colour per participant
color_map_individual = getIndividualMarkersPerParticipant(group_labels, ids, comparison_results_folder);

conditions = {'a', 'b', 'a_simple', 'b_simple'};
condition_labels = {'a', 'b', 'a simple',  'b simple'};
groups = {'ADHD', 'nonADHD'};

% experiment screen size
screenXpixels = 3240;
screenYpixels = 2160;

data_struct = struct();
% Loop through each folder and load the data
for i = 1:length(folders)
    folder = folders{i};
    %folder = strcat(folders{i}, '\results');
    group = group_labels{i};
    data_file   = dir(fullfile(folder, '\results\trial_results*.mat')); % original data
    eyetracking_file = dir(fullfile(folder, analysis_subfolder, '\cut_trials*.mat'));
    load(strcat(folder, analysis_subfolder, '\eye_rt.mat')); % calculated eye rt from analyseResults.m
    %pupil_file = dir(strcat(folder, analysis_subfolder, '\pupil_20before_31after_finding_stim.mat'));
    %pupil_diam_norm_file = dir(strcat(folder, analysis_subfolder, '\pupil_processed_found_stim_StartOfTrial.mat'));
    
    if ~isempty(data_file)
        load(fullfile(data_file.folder, data_file.name), 'trial_results');
        data_struct(i).id                = i;
        data_struct(i).group             = group;
        data_struct(i).Condition         = trial_results.Condition;
        data_struct(i).rt                = trial_results.rt;
        data_struct(i).accuracy          = trial_results.correct;
        data_struct(i).StimulusOnsetTime = trial_results.StimulusOnsetTime;
        data_struct(i).trialEndTime      = trial_results.trialEndTime;
        % eye reaction times
        %load(fullfile(eye_rt_file.folder, eye_rt_file.name), 'eye_rt');
        data_struct(i).rt_right          = eye_rt.RightEyeRT;
        data_struct(i).rt_left           = eye_rt.LeftEyeRT;
        % eyetracking data cut into trials
        try
            load(fullfile(eyetracking_file.folder, eyetracking_file.name), 'cutData');
        catch
            load(fullfile(eyetracking_file(1).folder, eyetracking_file(1).name), 'cutData');
        end
        
        data_struct(i).eyeTrial           = cutData.eyeTrial;
        data_struct(i).stimulusTrial      = cutData.stimulusTrial; % cut to start stim
        data_struct(i).TargetPosition     = cutData.TargetPosition;
        data_struct(i).x_centers          = cutData.x_centers;
        data_struct(i).y_centers          = cutData.y_centers;
        %         % pupil dilation data
        %         load(fullfile(pupil_file.folder, pupil_file.name), 'result_table');
        %         data_struct(i).pupilDiam          = result_table;
        %         load(fullfile(pupil_diam_norm_file.folder, pupil_diam_norm_file.name), 'condition_diam_avg');
        %         data_struct(i).pupilDiamNorm      = condition_diam_avg;
        
        try
            load(strcat(folder, subfolder_fixation, '\diam_t0_ntse.mat'));
            load(strcat(folder, subfolder_fixation, '\diam_t0_ntss.mat'));
            load(strcat(folder, subfolder_fixation, '\diam_t0_tse.mat'));
            load(strcat(folder, subfolder_fixation, '\diam_t0_tss.mat'));
            load(strcat(folder, subfolder_fixation, '\diam_t0_tss_first.mat'));
        catch
            disp('no data for id:', i)
        end
        load(strcat(folder, subfolder_fixation, '\trial_metrics.mat')); % TODO change path
        data_struct(i).ntse         = diam_t0_ntse;
        data_struct(i).ntss         = diam_t0_ntss;
        data_struct(i).tse          = diam_t0_tse;
        data_struct(i).tss          = diam_t0_tss;
        data_struct(i).ntse_first   = diam_t0_tss_first;
        data_struct(i).trial_metrics= trial_metrics;
        
        load(fullfile(folder, '\analysis_fixation_100ms\first_sacc.mat'));
        data_struct(i).first_sacc = first_sacc;
    else
        warning(['Data file not found in folder: ', folder]);
    end
end


% get the faster of the two eye RTs; ignore 0s where target was not looked at
data_struct = processEyeRTs(data_struct);
% normalized RT per participan, take mean or median for rt/noramlisation
average = 'median'; % 'mean' or 'median'
data_struct_norm_median = normalizeMeanRTsBySimpleConditions(data_struct, average, conditions, comparison_results_folder);

% set data to normalized data
data = data_struct_norm_median;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS/ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REACTION TIME AND ACCURACY
%% Compare RT_ButtonPress and RT_Eye across the four conditions ("a", "a simple", "b", "b simple").
%% Hyp: Easier conditions ("a simple", "b simple") might have faster RTs than harder conditions ("a", "b").
%% Hyp: ADHD participants might show different reaction time patterns compared
%%      to non-ADHD participants, especially in harder conditions.

%% Reaction Time
% 3x3 subplots: rt button press vs rt eyes vs diff; means of groups per condition
plotMeanRTButtonPressVsEyecomparison(data, color_map, conditions, comparison_results_folder, safe)

%% RT variability
measure = 'Standard Deviation';  % Use "Standard Deviation" or "Mean Squared Successive Difference" (or "Ex-Gaussian Analysis")
[rt_median_eye,rt_median_button]  = getRTvariability(data, group_labels, ids, conditions, condition_labels, measure, color_map, color_map_individual, comparison_results_folder);
%rt_variability_2 = getAndPlotRTV_ParticipantLines(data, conditions,condition_labels, 'std', groups, color_map, comparison_results_folder);

%% Confusion/ RTa - RTb or RTa/RTb for both button press/eye RT
% Mean confusion of participants for each condition plus p values in names
data = plotConfusionPerGroup(data, average, color_map, comparison_results_folder, safe);

%% Accuracy
% H: Accuracy might be lower in harder conditions for ADHD participants compared to non-ADHD participants
% Difference in correct trials
[mistakes, accuracy] = plotMeanAccuracyPerGroupCondition(data, group_labels, ids, conditions, condition_labels, color_map,color_map_individual, comparison_results_folder, safe);

                        %plotBarAccuracyGroupCondition(data, color_map, comparison_results_folder, safe)
% accuracy vs rt button press per condition and group
% anova and lme
plotAccuracyVsButtonPressRT(data, mistakes, accuracy, rt_median_eye, rt_median_button, group_labels, ids, conditions, condition_labels, color_map,color_map_individual, comparison_results_folder, safe)

% plot points per participant and connect with line
% wrong trials: simple vs nonsimple conditions per group
% TODO FIX
%plotAccuracyVsRTandErrorRates(data, color_map)



            %% correlation analysis: RTeye, RTbuttonpress and accuracy
            %[R, P] = correlationRTbuttonRTeyeAccuracy(data, comparison_results_folder);


%% condition-Specific Correlation Analysis
%[R_condition, P_condition] = correlationByCondition(data, conditions, comparison_results_folder);
% ADHD vs. non-ADHD Group Correlation Analysis
%[R_group, P_group] = correlationByGroup(data, groups, comparison_results_folder);

% correlation has nans in b_simple since theres only 1 incorrect trial for
% all participants/trials
checkAccuracyDistributionByCondition(data, conditions);


            %% BOTH TODO look up some values
            % [R_group_condition, P_group_condition] = correlationByGroupCondition(data, groups, conditions, comparison_results_folder);
            % % TODO add Rtbutton-rteye (lapse?)
            % plotGroupConditionCorrelations(R_group_condition, groups, conditions)
            % compareGroupCorrelations(R_group_condition, groups, data, conditions)


%% BAR and LINE plot: mean RT and error bars (~variability) per group and condition
%plotBarSEMMeanRTConditonGroup(groups, conditions, data, color_map, safe, comparison_results_folder)

% rt over time per participant
%% TODO change to subplots per condition
plotRTperParticipantOverTimeColouredByGroup(data, color_map, safe, comparison_results_folder)

% NOT NEEDED normalized RT across groups
%plotNormalizedRTCcomparison(data, color_map, comparison_results_folder, safe);



            %% --- parallel advantage -2022 zhaoping paper plots ---
            % Mean reaction times of participants for each condition
            plotBarRTmeanParticipantsCondition(data, color_map, comparison_results_folder);
            %same but per group
            plotBarRTmeanParticipantsConditionPerGroup(data, color_map, comparison_results_folder);

            % TODO RT for central vs. peripheral targets (per condition) TODO
            % " target that is more or less than its average distance (14 degrees) from the
            % display center is termed “peripheral” or “central,” respectively." -zp22


            %  TODO Learning effects 
            % "Each observer had several practice trials right before the 400 testing trials.
            % The saliency effect, manifested as a shorter RT for the parallel than the perpendicular target,
            % was significant for the peripheral but not the central targets in the first 200 testing trials,
            % and became significant for the central targets in the second 200 testing
            % trials."-zp22
            % Central vs. Peripheral Classification:
            %     The function calculates the Euclidean distance from the target to the center of the screen for each trial.
            %     Based on a threshold (half of the screen diagonal), targets are classified as "central" or "peripheral."
            % RT Collection by Block:
            %     The function divides trials into blocks of 100.
            %     It groups reaction times into "central" and "peripheral" categories for each block.
            % Plotting:
            %     For each participant, mean reaction times for central and peripheral targets are plotted across blocks.
            %     This should reveal any learning effect by showing changes in RT over time for each category.
%             distance_threshold = 14; % degrees
%             viewing_distance = 30; % cm
%             screen_width_cm = 28.2; %13.5 zoll
%             screen_height_cm = 18.8; % cm
%             %plotlearningEffects(data, screenXpixels, screenYpixels, comparison_results_folder, color_map, viewing_distance, screen_width_cm);
%             plotlearningEffects(data, screen_width_cm, screen_height_cm, screenXpixels, screenYpixels, color_map, comparison_results_folder, distance_threshold)
% 


%% Interaction Between Group, Condition (Difficulty), and RT:
% H: ADHD participants might have a weaker correlation between RT_Eye and RT_ButtonPress
%    possibly indicating a delay in motor response after visual detection.
rt_button_press_conditions = cell(length(groups), length(conditions));
rt_eye_conditions = cell(length(groups), length(conditions));
accuracy_conditions = cell(length(groups), length(conditions));

for i = 1:length(data)
    group_idx = find(strcmp(groups, data(i).group));
    if isempty(group_idx)
        warning('Group label "%s" not found in the predefined groups list.', data(i).group);
        continue;
    end
    for t = 1:length(data(i).Condition)
        condition = data(i).Condition{t};
        condition_idx = find(strcmp(conditions, condition));
        if isempty(condition_idx)
            warning('Condition label "%s" not found in the predefined conditions list.', condition);
            continue;
        end
        rt_button_press_conditions{group_idx, condition_idx} = [rt_button_press_conditions{group_idx, condition_idx}; data(i).rt(t)];
        rt_eye_conditions{group_idx, condition_idx} = [rt_eye_conditions{group_idx, condition_idx}; data(i).rt_eye(t)];
        accuracy_conditions{group_idx, condition_idx} = [accuracy_conditions{group_idx, condition_idx}; data(i).accuracy(t)];
    end
end

% Visualize Reaction Times Across Conditions and Groups
% todo change condition names
% First Figure: Each group and condition in separate subplots
plotRTScatterGroupCondition(groups, conditions, rt_eye_conditions, rt_button_press_conditions, color_map, safe, comparison_results_folder);

plotRTscatterConditionAcrossGroups(groups, conditions, rt_eye_conditions, rt_button_press_conditions, color_map, safe, comparison_results_folder)



%% --- EYE-TRACKING DATA ---
dist_threshold = 50; % max distance between consecutive points (in pixels) for a fixations cluster
plot_these = [2, 119]; % just plot two for checking
comp_results_fix = strcat(comparison_results_folder, subfolder_fixation, '\');
try
    mkdir(comp_results_fix);
catch % folders already exists
end
% -01 load/detect fixations -
all_fixations = analyseFixationAllWrapper(data, plot_these, dist_threshold, screenXpixels, screenYpixels, safe, folders, subfolder_fixation,comp_results_fix,  fixation_duration);
   
% -02 load/reverse engineer "saccades" from fixation clusters -
all_saccades = analyseSaccadesAllWrapper(data, all_fixations, plot_these, screenXpixels, screenYpixels, safe, folders, subfolder_fixation, comp_results_fix, fixation_duration);

% -03 get num saccades, distance, angle, diff, cosine from trial start to target center/saccades,... -
tolerance = 200; % num pixels away to count as target found
trial_metrics = calcGazeDistanceDirection(data, all_fixations, all_saccades, tolerance, screenXpixels, screenYpixels);
% todo start at target high for some- mistake maybe??








                %% --- hamits ideas TODO ---
                % Step 1: Epoch Data Around Saccades
                % For each trial, we will create time-locked segments around each saccade:
                % Baseline (before saccade): Take pupil data for 250ms before the start of the saccade.
                % Saccade Start: The moment when the saccade begins.
                % Post-saccade: Record pupil data for 2000ms after the saccade.

                % get data
                min_length = 32;  % Datapoints; 0.016ms sampling rate -> 0.5/0.016=31.3
                num_before = 20;  % Number of datapoints before target found
                num_after = 31;   % Number of datapoints after target found
                sr = 1/60;        % Sampling rate of eyetracker
                % [result_table, tnf, tts] = findTargetAndExtractDataSaccadeStart(trial_metrics, data, comparison_results_folder, num_before, num_after, min_length);

                % baseline_window = 0.250;
                % post_saccade_window = 0.500;
                % gaze_shift_threshold = 200;
                % [epoch_data, confirmatory_epoch_data] = analyzePupilAndSaccades(data, participant_saccades, screenXpixels, screenYpixels, baseline_window, post_saccade_window, gaze_shift_threshold);

%% --- 01_fixation: Fixation Duration and Distribution ---
% get fixation durations, counts per trial ,  mean per condition
fixation_stats = getFixationDuration(all_fixations, data,conditions, comp_results_fix);

% Fixation Duration: Compare the average duration of fixations between ADHD
% and non-ADHD participants. ADHD individuals often have shorter fixation durations,
% indicating less stable attention on specific areas of interest.
plotGroupFixationDurations(fixation_stats, ids, group_labels, condition_labels, conditions, color_map_individual, color_map, comp_results_fix, fixation_duration, safe);

% Fixation Number: Analyze the number (distribution?) of fixations.
% ADHD participants might have a more scattered fixation pattern, with less
% focus on task-relevant areas compared to non-ADHD participants
%plotViolinNumFixations(fixation_stats, group_labels, conditions, color_map, comparison_results_folder, safe)
plotBarNumFixations(fixation_stats, group_labels, ids, conditions,condition_labels, color_map,color_map_individual, comp_results_fix, safe)
% -> GEGENTEIL: less fixation in condition a for adhd
% target found faster- more salient- more bottom up?

% Arrive Abandon Return Trials
% 1. If a target fixation is followed by a non-target fixation -> AARTs
% 2. If the target fixation is the last fixation -> no abandon
% 3. combo of the two
participant_data = plotAbandonReturnTrials(data_struct, all_fixations, conditions, condition_labels, group_labels, color_map, comp_results_fix, safe);
% -> no real difference in those 3 categories
% -> slightly less target fixation abandonments in condition a for adhd
% this matches with less fixations from above

% RThand - RTeye for AARTs
plotAARtrialsRThandeye(participant_data, data_struct, ids, conditions, condition_labels, group_labels, color_map, color_map_individual, comp_results_fix, safe);

                    % todo vlads thing with eccentricity
                    %     h = 25; % Monitor height in cm
                    %     d = 30; % Distance between monitor and participant in cm
                    %     r = 768; % Vertical resolution of the monitor 34.3 cm
                    %     size_in_px = 100; % The stimulus size in pixels
                    %     % Calculate the number of degrees that correspond to a single pixel. This will
                    %     % generally be a very small value, something like 0.03.
                    %     deg_per_px = degrees(atan2(.5*h, d)) / (.5*r);

%% --- 01_saccade: Saccade Metrics ---
% Saccade Amplitude and Velocity: ADHD participants may exhibit larger saccade
% amplitudes and higher velocities, reflecting more erratic or less controlled eye movements.
saccadeStats = analyzeSaccadeAmplitudeAndVelocity(all_saccades, data, conditions);
plotSaccadeDifferences(saccadeStats,data, group_labels, ids, conditions,condition_labels, color_map, color_map_individual, comp_results_fix, safe);
% -> looks basically the same, just as zhaoping said it would be
% TODO check if consistent with raw path data: distance per second 


            % Saccade Frequency: Higher saccade frequency might indicate more frequent
            % shifts in attention, which is often characteristic of ADHD.
            % -> is basically plotBarNumFixations since between each 2 fixations i assume a saccade
            %plotSaccadeCount(saccadeStats, group_labels, conditions, color_map, ids, color_map_individual, comp_results_fix, safe); 
            % -> less saccades in a/b in adhd: matches the ARTs

     
% Time to First Fixation on Target: Assess the time taken for the first
% fixation on the target. Delays in this metric could suggest differences
% in how quickly participants orient their attention to relevant stimuli.
% in easy conditions bottom up attention in harder condition maybe no
% as big influence of BU since its less salient?
plotAverageStartFirstSaccadeGroup(data, ids, conditions, condition_labels, group_labels, color_map, color_map_individual, comp_results_fix, safe);
% -> adhd slighly faster everywhere but b


%% --- Search Efficiency/Strategy/.... ---
% Path Efficiency: Measure the efficiency of the visual search path.
% Non-ADHD participants might display more direct paths towards the target,
% whereas ADHD participants could show more erratic paths with unnecessary movements.
plotCosineSimilarityGroup(data, group_labels, ids, all_saccades, all_fixations, ...
                          conditions, condition_labels, color_map, color_map_individual, comp_results_fix);
% -> distance of gaze per trial in saccade amplitude/ distance per trial plot 




%% pupil diam lets go
num_before = 20;        % Number of datapoints before target found
num_after = 31;
sr = 1/60;
time_vector = (-num_before:num_after - 1) * sr; % Time in seconds
x_label_text = 'Time from t0 (s)';
avg_result_table = [];

which_diam_before = 'tss_diam_before';
which_diam_after = 'tss_diam_after';
s_start = 'ts_start';
    
%     which_diam_before = 'ntss_diam_before';
%     which_diam_after = 'ntss_diam_after';
%     s_start = 'nts_start';
%     
    
for p=1:17
    trial_metrics = data_struct(p).trial_metrics;     
    %cut_data = data_struct(p).tss;   % only needed for conditions
    id = data_struct(p).id;
    %  --- 4.1 plot saccades TOWARDS target aligned by START of saccade [raw, 0-aligned, cleaned] ---
%     diam_t0_tss = getPupilDiamsBasedOnt0(cut_data, current_data, 'tss_diam_before','tss_diam_after', 'ts_start', comparison_results_folder, comp_results_fix, ...
%                         safe, conditions, condition_labels, time_vector, x_label_text, color_map, ...
%                         'aligned_diams', strcat('ID ', num2str(id), ': Pupil Response Around Target Saccade Starts [zero-aligned]'), ...
%                         't0_s', strcat('diam_t0_tss_zero_aligned_subj',num2str(id)), id);

    result_table = [];result_table_first= [];
    num_trials = height(trial_metrics.(which_diam_before));
    for ts = 1:num_trials
        current_condition = cut_data.Condition(ts);
        num_saccades = height(trial_metrics.(which_diam_before){ts});
        for s = 1:num_saccades
            
            %diam_combined =  
            saccade_idx = trial_metrics.(s_start){ts}(s,:);
            result_table       = [result_table; {ts, current_condition, diam_combined, saccade_idx}];
            
        end
    end
    diam_t0 = cell2table(result_table, 'VariableNames', {'TrialNumber', 'Condition', 'DataPointsBefore', 'DataPointsAfter', 'DiamCombined', 'SaccadeTrialIndex'});

    % condition avg
    for i = 1:4
        PLOTDATA = diam_t0; 
        current_condition = conditions{i}; 
        condition_trials = PLOTDATA(strcmp(PLOTDATA.Condition, current_condition), :);
        
        % Convert cell array to matrix, but first remove NaN rows
        valid_rows = ~cellfun(@(x) all(isnan(x)), condition_trials.DataPointsBefore); % Find valid rows
        before_matrix = cell2mat(condition_trials.DataPointsBefore(valid_rows)); % Convert only valid rows
        diam_before = nanmedian(before_matrix, 1); % Compute column-wise median ignoring NaNs

        valid_rows = ~cellfun(@(x) all(isnan(x)), condition_trials.DataPointsAfter); 
        after_matrix = cell2mat(condition_trials.DataPointsAfter(valid_rows));
        diam_after = nanmedian(after_matrix, 1);       
        
        diam_combined = [diam_before, diam_after];
             %saccade_idx = trial_metrics.(s_start){ts}(s,:);
        avg_result_table = [avg_result_table; {id, current_condition, diam_before, diam_after, diam_combined, 1}]; % 1 bc avg
            
    end
end
diam_t0 = cell2table(avg_result_table, 'VariableNames', {'TrialNumber', 'Condition', 'DataPointsBefore', 'DataPointsAfter', 'DiamCombined', 'SaccadeTrialIndex'});
for ts = 1:size(diam_t0, 1)
        if ~isnan(diam_t0.SaccadeTrialIndex(ts))
            try current_diam = diam_t0.DiamCombined{ts,:};
            catch, current_diam = diam_t0.DiamCombined(ts,:);
            end
            try diam_t0.aligned_diams(ts,:) = current_diam - current_diam(1);  % Aligning each trial to start at zero
            catch, diam_t0.aligned_diams(ts,:) = NaN(1,51);
            end
        else
            diam_t0.aligned_diams(ts,:) = NaN(1,51);
        end
end

adhd_ids= find(strcmp(group_labels, 'ADHD'));
trial_numbers = diam_t0.TrialNumber;
is_ADHD = ismember(trial_numbers, adhd_ids); % True if TrialNumber is in adhd_ids
is_NonADHD = ~is_ADHD; % The rest are Non-ADHD

diam_t0_ADHD = diam_t0(is_ADHD, :);     % ADHD table
diam_t0_nonADHD = diam_t0(is_NonADHD, :); % Non-ADHD table
 
%         plotDiamsPerConditionAndAverage(diam_t0_ADHD, 'aligned_diams', conditions, condition_labels, ...
%                                         time_vector, x_label_text, color_map, 'ADHD: Avg Pupil Diam After Start of Saccade to Target', ...
%                                         't0_s', comparison_results_folder, comparison_results_folder, 't0_s', id);
%
diese = 1;
if diese == 1
    PLOTDATA = diam_t0_ADHD;
    data_type =  'aligned_diams';
    xlinelabel = 't0_s';
    bigtitle  = 'ADHD: Avg Pupil Diam After Start of Saccade to Target';
    t0 = 't0_s';
    safe_name = strcat('mean_ADHD_ts_t0_diam',s_start,'.png');
else
    PLOTDATA = diam_t0_nonADHD;
    data_type =  'aligned_diams';
    xlinelabel = 't0_s';
    bigtitle  = 'nonADHD: Avg Pupil Diam After Start of Saccade to Target';
    t0 = 't0_s';
    safe_name = strcat('mean_nonADHD_ts_t0_diam',s_start,'.png');
end

num_rows = 3; num_col = 2; %figure; 
all_y_values = [];t = tiledlayout(num_rows, num_col, 'TileSpacing', 'Compact', 'Padding', 'Compact');
% Add shared title and axis labels
title(t,bigtitle); xlabel(t,x_label_text); ylabel(t,'Change in Pupil Diameter (mm)')
for i = 1:4
    ax = nexttile;     hold on;
    condition = conditions{i}; 
    title(strcat('Averages of Participants: ', condition_labels{i}));
    condition_trials = PLOTDATA(strcmp(PLOTDATA.Condition, condition), :);
    for j = 1:height(condition_trials)
        try
            y_data = condition_trials.(data_type)(j,:);
            plot(time_vector, y_data, 'Color', color_map(condition));
            all_y_values = [all_y_values; [max(y_data), min(y_data)]];  % Collect y-values for unified axis
        catch 
        end
    end
    xline(0, '--', xlinelabel,  'Color', [0.3 0.3 0.3], 'LabelOrientation', 'horizontal', 'HandleVisibility', 'off');
    hold off;
end
y_limits = [min(min(all_y_values)), max(max(all_y_values))];
for i = 1:4% Apply unified y-limits and x-limits across individual trial subplots
    ax = nexttile(i);
    ylim(ax, y_limits);
    xlim(ax, [min(time_vector), max(time_vector)]);
end
% Plot the condition average in a larger tile (spanning width of layout)
ax_avg = nexttile([1, 2]);  % Span across 2 columns
hold on; title('Condition Average');
for i = 1:4
    condition = conditions{i}; 
    condition_trials = PLOTDATA(strcmp(PLOTDATA.Condition, condition), :);
    plotConditionDataWithShading(condition_trials.(data_type), time_vector, color_map(condition), condition);
end
xline(0, '--', xlinelabel, 'LabelHorizontalAlignment', 'right', 'LabelOrientation', 'horizontal', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
lgd = legend(condition_labels); 
set(lgd, 'Orientation', 'horizontal', 'Location', 'southoutside', 'Box', 'off');
xlim(ax_avg, [min(time_vector), max(time_vector)]);
hold off;                           
 saveas(gcf, fullfile(comparison_results_folder , safe_name));
                                     
        
%         plotDiamsPerConditionAndAverage(diam_t0_NonADHD, 'aligned_diams', conditions, condition_labels, ...
%                                         time_vector, x_label_text, color_map, 'nonADHD: Avg Pupil Diam After Start of Saccade to Target', ...
%                                         't0_s', comparison_results_folder, comparison_results_folder, 't0_s', id);
%         saveas(gcf, fullfile(comparison_results_folder, strcat('nonADHD_nts_t0_diam',s_start,'.png')));
%                                     
%%  --- 4.5 plot FIRST START target saccade of trial ---
result_table = [];  
for trial=1:num_trials
    current_condition = cut_data.Condition(trial);
    current_saccades = diam_t0_tss(ismember(diam_t0_tss.TrialNumber, trial), :);
    if ~isempty(current_saccades)
        result_table = [result_table; {trial, current_condition, current_saccades.SaccadeTrialIndex(1),[current_saccades.aligned_diams(1,:)]}];
    else % if no target saccades in trial skip it
        %result_table = [result_table; {trial, NaN, NaN, [NaN(1,51)], [NaN(1,51)]}];
    end
end
diam_t0_tss_first = cell2table(result_table, 'VariableNames', {'TrialNumber', 'Condition', 'SaccadeTrialIndex', 'DiamCombined'});








for ts = 1:size(diam_t0_tss_first, 1)
    if ~isnan(diam_t0_tss_first.SaccadeTrialIndex(ts))
        try current_diam = diam_t0_tss_first.DiamCombined{ts,:};
        catch, current_diam = diam_t0_tss_first.DiamCombined(ts,:);
        end
        try diam_t0_tss_first.aligned_diams(ts,:) = current_diam - current_diam(1);  % Aligning each trial to start at zero
        catch, diam_t0_tss_first.aligned_diams(ts,:) = NaN(1,51);
        end
    else
        diam_t0_tss_first.aligned_diams(ts,:) = NaN(1,51);
    end
end

diam_t0_tss_first.cleanedData = filloutliers(diam_t0_tss_first.aligned_diams,'nearest','percentiles',[10 90]); %remove outliers
save(fullfile(analysis_folder, strcat('\diam_t0_tss_first.mat')), 'diam_t0_tss_first');
if do_plots == 1
    plotCompareTwoStepsDiam(diam_t0_tss_first, 'aligned_diams', 'ZeroAligned', ...
                            diam_t0_tss_first, 'cleanedData',   '80Percentile', ...
                            strcat('ID ', num2str(id), ': Pupil Response Around First Target Saccade Start'), ...
                            color_map, conditions, condition_labels, time_vector, x_label_text, 't0_s',analysis_folder, compare_folder, 'tss_first',id);
    plotDiamsPerConditionAndAverage(diam_t0_tss_first, 'aligned_diams', conditions,  condition_labels, time_vector, x_label_text, color_map, ...
            strcat('ID ', num2str(id), ': Pupil Response Around First Target Saccade Start [ZeroAligned]'), 't0_s',...
            analysis_folder, compare_folder, 'tss_first',id);
else
end                                    
        

%% averages for both groups
safe_name_group = strcat('ts_t0_diam',s_start,'_groups_mean.png');

num_rows = 1; 
num_col = 2; 
figure; 
t = tiledlayout(num_rows, num_col, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Common X-axis label
x_label_text = 'Time (s)';
y_label_text = 'Change in Pupil Diameter (mm)';

% ADHD Plot
PLOTDATA = diam_t0_ADHD;
bigtitle_ADHD  = 'ADHD: Avg Pupil Diam After Start of Saccade to Target';
nexttile; hold on;
title(bigtitle_ADHD);
xlabel(x_label_text); ylabel(y_label_text);

for i = 1:4
    condition = conditions{i}; 
    condition_trials = PLOTDATA(strcmp(PLOTDATA.Condition, condition), :);
    plotConditionDataWithShading(condition_trials.(data_type), time_vector, color_map(condition), condition);
end

xline(0, '--', xlinelabel, 'LabelHorizontalAlignment', 'right', 'LabelOrientation', 'horizontal', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
hold off;

% Non-ADHD Plot
PLOTDATA = diam_t0_NonADHD;
bigtitle_NonADHD  = 'Non-ADHD: Avg Pupil Diam After Start of Saccade to Target';
nexttile; hold on;
title(bigtitle_NonADHD);
xlabel(x_label_text); ylabel(y_label_text);

for i = 1:4
    condition = conditions{i}; 
    condition_trials = PLOTDATA(strcmp(PLOTDATA.Condition, condition), :);
    plotConditionDataWithShading(condition_trials.(data_type), time_vector, color_map(condition), condition);
end

xline(0, '--', xlinelabel, 'LabelHorizontalAlignment', 'right', 'LabelOrientation', 'horizontal', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
lgd = legend(condition_labels, 'Orientation', 'horizontal', 'Location', 'southoutside', 'Box', 'off');
hold off;

title(t, 'Comparison of Pupil Diameter Changes: ADHD vs nonADHD');
set(gcf, 'Position', [50, 50, 1400, 700]); % Resize the figure window (x, y, width, height)
saveas(gcf, fullfile(comparison_results_folder , safe_name_group));

                    
          





















    %% Fill outliers
    diam_t0_tss.cleanedData = filloutliers(diam_t0_tss.aligned_diams,'nearest', 'percentiles',[10 90]); 
    if do_plots == 1
        plotCompareTwoStepsDiam(diam_t0_tss, 'aligned_diams', 'ZeroAligned',  ...
                                diam_t0_tss, 'cleanedData',   '80Percentile', ...
                                strcat('ID ', num2str(id), ': Pupil Response Around Target Saccade Starts'), ...
                                color_map, conditions, condition_labels, time_vector, x_label_text, 't0_s',analysis_folder, compare_folder, 'tss',id);
        plotHistNaNsBufferedDiamSacc(diam_t0_tss, analysis_folder, 'Target Saccade NaN Counts before/after t0start', 'diam_t0_tss_nan_counts_histogram');
else
    end

%save(fullfile(analysis_folder, strcat('\diam_t0_tss.mat')), 'diam_t0_tss');

        % Reverse engineer Saccades from Fixations and take as t0:
        % 1. [Start of Saccade that Reached Target]
        % 2. [End of Saccade that Reached Target]
        % 3. [Start of Saccade that didnt Reach Target]
        % 4. [End of Saccade that didnt Reach Target]
        % 5. [FIRST Start Target Saccade of trial]
        % 6. [FIRST End of Target Saccade of trial]
        % 7. [LAST Start of Target Saccades per trial] (uneccessary bc fixation as marker for fixation ensures the saccade was goal oriented)
        analyseChangePupilDilation(id, cut_data, trial_metrics, conditions, condition_labels, analysis_folder, compare_folder, do_plots, time_vector, x_label_text, color_map);
        close all; % close plots
        
        
        
        

    %% --- Normalized and Max Pupil Diameter  ---
    % get pupil diam of all participants by group

    % --- plot nRT versus maxPupilDiam after target reaches stim ---
    % Set threshold for outlier detection
    % 2 plots: per group and per participant
    outlier_threshold = 4; % remove outliers from rt/rtv, already done for pupil dilation in analyseResults
    [avg_adhd, avg_nonadhd] = extractPupilData(data, outlier_threshold, 'ADHD');
    plotPupilRTRelationships(avg_adhd, avg_nonadhd, conditions, color_map, comparison_results_folder);

    % --- pupil diameter over time per group with stdabweichung ---
    % get data into format needed for plot
    outlier_threshold = 4;
    [avg_adhd, avg_nonadhd] = extractPupilData(data, outlier_threshold, 'ADHD');
    before = size(data(1).pupilDiamNorm.a.beforeMean,2);
    after = size(data(1).pupilDiamNorm.a.afterMean,2);
    [avg_adhd_struct, avg_nonadhd_struct] = changeFormatPupilRT(avg_adhd, avg_nonadhd, before, after);
    plotAvgPupilDiamsBeforeAfterStimFoundByGroupByCondition('Mean', 'Mean Change in Pupil Diameter (z-score)', ...
        avg_adhd_struct, avg_nonadhd_struct, conditions, color_map, before, after, outlier_threshold, comparison_results_folder);
    avg_both = [avg_adhd; avg_nonadhd];

                %% MS participants
                outlier_threshold = 3;
                [avg_ms, avg_nonms] = extractPupilData(data, outlier_threshold, 'MS');
                avg_both = [avg_ms; avg_nonms];
                before = size(data(1).pupilDiamNorm.a.beforeMean,2);
                after = size(data(1).pupilDiamNorm.a.afterMean,2);
                [avg_ms_struct, avg_nonms_struct] = changeFormatPupilRT(avg_ms, avg_nonms, before, after);
                plotAvgPupilDiamsBeforeAfterStimFoundByGroupByCondition('Mean', 'Mean Change in Pupil Diameter (z-score)', avg_ms_struct, avg_nonms_struct, conditions, color_map, before, after, comparison_results_folder);
%%
% --- max pupil diam after stim found per group ---
plotPupilMaxStimFound(avg_adhd, avg_nonadhd, conditions, color_map, comparison_results_folder);
% --- max min and mean pupil diam during stimulus presentation ---
plotPupilMagnitudesTwoGroups(avg_adhd, avg_nonadhd, conditions, color_map, comparison_results_folder)
% --- anova ---
% TODO FIX
%repeatedMeasuresANOVA(avg_adhd, avg_nonadhd, conditions, comparison_results_folder)








% Search Latency: Analyze the latency before initiating a search, and how quickly
% participants give up on a search (indicative of persistence and focus).
%searchLatencyStats = analyzeSearchLatency(all_fixations, data, conditions, comparison_results_folder);
%nanmean(searchLatencyStats(1).conditions(1).searchLatencies)
% TODO something is wrong here

%% --- Attention Shifts and Disengagement ---
% Task Engagement Over Time: Evaluate how task engagement changes over time.
% ADHD participants might show a decline in performance or more frequent disengagement from the task as it progresses.

% Revisits: Track how often participants revisit previously searched areas.
% Frequent revisits might indicate memory deficits or difficulties in sustaining attention.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QUESTIONAIRE RESULTS
% Define the response mapping
response_mapping = containers.Map({'Never', 'Rarely', 'Sometimes', 'Often', 'Very Often'}, [1, 2, 3, 4, 5]);
% Create a generic set of column names including the ID column
num_questions = 18;
column_names = [{'id'}, arrayfun(@(x) sprintf('Question%d', x), 1:num_questions, 'UniformOutput', false)];
question_data = loadQuesionnaires(folders, group_labels);% load questionnaire results
adhd_data = question_data.ADHD;
non_adhd_data = question_data.nonADHD;
quest_table = [adhd_data; non_adhd_data];
% Convert responses to numeric for easier analysis and plotting
adhd_numeric = adhd_data{:, 2:end};  % Skip the ID column
non_adhd_numeric = non_adhd_data{:, 2:end};  % Skip the ID column
fprintf('Loaded %d ADHD entries and %d nonADHD entries.\n', height(adhd_data), height(non_adhd_data));

% get questionnaire results/mean answers
quest_scores = questionnaireScale(quest_table, comparison_results_folder);

% Pie Chart for overall response proportions for (non)ADHD group
plotQPieChart(adhd_numeric, 'ADHD', color_map, 'all', safe, comparison_results_folder);
plotQPieChart(non_adhd_numeric, 'nonADHD', color_map, 'all', safe, comparison_results_folder);

% Pie Chart for first 6 responses: proportions for (non)ADHD group
plotQPieChart(adhd_numeric(:,1:6), 'ADHD', color_map, '6', safe, comparison_results_folder);
plotQPieChart(non_adhd_numeric(:,1:6), 'nonADHD', color_map, '6', safe, comparison_results_folder);

% Answers per question with a unique color for each participant
plotQParticipantAnswers(quest_table, color_map, safe, comparison_results_folder);

% Means of answers per question per group
plotQMeanAndStd(adhd_numeric, non_adhd_numeric, color_map, safe, comparison_results_folder)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INFLUENCE OF (IN)ATTENTION
%% combine questionnaire results with RT and accuracy
% Convert your data structure to a table for easier manipulation
data_table = struct2table(data);

% Merge questionnaire data with RT/accuracy/saccade data based on ID
merged_data = innerjoin(data_table, quest_table, 'Keys', 'id');
merged_data = innerjoin(merged_data, quest_scores, 'Keys', 'id');
%merged_data = innerjoin(merged_data, struct2table(rt_variability), 'Keys', 'id');
merged_data = innerjoin(merged_data, struct2table(trial_metrics), 'Keys', 'id');


%% define attentional levels based on ASRS guidelines
% define thresholds for categorizing attention levels (based on guidelines)
thresholds = [0, 4, 6];  % Thresholds between categories
%attention_levels = discretize(merged_data.PartASymptomsCount, thresholds, 'categorical', {'Low', 'Moderate', 'High'});
attention_levels = discretize(merged_data.PartASymptomsCount, thresholds, 'IncludedEdge', 'right');
% Add the attention level to your data
merged_data.InattentionLevel = attention_levels;

% Separate data into two groups
high_indices = ismember(merged_data.InattentionLevel, [2,3]);
low_indices = merged_data.InattentionLevel == 1;

adhd_indices = strcmp(merged_data.group, 'ADHD');
nonadhd_indices = strcmp(merged_data.group, 'nonADHD');


                %% correlation analyis (unpassend zum vergleich der gruppen, wenn dann innerhalb der gruppen)
                % ttest between AHD/nonADHD groups AND
                % ttestResultsADHD = cleanADHDnonADHDattRT(merged_data, conditions);
                % ttestResults = cleanAttRTcorrTtest(merged_data, conditions);


%% --- lieber ttest ---
% --- ADHD-nonADHD ---
% % RT
% [~, p_rt_a] = ttest2(merged_data.nRTa(adhd_indices), merged_data.nRTa(nonadhd_indices));
% [~, p_rt_b] = ttest2(merged_data.nRTb(adhd_indices), merged_data.nRTb(nonadhd_indices));
% [~, p_rt_asimple] = ttest2(merged_data.nRTasimple(adhd_indices), merged_data.nRTasimple(nonadhd_indices));
% [~, p_rt_bsimple] = ttest2(merged_data.nRTbsimple(adhd_indices), merged_data.nRTbsimple(nonadhd_indices));
% fprintf('T-test p-values for RT comparisons between ADHD and non-ADHD:\n');
% fprintf('Condition A: p = %.3f\n', p_rt_a);
% fprintf('Condition B: p = %.3f\n', p_rt_b);
% fprintf('Condition A Simple: p = %.3f\n', p_rt_asimple);
% fprintf('Condition B Simple: p = %.3f\n', p_rt_bsimple);
% RTV
[~, p_rtv_a] = ttest2(cell2mat({merged_data.a(adhd_indices,1).RT_ButtonPress_var}), cell2mat({merged_data.a(nonadhd_indices,1).RT_ButtonPress_var}));
[~, p_rtv_b] = ttest2(cell2mat({merged_data.b(adhd_indices,1).RT_ButtonPress_var}), cell2mat({merged_data.b(nonadhd_indices,1).RT_ButtonPress_var}));
[~, p_rtv_asimple] = ttest2(cell2mat({merged_data.a_simple(adhd_indices,1).RT_ButtonPress_var}), cell2mat({merged_data.a_simple(nonadhd_indices,1).RT_ButtonPress_var}));
[~, p_rtv_bsimple] = ttest2(cell2mat({merged_data.b_simple(adhd_indices,1).RT_ButtonPress_var}), cell2mat({merged_data.b_simple(nonadhd_indices,1).RT_ButtonPress_var}));
fprintf('T-test p-values for RTV comparisons between ADHD and non-ADHD:\n');
fprintf('Condition A: p = %.3f\n', p_rtv_a);
fprintf('Condition B: p = %.3f\n', p_rtv_b);
fprintf('Condition A Simple: p = %.3f\n', p_rtv_asimple);
fprintf('Condition B Simple: p = %.3f\n', p_rtv_bsimple);
% inattention levels
[~, p_a_inatt_count] = ttest2(merged_data.PartASymptomsCount(adhd_indices), merged_data.PartASymptomsCount(nonadhd_indices));
fprintf('T-test p-values comparisons ADHD/ non-ADHD: Inattention - PartASymptomsCount: p = %.3f\n', p_a_inatt_count);
% confusion
[~, p_conf_diff_button] = ttest2(merged_data.confusion_diff_press(adhd_indices), merged_data.confusion_diff_press(nonadhd_indices));
[~, p_conf_ratio_button] = ttest2(merged_data.confusion_ratio_press(adhd_indices), merged_data.confusion_ratio_press(nonadhd_indices));
fprintf('T-test p-values comparisons ADHD/non-ADHD: Confusion Difference RT Button Press: p = %.3f\n', p_conf_diff_button);
fprintf('T-test p-values comparisons ADHD/non-ADHD: Confusion Ratio RT Button Press: p = %.3f\n', p_conf_ratio_button);
[~, p_conf_diff_button] = ttest2(merged_data.confusion_diff_eye(adhd_indices), merged_data.confusion_diff_eye(nonadhd_indices));
[~, p_conf_ratio_button] = ttest2(merged_data.confusion_ratio_eye(adhd_indices), merged_data.confusion_ratio_eye(nonadhd_indices));
fprintf('T-test p-values comparisons ADHD/non-ADHD: Confusion Difference RT Eye: p = %.3f\n', p_conf_diff_button);
fprintf('T-test p-values comparisons ADHD/non-ADHD: Confusion Ratio RT Eye: p = %.3f\n', p_conf_ratio_button);

% accuracy

% saccade metrics
ttestSaccadeMetrics(merged_data, 'a', 'ADHD/non-ADHD', adhd_indices, nonadhd_indices);
ttestSaccadeMetrics(merged_data, 'b', 'ADHD/non-ADHD', adhd_indices, nonadhd_indices);
ttestSaccadeMetrics(merged_data, 'a_simple', 'ADHD/non-ADHD', adhd_indices, nonadhd_indices);
ttestSaccadeMetrics(merged_data, 'b_simple', 'ADHD/non-ADHD', adhd_indices, nonadhd_indices);
%ttestSaccadeMetrics(merged_data, 'a', 'high inattention/low inattention', high_indices, low_indices)


% Extract maxMean from .a field for each ADHD participant
pupilDiamNorm_adhd = arrayfun(@(x) nanmean(x.a.maxMean), merged_data.pupilDiamNorm(adhd_indices));
pupilDiamNorm_nonadhd = arrayfun(@(x) nanmean(x.a.maxMean), merged_data.pupilDiamNorm(nonadhd_indices));
[~, p_pupilDiamNorm] = ttest2(pupilDiamNorm_adhd, pupilDiamNorm_nonadhd);
fprintf('Condition A p_pupilDiamNorm: p = %.3f\n', p_pupilDiamNorm);



% --- high inattention vs low inattention levels ---
% RTV
[~, p_rtv_a_inatt] = ttest2(cell2mat({merged_data.a(high_indices,1).RT_ButtonPress_var}), cell2mat({merged_data.a(low_indices,1).RT_ButtonPress_var}));
[~, p_rtv_b_inatt] = ttest2(cell2mat({merged_data.b(high_indices,1).RT_ButtonPress_var}), cell2mat({merged_data.b(low_indices,1).RT_ButtonPress_var}));
[~, p_rtv_asimple_inatt] = ttest2(cell2mat({merged_data.a_simple(high_indices,1).RT_ButtonPress_var}), cell2mat({merged_data.a_simple(low_indices,1).RT_ButtonPress_var}));
[~, p_rtv_bsimple_inatt] = ttest2(cell2mat({merged_data.b_simple(high_indices,1).RT_ButtonPress_var}), cell2mat({merged_data.b_simple(low_indices,1).RT_ButtonPress_var}));
fprintf('T-test p-values for RTV comparisons between high and low attention groups:\n');
fprintf('Condition A: p = %.3f\n', p_rtv_a_inatt);
fprintf('Condition B: p = %.3f\n', p_rtv_b_inatt);
fprintf('Condition A Simple: p = %.3f\n', p_rtv_asimple_inatt);
fprintf('Condition B Simple: p = %.3f\n', p_rtv_bsimple_inatt);





%%
function sig_level = getSignificanceStars(p_value)
% Define significance stars based on p-value
if p_value < 0.001
    sig_level = '***';
elseif p_value < 0.01
    sig_level = '**';
elseif p_value < 0.05
    sig_level = '*';
else
    sig_level = '';
end
end

function addSigStars(x, sig_level, y_offset)
% Add significance stars to plot at given x positions
if ~isempty(sig_level)
    text(mean(x), max(ylim) + y_offset, sig_level, 'HorizontalAlignment', 'center', 'FontSize', 12);
end
end


function [y_upper,y_lower] = plotConditionDataWithShading(data, x_values, color, condition)
    % Plot median ± SEM with shading for each condition
    avg = mean(data, 1, 'omitnan');
    std_error = std(data, 0, 1, 'omitnan') / sqrt(size(data, 1));
    nan_mask = ~isnan(avg) & ~isnan(std_error);
    x_shaded = x_values(nan_mask);
    y_upper = avg(nan_mask) + std_error(nan_mask);
    y_lower = avg(nan_mask) - std_error(nan_mask);
    fill([x_shaded, fliplr(x_shaded)], [y_upper, fliplr(y_lower)], color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(x_values, avg, 'Color', color, 'LineWidth', 2, 'DisplayName', condition);
end