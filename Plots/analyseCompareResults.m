%% load all data & define variables
% save path
comparison_results_folder = 'C:\Users\flohrmann\Documents\Analysis\';
% save plots yes/no 1/0
safe = 1;
% colour scheme
color_map = containers.Map({'a', 'a_simple', 'b', 'b_simple', 'ADHD', 'nonADHD'}, {
    [0.9, 0.5, 0  ]  % orange
    [1.0, 0  , 0  ]  % red
    [0  , 0  , 1  ]  % blue
    [0  , 0.5, 0.1]  % green
    [0.5, 0.7, 0.5]  % orange ish
    [0.7, 0.2, 0.2]  % green ish
    });

folders = {
    'C:\Users\flohrmann\Documents\Results\1_20240714_140048'; % alex
    'C:\Users\flohrmann\Documents\Results\2_20240726_191755'; % mara adhd
    'C:\Users\flohrmann\Documents\Results\3_20240805_105213'; % tilo adhd
    'C:\Users\flohrmann\Documents\Results\4_20240811_131601'; % anu adhd
    'C:\Users\flohrmann\Documents\Results\5_20240813_114700'; % kieran adhd
    };

% Define the classification of the folders (ADHD or non-ADHD)
group_labels = {
    'nonADHD';  % alex
    'ADHD';     % mara adhd
    'ADHD';     % tilo adhd
    'ADHD';     % anu adhd
    'ADHD';     % kieran adhd
    };

conditions = {'a', 'a_simple', 'b', 'b_simple'};
groups = {'ADHD', 'nonADHD'};


data_struct = struct();
% Loop through each folder and load the data
for i = 1:length(folders)
    folder = folders{i};
    %folder = strcat(folders{i}, '\results');
    group = group_labels{i};
    data_file   = dir(fullfile(folder, '\results\trial_results*.mat')); % original data
    eye_rt_file = dir(strcat(folder, '\analysis\eye_rt.mat')); % calculated eye rt from analyseResults.m
    eyetracking_file = dir(fullfile(folder, '\analysis\cut_trials*.mat'));
    if ~isempty(data_file)
        load(fullfile(data_file.folder, data_file.name), 'trial_results');
        data_struct(i).id                = i;
        data_struct(i).group             = group;
        data_struct(i).Condition         = trial_results.Condition;
        data_struct(i).rt                = trial_results.rt;
        data_struct(i).accuracy          = trial_results.correct;
        data_struct(i).StimulusOnsetTime = double(trial_results.StimulusOnsetTime/ 1e6);
        data_struct(i).trialEndTime      = double(trial_results.trialEndTime/ 1e6);
        % eye reaction times
        load(fullfile(eye_rt_file.folder, eye_rt_file.name), 'eye_rt');
        data_struct(i).rt_right          = eye_rt.RightEyeRT;
        data_struct(i).rt_left           = eye_rt.LeftEyeRT;
        % eyetracking data cut into trials
        load(fullfile(eyetracking_file.folder, eyetracking_file.name), 'cutData');
        data_struct(i).eyeTrial           = cutData.eyeTrial;
        data_struct(i).stimulusTrial      = cutData.stimulusTrial;
        data_struct(i).TargetPosition     = cutData.TargetPosition;
        data_struct(i).x_centers          = cutData.x_centers;
        data_struct(i).y_centers          = cutData.y_centers;
    else
        warning(['Data file not found in folder: ', folder]);
    end
end


% get the faster of the two eye RTs; ignore 0s where target was not looked at
data_struct = processEyeRTs(data_struct);
% normalized RT per participant
data_struct_norm = normalizeRTsBySimpleConditions(data_struct, 'a_simple', 'b_simple', comparison_results_folder);
data_struct_norm_mean = normalizeMeanRTsBySimpleConditions(data_struct, comparison_results_folder);

% set data to normalized data
data = data_struct_norm_mean;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REACTION TIME AND ACCURACY
%% Compare RT_ButtonPress and RT_Eye across the four conditions ("a", "a simple", "b", "b simple").
%% Hyp: Easier conditions ("a simple", "b simple") might have faster RTs than harder conditions ("a", "b").
%% Hyp: ADHD participants might show different reaction time patterns compared
%%      to non-ADHD participants, especially in harder conditions.

%% Compare Reaction Times Across Conditions and Groups
% H: ADHD participants might show different reaction time patterns compared
% to non-ADHD participants, especially in harder conditions

% correlation analysis: RTeye, RTbuttonpress and accuracy
[R, P] = correlationRTbuttonRTeyeAccuracy(data, comparison_results_folder);

% BAR and LINE plot: mean RT and error bars (~variability) per group and condition
plotBarSEMMeanRTConditonGroup(groups, conditions, data, color_map, safe, comparison_results_folder)

% rt over time per participant
%% TODO change to subplots per condition
plotRTperParticipantOverTimeColouredByGroup(data, color_map, safe, comparison_results_folder)

% violin plot rt time across groups and conditions
plotViolinMeanRTGroupCondition(data, color_map, comparison_results_folder, safe)
plotViolinRTGroupCondition(data, color_map, comparison_results_folder, safe)

% NOT NEEDED normalized RT across groups
%plotNormalizedRTCcomparison(data, color_map, comparison_results_folder, safe);

% rt button press vs rt eyes; scatter plot
plotRTButtonPressVsEye(data, comparison_results_folder, safe)

% 3x3 subplots: rt button press vs rt eyes; means of groups per condition
plotMeanRTButtonPressVsEyecomparison(data, color_map, comparison_results_folder, safe)

%% Confusion/ RTa - RTb or RTa/RTb
plotConfusion(data, color_map, comparison_results_folder, safe)


%% Compare Accuracy Across Conditions and Groups
% H: Accuracy might be lower in harder conditions for ADHD participants
% compared to non-ADHD participants
% TODO add error bars
plotMeanAccuracyPerGroupCondition(data, color_map, safe, comparison_results_folder)

% mean accuracy comparison group and condition + mean points
%plotBarAccuracyGroupCondition(data, color_map, comparison_results_folder, safe)
% TODO: change bars to violin!

% accuracy vs rt button press per condition and group
plotAccuracyVsButtonPressRT(data, color_map, comparison_results_folder, safe);

% scatter plot accuracy vs RT all trials per person
figure; hold on;
for i = 1:length(data)
    scatter(data(i).rt, data(i).accuracy, 'DisplayName', data(i).group);
end
xlabel('Reaction Time'); ylabel('Accuracy'); title('Accuracy vs. Reaction Time');
legend('Location', 'northeastoutside'); hold off;

% TODO  outlier analysis
% all_rt = [];
% all_conditions = [];
% all_group = [];
% for i = 1:length(data)
%     n = length(data(i).rt);
%     all_rt = [all_rt; data(i).rt];
%     all_conditions = [all_conditions; data(i).Condition];
%     all_group = [all_group; repmat({data(i).group}, n, 1)];
% end
%
% % Identify outliers
% is_outlier = isoutlier(all_rt);
% outlier_data = table(all_group(is_outlier), all_conditions(is_outlier), all_rt(is_outlier), ...
%     'VariableNames', {'Group', 'Condition', 'ReactionTime'});
% disp(outlier_data);


%% Interaction Between Group, Condition (Difficulty), and RT:
% H: ADHD participants might have a weaker correlation between RT_Eye and RT_ButtonPress
%    possibly indicating a delay in motor response after visual detection.
rt_button_press_conditions = cell(length(groups), length(conditions));
rt_eye_conditions = cell(length(groups), length(conditions));
accuracy_conditions = cell(length(groups), length(conditions));

% Loop through each participant
for i = 1:length(data)
    % Find the group index
    group_idx = find(strcmp(groups, data(i).group));
    
    if isempty(group_idx)
        warning('Group label "%s" not found in the predefined groups list.', data(i).group);
        continue;
    end
    
    % Loop through each trial for the participant
    for t = 1:length(data(i).Condition)
        condition = data(i).Condition{t};
        condition_idx = find(strcmp(conditions, condition));
        
        if isempty(condition_idx)
            warning('Condition label "%s" not found in the predefined conditions list.', condition);
            continue;
        end
        
        % Store data based on condition and group
        rt_button_press_conditions{group_idx, condition_idx} = [rt_button_press_conditions{group_idx, condition_idx}; data(i).rt(t)];
        rt_eye_conditions{group_idx, condition_idx} = [rt_eye_conditions{group_idx, condition_idx}; data(i).rt_eye(t)];
        accuracy_conditions{group_idx, condition_idx} = [accuracy_conditions{group_idx, condition_idx}; data(i).accuracy(t)];
    end
end

% Visualize Reaction Times Across Conditions and Groups
plotScatterRTeyeRTbuttonGroupCondition(groups, conditions, rt_eye_conditions, rt_button_press_conditions, color_map, safe, comparison_results_folder)

%% ANOVA
% Two-way ANOVA, where one factor is the group (ADHD vs. non-ADHD)
% and the other factor is the condition (e.g., "a", "a simple", "b", "b simple")
% Assuming your data is already extracted and organized
% Group data
%conditions = {'a', 'a_simple', 'b', 'b_simple'};
all_group_labels = [];   % Vector to store group labels (ADHD or nonADHD)
all_condition_labels = []; % Vector to store condition labels (a, a simple, b, b simple)
reaction_times = []; % Vector to store reaction times

% Loop through each participant
for i = 1:length(data)
    num_trials = length(data(i).Condition);
    % Extend vectors with data from this participant
    all_group_labels = [all_group_labels; repmat({data(i).group}, num_trials, 1)];
    all_condition_labels = [all_condition_labels; data(i).Condition];
    reaction_times = [reaction_times; data(i).rt]; % or use data(i).accuracy for accuracy analysis
end

% Convert group and condition labels to categorical
all_group_labels = categorical(all_group_labels);
all_condition_labels = categorical(all_condition_labels);

% Perform two-way ANOVA
[p, tbl, stats] = anovan(reaction_times, {all_group_labels, all_condition_labels}, ...
    'model', 'interaction', ...
    'varnames', {'Group', 'Condition'});

% Display the ANOVA table
disp(tbl);

% Perform post-hoc tests if needed
c = multcompare(stats, 'Dimension', [1 2]);

%% Repeated measures ANOVA
rt_table = table;
% Loop through each participant
for i = 1:length(data)
    participant_rt = data(i).rt;
    conditions = data(i).Condition;
    % Organize reaction times by condition
    rt_a = participant_rt(strcmp(conditions, 'a'));
    rt_asimple = participant_rt(strcmp(conditions, 'a_simple'));
    rt_b = participant_rt(strcmp(conditions, 'b'));
    rt_bsimple = participant_rt(strcmp(conditions, 'b_simple'));
    rt_table = [rt_table; table({data(i).group}, mean(rt_a), mean(rt_asimple), mean(rt_b), mean(rt_bsimple), ...
        'VariableNames', {'Group', 'A', 'ASimple', 'B', 'BSimple'})];
end

rt_table.Group = categorical(rt_table.Group);

% Fit repeated measures model
rm = fitrm(rt_table, 'A-BSimple ~ Group', 'WithinDesign', table({'A'; 'ASimple'; 'B'; 'BSimple'}, 'VariableNames', {'Condition'}));

% Perform repeated measures ANOVA
% Sphericity Assumption: Repeated measures ANOVA assumes sphericity
% (the variances of the differences between all combinations of related groups are equal)
% MATLAB will automatically test this and apply corrections if needed.
ranova_results = ranova(rm, 'WithinModel', 'Condition');
disp(ranova_results);

% Post-hoc comparisons if needed
multcompare(rm, 'Condition');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QUESTIONAIRE RESULTS
%%
% Define the response mapping
response_mapping = containers.Map({'Never', 'Rarely', 'Sometimes', 'Often', 'Very Often'}, [1, 2, 3, 4, 5]);
% Create a generic set of column names including the ID column
num_questions = 18;
column_names = [{'id'}, arrayfun(@(x) sprintf('Question%d', x), 1:num_questions, 'UniformOutput', false)];

% load questionnaire results
question_data = loadQuesionnaires(folders, group_labels);
adhd_data = question_data.ADHD;
non_adhd_data = question_data.nonADHD;
quest_table = [adhd_data; non_adhd_data];
% Convert responses to numeric for easier analysis and plotting
adhd_numeric = adhd_data{:, 2:end};  % Skip the ID column
non_adhd_numeric = non_adhd_data{:, 2:end};  % Skip the ID column
fprintf('Loaded %d ADHD entries and %d nonADHD entries.\n', height(adhd_data), height(non_adhd_data));

% get questionnaire results/mean answers
quest_scores = questionnaireScale(quest_table, comparison_results_folder);


%% plots
% Pie Chart for overall response proportions for (non)ADHD group
plotQPieChart(adhd_numeric, 'ADHD', color_map, safe, comparison_results_folder);
plotQPieChart(non_adhd_numeric, 'nonADHD', color_map, safe, comparison_results_folder);

% Answers per question with a unique color for each participant
plotQParticipantAnswers(quest_table, color_map, safe, comparison_results_folder);

% Means of answers per question per group
plotQMeanAndStd(adhd_numeric, non_adhd_numeric, color_map, safe, comparison_results_folder)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INFLUENCE OF (IN)ATTENTION
%% combine questionnaire results with RT and accuracy
% Convert your data structure to a table for easier manipulation
data_table = struct2table(data);

% Merge questionnaire data with RT/accuracy data based on ID
merged_data = innerjoin(data_table, quest_table, 'Keys', 'id');
merged_data = innerjoin(merged_data, quest_scores, 'Keys', 'id');

%% define if adhd or not (make function based on questionnaire!!)
% Define thresholds for categorizing attention levels (these are arbitrary and can be adjusted)
thresholds = [0.5, 3.5, 4.5, 5.5];  % Thresholds between categories
attention_levels = discretize(merged_data.PartASymptomsCount, thresholds, 'categorical', {'Low', 'Moderate', 'High'});

% Add the attention level to your data
merged_data.InattentionLevel = attention_levels;

% Initialize empty arrays to hold mean RTs and accuracies
mean_rt_adhd = [];
mean_accuracy_adhd = [];
mean_rt_nonadhd = [];
mean_accuracy_nonadhd = [];
mean_inatt_adhd = [];
mean_inatt_nonadhd = [];
mean_con_adhd = [];
mean_con_nonadhd = [];

% Loop through each participant and calculate mean RT and accuracy
for i = 1:size(merged_data,1)
    % Check if the participant belongs to the ADHD group
    %if strcmp(merged_data.group{i}, 'ADHD')
    % the mean RT
    rt_values = [merged_data.nRTa(i), merged_data.nRTb(i), merged_data.nRTasimple(i), merged_data.nRTbsimple(i)];
    mean_rt_adhd = [mean_rt_adhd; mean(rt_values, 'omitnan')];
    % mean accuracy
    mean_accuracy_adhd = [mean_accuracy_adhd; mean(merged_data.accuracy{i}, 'omitnan')];
    % confusion
    mean_con_adhd  = [mean_con_adhd; merged_data.nRTa(i)/merged_data.nRTb(i)];
    % inattention level
    mean_inatt_adhd = [mean_inatt_adhd; merged_data.PartAMeanSymptoms(i)];
    
    % Check if the participant belongs to the non-ADHD group
    %elseif strcmp(merged_data.group{i}, 'nonADHD')
    % mean RT
    %         rt_values = [merged_data.nRTa(i), merged_data.nRTb(i), merged_data.nRTasimple(i), merged_data.nRTbsimple(i)];
    %         mean_rt_nonadhd = [mean_rt_nonadhd; mean(rt_values, 'omitnan')];
    %         % confusion
    %         mean_con_nonadhd  = [mean_con_nonadhd; merged_data.nRTa(i)/merged_data.nRTb(i)];
    %         % mean accuracy
    %         mean_accuracy_nonadhd = [mean_accuracy_nonadhd; mean(merged_data.accuracy{i}, 'omitnan')];
    %         % inattention level
    %         mean_inatt_nonadhd = [mean_inatt_nonadhd; merged_data.PartAMeanSymptoms(i)];
    %     end
end

%% Correlation Analysis
% Correlation Analysis for BOTH group
[r_adhd_rt, p_adhd_rt] = corr(mean_inatt_adhd, mean_rt_adhd);
[r_adhd_con, p_adhd_con] = corr(mean_inatt_adhd, mean_con_adhd);
[r_adhd_acc, p_adhd_acc] = corr(mean_inatt_adhd, mean_accuracy_adhd);

% Display results for BOTH group
fprintf('Correlation between inattention levels and RT: r = %.2f, p = %.3f\n', r_adhd_rt, p_adhd_rt);
fprintf('Correlation between inattention levels and Accuracy: r = %.2f, p = %.3f\n', r_adhd_acc, p_adhd_acc);
fprintf('Correlation between inattention levels and Confusion: r = %.2f, p = %.3f\n', r_adhd_con, p_adhd_con);


% % Correlation Analysis for non-ADHD group
% [r_nonadhd_rt, p_nonadhd_rt] = corr(mean_inatt_nonadhd, mean_rt_nonadhd);
% [r_nonadhd_con, p_nonadhd_con] = corr(mean_inatt_nonadhd, mean_con_nonadhd);
% [r_nonadhd_acc, p_nonadhd_acc] = corr(mean_inatt_nonadhd, mean_accuracy_nonadhd);
%
% % Display results for non-ADHD group
% fprintf('Correlation between non-ADHD questionnaire scores and RT: r = %.2f, p = %.3f\n', r_nonadhd_rt, p_nonadhd_rt);
% fprintf('Correlation between non-ADHD questionnaire scores and Accuracy: r = %.2f, p = %.3f\n', r_nonadhd_acc, p_nonadhd_acc);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Saccades


%% interpolated saccade data

% Example timestamps vector (your actual data will be different)
timestamps = double(data(1).eyeTrial(1).systemTimeStamp) / 1000;  % Convert to milliseconds

% Check the range of your timestamps
start_time = timestamps(1);
end_time = timestamps(end);
time_range = end_time - start_time;

% Decide on step size based on range
if time_range < 8
   % step_size = 1;  % Use 1 ms steps if the range is too small
   disp('not enough data/too short interval/trial');
else
    step_size = 8;  % Otherwise, use 8 ms steps
end

% Generate new timestamps for interpolation
new_timestamps = start_time:step_size:end_time;

% Display the size and contents of new_timestamps to confirm
disp(['new_timestamps has ', num2str(length(new_timestamps)), ' elements.']);

% Interpolate the mean gaze coordinates to match the new temporal resolution
x_mean_interp = interp1(timestamps, x_mean, new_timestamps, 'linear');
y_mean_interp = interp1(timestamps, y_mean, new_timestamps, 'linear');


%% plot interpolated gaze data
figure; hold on;
plot(x_mean, y_mean, 'bo-', 'MarkerSize', 3, 'DisplayName', 'Original Gaze Path');
plot(x_mean_interp, y_mean_interp, 'r.-', 'MarkerSize', 10, 'DisplayName', 'Interpolated Gaze Path');
xlabel('Horizontal Gaze Position (pixels)');
ylabel('Vertical Gaze Position (pixels)');
title('Interpolated Gaze Path Visualization');
legend('Original', 'Interpolated');
grid on; hold off;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Saccades (try 2)
% size of stimulation screen for plotting
%% TODO make loop instead of placeholder
trial = 1;
screenXpixels = 3240;
screenYpixels = 2160;

% 1. Preprocessing the Eye-Tracking Data
% Interpolate NaN values for left eye
x_left = data(1).eyeTrial(trial).left.gazePoint.onDisplayArea(1, :) * screenXpixels;
y_left = screenYpixels -  data(1).eyeTrial(trial).left.gazePoint.onDisplayArea(2, :) * screenYpixels; % Inverting y-axis
x_left_interpolated = fillmissing(x_left, 'linear');
y_left_interpolated = fillmissing(y_left, 'linear');

% Interpolate NaN values for right eye
x_right = data(1).eyeTrial(trial).right.gazePoint.onDisplayArea(1, :) * screenXpixels;
y_right = screenYpixels -  data(1).eyeTrial(trial).right.gazePoint.onDisplayArea(2, :) * screenYpixels; % Inverting y-axis
x_right_interpolated = fillmissing(x_right, 'linear');
y_right_interpolated = fillmissing(y_right, 'linear');

% mean gaze coordinates between the left and right eyes
x_mean = mean([x_left_interpolated; x_right_interpolated], 1);
y_mean = mean([y_left_interpolated; y_right_interpolated], 1);

% interpolate to double the data
%timestamps = double(data(1).eyeTrial(1).systemTimeStamp) / 1000000; % Convert to milliseconds
timestamps = double(data(1).eyeTrial(trial).systemTimeStamp) / 1000;
% Generate new timestamps that are halfway between each of the original timestamps
new_timestamps = linspace(timestamps(1), timestamps(end), 2 * length(timestamps) - 1);
% Interpolate the gaze coordinates to match the new temporal resolution
x_mean_interp = interp1(timestamps, x_mean, new_timestamps, 'linear');
y_mean_interp = interp1(timestamps, y_mean, new_timestamps, 'linear');







%% 2. Saccade detection
% Define a velocity threshold for saccade detection
% For 30 Hz data, a good detection threshold might be around 30°/s to 50°/s.
% Saccades at this sampling rate will have larger apparent velocities due to fewer data points representing each movement.
threshold_vel = 100;
threshold_dist = 30;


% --- smooth vs non-smoothed data (unused) --- 
% timestamps = double(data(1).eyeTrial(1).systemTimeStamp) / 1000;
% step_size = 8; threshold_vel = 100; threshold_dist = 30;
% % --- unsmoothed data ---
% [x_mean_interp, y_mean_interp, new_timestamps] = interpolateGazeData(x_mean, y_mean, timestamps, step_size);
% % Original data saccades detection
% [saccade_onsets_dist, saccade_offsets_dist] = detectSaccades(x_mean_interp, y_mean_interp, new_timestamps, 'distance', threshold_dist);
% [saccade_onsets_vel, saccade_offsets_vel] = detectSaccades(x_mean_interp, y_mean_interp, new_timestamps, 'velocity', threshold_vel);
% % --- smoothed data ---
% window_size = 5;  % Define the Gaussian window size
% [x_mean_smoothed, y_mean_smoothed] = smoothGazeData(x_mean, y_mean, window_size);
% % smoothed and interpolated
% [x_mean_smooth_interp, y_smooth_mean_interp, new_smooth_timestamps] = interpolateGazeData(x_mean_smoothed, y_mean_smoothed, timestamps, step_size);
% % Smoothed data saccades detection
% [saccade_onsets_dist_smoothed, saccade_offsets_dist_smoothed] = detectSaccades(x_mean_smooth_interp, y_smooth_mean_interp, new_smooth_timestamps, 'distance', threshold_dist);
% [saccade_onsets_vel_smoothed, saccade_offsets_vel_smoothed] = detectSaccades(x_mean_smooth_interp, y_smooth_mean_interp, new_smooth_timestamps, 'velocity', threshold_vel);
% figure;
% subplot(2, 2, 1); % Distance-Based Saccades (Original Data) 
% plotSaccades(x_mean_interp, y_mean_interp, saccade_onsets_dist, saccade_offsets_dist, sprintf('Original Data: Distance-Based Saccades (Threshold: %d pixels)', threshold_dist));
% subplot(2, 2, 2); % Velocity-Based Saccades (Original Data)
% plotSaccades(x_mean_interp, y_mean_interp, saccade_onsets_vel, saccade_offsets_vel, sprintf('Original Data: Velocity-Based Saccades (Threshold: %d deg/s)', threshold_vel));
% subplot(2, 2, 3); % Distance-Based Saccades (Smoothed Data)
% plotSaccades(x_mean_smooth_interp, y_smooth_mean_interp, saccade_onsets_dist_smoothed, saccade_offsets_dist_smoothed, sprintf('Smoothed Data: Distance-Based Saccades (Threshold: %d pixels)', threshold_dist));
% subplot(2, 2, 4); % Velocity-Based Saccades (Smoothed Data)
% plotSaccades(x_mean_smooth_interp, y_smooth_mean_interp, saccade_onsets_vel_smoothed, saccade_offsets_vel_smoothed, sprintf('Smoothed Data: Velocity-Based Saccades (Threshold: %d deg/s)', threshold_vel));
% sgtitle('Comparison of Saccade Detection Methods for Original and Smoothed Data'); saveas(gcf, 'saccade_comparison_smoothed_plot.png');  % Save the figure as a PNG file


%% --- distance and velocity based saccades --- 

% --- Distances and velocities for non-interpolated data ---
dx = diff(x_mean);  % Difference in x-coordinates (non-interpolated)
dy = diff(y_mean);  % Difference in y-coordinates (non-interpolated)
distance_non_interp = sqrt(dx.^2 + dy.^2);  % Euclidean distance

dt_non_interp = diff(timestamps) / 1000;  % Convert time differences to seconds
velocity_non_interp = distance_non_interp ./ dt_non_interp;  % Velocity

% --- Distances and velocities for interpolated data ---
dx_interp = diff(x_mean_interp);  % Difference in x-coordinates (interpolated)
dy_interp = diff(y_mean_interp);  % Difference in y-coordinates (interpolated)
distance_interp = sqrt(dx_interp.^2 + dy_interp.^2);  % Euclidean distance

dt_interp = diff(new_timestamps) / 1000;  % Convert time differences to seconds
velocity_interp = distance_interp ./ dt_interp;  % Velocity



% ---  Create histograms and log y-axis ---
figure;
% Subplot 1: Histogram of Distances (Non-Interpolated Data)
subplot(2, 2, 1);
histogram(distance_non_interp, 'Normalization', 'probability');
xlabel('Distance (pixels)');
ylabel('Probability');
title('Non-Interpolated Data: Distance');
set(gca, 'YScale', 'log');  % Set y-axis to logarithmic scale
% Subplot 2: Histogram of Velocities (Non-Interpolated Data)
subplot(2, 2, 2);
histogram(velocity_non_interp, 'Normalization', 'probability');
xlabel('Velocity (pixels/second)');
ylabel('Probability');
title('Non-Interpolated Data: Velocity');
set(gca, 'YScale', 'log');  % Set y-axis to logarithmic scale
% Subplot 3: Histogram of Distances (Interpolated Data)
subplot(2, 2, 3);
histogram(distance_interp, 'Normalization', 'probability');
xlabel('Distance (pixels)');
ylabel('Probability');
title('Interpolated Data: Distance');
set(gca, 'YScale', 'log');  % Set y-axis to logarithmic scale
% Subplot 4: Histogram of Velocities (Interpolated Data)
subplot(2, 2, 4);
histogram(velocity_interp, 'Normalization', 'probability');
xlabel('Velocity (pixels/second)');
ylabel('Probability');
title('Interpolated Data: Velocity');
set(gca, 'YScale', 'log');  % Set y-axis to logarithmic scale
sgtitle('Histograms of Distances and Velocities');
saveas(gcf, 'histogram_comparison_plot_log_y.png');  % Save the figure as a PNG file



%% get threshholds for distances/velocities
% --- Calculate the KDE/global minimum thresholds ---
% distributions are bimodal in both cases (slow non saccades vs fast saccades)
% find global minimum and use it as threshhold for detecting saccades
smoothness = 0.4;
smoothness_dist = 0.1;
% Non-interpolated data - Distance
% Fit a Kernel Density Estimate (KDE)
[dni_f, dni_xi] = ksdensity(distance_non_interp, 'Bandwidth', smoothness_dist); % Adjust bandwidth for smoothness
% Find the global minimum in the KDE
[dni_global_min_value, dni_idx] = min(dni_f); % Find the minimum value of f and its index
threshold_dist_non_interp = dni_xi(dni_idx);  % The corresponding x value (threshold)

% Non-interpolated data - Velocity
[vni_f, vni_xi] = ksdensity(velocity_non_interp, 'Bandwidth', smoothness); 
[vni_global_min_value, vni_idx] = min(vni_f);  % Find the global minimum in the KDE
threshold_vel_non_interp = vni_xi(vni_idx);  % The corresponding x value (threshold)

% Interpolated data - Distance
[di_f, di_xi] = ksdensity(distance_interp, 'Bandwidth', smoothness_dist); 
[di_global_min_value, di_idx] = min(di_f);  % Find the global minimum in the KDE
threshold_dist_interp = di_xi(di_idx);  % The corresponding x value (threshold)

% Interpolated data - Velocity
[vi_f, vi_xi] = ksdensity(velocity_interp, 'Bandwidth', smoothness); 
[vi_global_min_value, vi_idx] = min(vi_f);  % Find the global minimum in the KDE
threshold_vel_interp = vi_xi(vi_idx);  % The corresponding x value (threshold)

fprintf('Non-Interpolated Data: Global Minimum Thresholds:\nDistance: %.2f pixels\nVelocity: %.2f pixels/second\n', threshold_dist_non_interp, threshold_vel_non_interp);
fprintf('Interpolated Data: Global Minimum Thresholds:\nDistance: %.2f pixels\nVelocity: %.2f pixels/second\n', threshold_dist_interp, threshold_vel_interp);

% 2x2 plot: KDEs and the global minimum
plotKDEWithGlobalMin(dni_f, dni_xi, threshold_dist_non_interp, dni_global_min_value, ...
                     vni_f, vni_xi, threshold_vel_non_interp, vni_global_min_value, ...
                     di_f, di_xi, threshold_dist_interp, di_global_min_value, ...
                     vi_f, vi_xi, threshold_vel_interp, vi_global_min_value, ...
                     safe, comparison_results_folder, trial);
% TODO change folder above to participant folder or sth 


% --- Detect saccades using calculated thresholds ---
% Non-interpolated data
[saccade_onsets_dist, saccade_offsets_dist] = detectSaccades(dx, dy, timestamps, 'distance', threshold_dist_non_interp);
[saccade_onsets_vel, saccade_offsets_vel] = detectSaccades(dx, dy, timestamps, 'velocity', threshold_vel_non_interp );

% Interpolated data
[saccade_onsets_dist_interp, saccade_offsets_dist_interp] = detectSaccades(dx_interp, dy_interp, new_timestamps, 'distance', threshold_dist_interp);
[saccade_onsets_vel_interp, saccade_offsets_vel_interp] = detectSaccades(dx_interp, dy_interp, new_timestamps, 'velocity', threshold_vel_interp);

% --- Step 6: Plot saccades detected for both datasets ---
figure;
% Subplot 1: Distance-Based Saccades (Non-Interpolated Data)
subplot(2, 2, 1);
plotSaccades(x_mean, y_mean, saccade_onsets_dist, saccade_offsets_dist, sprintf('Non-Interpolated Data: Distance-Based Saccades (Threshold: %.2f pixels)', threshold_dist));
% Subplot 2: Velocity-Based Saccades (Non-Interpolated Data)
subplot(2, 2, 2);
plotSaccades(x_mean, y_mean, saccade_onsets_vel, saccade_offsets_vel, sprintf('Non-Interpolated Data: Velocity-Based Saccades (Threshold: %.2f deg/s)', threshold_vel));
% Subplot 3: Distance-Based Saccades (Interpolated Data)
subplot(2, 2, 3);
plotSaccades(x_mean_interp, y_mean_interp, saccade_onsets_dist_interp, saccade_offsets_dist_interp, sprintf('Interpolated Data: Distance-Based Saccades (Threshold: %.2f pixels)', threshold_dist_interp));
% Subplot 4: Velocity-Based Saccades (Interpolated Data)
subplot(2, 2, 4);
plotSaccades(x_mean_interp, y_mean_interp, saccade_onsets_vel_interp, saccade_offsets_vel_interp, sprintf('Interpolated Data: Velocity-Based Saccades (Threshold: %.2f deg/s)', threshold_vel_interp));
sgtitle('Comparison of Saccade Detection Methods (Interpolated vs Non-Interpolated Data)');
saveas(gcf, 'saccade_comparison_plot.png');  % Save the figure as a PNG file







%% --- Step 4: Calculate the 95th percentile thresholds ---
% Interpolated data
%perc = 85;
%threshold_dist_interp = prctile(distance_interp, perc);
%threshold_vel_interp = prctile(velocity_interp, perc);
threshold_dist_interp = 200; threshold_vel_interp  = 3000;
fprintf('Interpolated Data: %.2ft h Percentile Thresholds:\nDistance: %.2f pixels\nVelocity: %.2f pixels/second\n', perc, threshold_dist_interp, threshold_vel_interp);

% Non-interpolated data
%threshold_dist = prctile(distance_non_interp, perc);
%threshold_vel = prctile(velocity_non_interp, perc);
threshold_dist = 400; threshold_vel  = 3000;
fprintf('Non-Interpolated Data: %.2f th Percentile Thresholds:\nDistance: %.2f pixels\nVelocity: %.2f pixels/second\n', perc, threshold_dist, threshold_vel);

% --- Step 5: Detect saccades using calculated thresholds ---
% Non-interpolated data
[saccade_onsets_dist, saccade_offsets_dist] = detectSaccades(x_mean, y_mean, timestamps, 'distance', threshold_dist);
[saccade_onsets_vel, saccade_offsets_vel] = detectSaccades(x_mean, y_mean, timestamps, 'velocity', threshold_vel);

% Interpolated data
[saccade_onsets_dist_interp, saccade_offsets_dist_interp] = detectSaccades(x_mean_interp, y_mean_interp, new_timestamps, 'distance', threshold_dist_interp);
[saccade_onsets_vel_interp, saccade_offsets_vel_interp] = detectSaccades(x_mean_interp, y_mean_interp, new_timestamps, 'velocity', threshold_vel_interp);

% --- Step 6: Plot saccades detected for both datasets ---
figure;
% Subplot 1: Distance-Based Saccades (Non-Interpolated Data)
subplot(2, 2, 1);
plotSaccades(x_mean, y_mean, saccade_onsets_dist, saccade_offsets_dist, sprintf('Non-Interpolated Data: Distance-Based Saccades (Threshold: %.2f pixels)', threshold_dist));
% Subplot 2: Velocity-Based Saccades (Non-Interpolated Data)
subplot(2, 2, 2);
plotSaccades(x_mean, y_mean, saccade_onsets_vel, saccade_offsets_vel, sprintf('Non-Interpolated Data: Velocity-Based Saccades (Threshold: %.2f deg/s)', threshold_vel));
% Subplot 3: Distance-Based Saccades (Interpolated Data)
subplot(2, 2, 3);
plotSaccades(x_mean_interp, y_mean_interp, saccade_onsets_dist_interp, saccade_offsets_dist_interp, sprintf('Interpolated Data: Distance-Based Saccades (Threshold: %.2f pixels)', threshold_dist_interp));
% Subplot 4: Velocity-Based Saccades (Interpolated Data)
subplot(2, 2, 4);
plotSaccades(x_mean_interp, y_mean_interp, saccade_onsets_vel_interp, saccade_offsets_vel_interp, sprintf('Interpolated Data: Velocity-Based Saccades (Threshold: %.2f deg/s)', threshold_vel_interp));
sgtitle('Comparison of Saccade Detection Methods (Interpolated vs Non-Interpolated Data)');
saveas(gcf, 'saccade_comparison_plot.png');  % Save the figure as a PNG file


