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
    };

% Define the classification of the folders (ADHD or non-ADHD)
group_labels = {
    'nonADHD';  % alex
    'ADHD';     % mara 
    'ADHD';     % tilo 
    'ADHD';     % anu 
    'ADHD';     % kieran 
    'ADHD';     % sura 
    'nonADHD';  % hamit
    'nonADHD';  % jannik
    'ADHD';     % farn 
    'nonADHD';  % julia 
    'ADHD';     % felix 
    'nonADHD';  % florian 
    };

conditions = {'a', 'a_simple', 'b', 'b_simple'};
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
    eye_rt_file = dir(strcat(folder, '\analysis\eye_rt.mat')); % calculated eye rt from analyseResults.m
    eyetracking_file = dir(fullfile(folder, '\analysis\cut_trials*.mat'));
    pupil_file = dir(strcat(folder, '\analysis\pupil_before_after_finding_stim.mat'));
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
        load(fullfile(eye_rt_file.folder, eye_rt_file.name), 'eye_rt');
        data_struct(i).rt_right          = eye_rt.RightEyeRT;
        data_struct(i).rt_left           = eye_rt.LeftEyeRT;
        % eyetracking data cut into trials
        load(fullfile(eyetracking_file.folder, eyetracking_file.name), 'cutData');
        data_struct(i).eyeTrial           = cutData.eyeTrial;
        data_struct(i).stimulusTrial      = cutData.stimulusTrial; % cut to start stim
        data_struct(i).TargetPosition     = cutData.TargetPosition;
        data_struct(i).x_centers          = cutData.x_centers;
        data_struct(i).y_centers          = cutData.y_centers;
        % pupil dilation data
        load(fullfile(pupil_file.folder, pupil_file.name), 'result_table');
        data_struct(i).pupilDiam          = result_table;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS/ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- EYE-TRACKING DATA ---


%% analyse saccades (unzufriedenstellend)
% sampling rate too low to really get the saccades this way - only kind of works
%analyseSaccades(data, screenXpixels, screenYpixels, safe, comparison_results_folder)

%% --- detect fixations ---
safe = 1;
% params to define fixations
dist_threshold = 50; % max distance between consecutive points (in pixels)
min_duration = 3; % min num of consecutive points to count as a fixation

all_fixations = struct();
for participant=1:size(data,2)
    analysis_folder = strcat(folders{participant}, '\analysis\'); % save per participant
    fixations = analyseFixation(data(participant), dist_threshold, min_duration, screenXpixels, screenYpixels, safe, analysis_folder);
    all_fixations(participant).id        = data(participant).id;    
    all_fixations(participant).fixations = fixations;    
end
save(strcat(comparison_results_folder, 'all_fixations.mat'), 'all_fixations'); % save for all participants

% plot some fixations for checking 
plot_these = [1, 160]; % just plot two for checking
for participant = 1:size(data, 2) % size(plot_these, 2)
    participant_data = data(participant);
    participant_fixations = all_fixations(participant).fixations; %.stimulusFixations;
    analysis_folder = strcat(folders{participant}, '\analysis\'); 
    plotFixations(participant_data, participant_fixations, screenXpixels, screenYpixels, plot_these, analysis_folder, safe);
    close all
end

% get fixation durations, counts per trial and average over conditions
fixation_stats = getFixationDuration(all_fixations, data, comparison_results_folder);


%% --- reverse engineer "saccades" from fixation clusters ---
all_saccades = struct();
for participant = 1:size(data, 2)
    saccades = struct();
    for trial = 1:size(data(participant).stimulusTrial, 1)
        saccades(trial).trial = trial;
        participant_data = data(participant).stimulusTrial(trial);
        participant_fixations = all_fixations(participant).fixations.stimulusFixations(trial).fixations; %.stimulusFixations;
        saccades(trial).saccades = detectSaccadesFromFixations(participant_fixations, screenXpixels, screenYpixels);
        saccades(trial).saccadeDurations = getSaccadeTimes(participant_fixations, participant_data);
    end
    all_saccades(participant).id = data(participant).id; 
    all_saccades(participant).saccades = saccades;
end
save(fullfile(comparison_results_folder, 'all_saccades.mat'), 'saccades');% save the saccades data

% plot saccades
plot_these = [1,160];
for participant = 1:size(data, 2)
    analysis_folder = strcat(folders{participant}, '\analysis\'); 
    participant_data = data(participant).stimulusTrial;
    participant_fixations = all_fixations(participant).fixations.stimulusFixations;
    plotFixationSaccades(data(participant), all_saccades(participant), all_fixations(participant), screenXpixels, screenYpixels, plot_these, analysis_folder, safe);
end
close all

%% --- get distance from starting point in trial to target center --- 
for participant = 1:size(data, 2)
    for trial=1:size(data(participant).stimulusTrial, 1)
        participant_data = data(participant).stimulusTrial(trial);
        participant_fixations = all_fixations(participant).fixations.stimulusFixations(trial).fixations; %.stimulusFixations;
        participant_saccades = all_saccades(participant).saccades(trial).saccades;
        
        % coordinates first gaze of trial
        first_gaze_x = participant_data.left.gazePoint.onDisplayArea(1, 1) * screenXpixels;
        first_gaze_y = screenYpixels - (participant_data.left.gazePoint.onDisplayArea(2, 1) * screenYpixels);
        % coordinates target
        target_x = all_fixations(participant).fixations.targetCenters(1, trial);
        target_y = all_fixations(participant).fixations.targetCenters(2, trial);
        % number of saccades in trial 
        n_sac = size(participant_saccades, 2);
        % sum distance of saccades
        distance_sac = sum([participant_saccades.distance]);
        mean_distance_sac = mean([participant_saccades.distance]);
        % sum of entire gaze path
        %distance_trial = sum();
    end
end


%% --- Fixation Duration and Distribution ---
%fixation_stats = getFixationDuration(all_fixations, data, conditions, comparison_results_folder);
% get fixation durations, counts per trial 
fixation_stats = getFixationDuration(all_fixations, data, comparison_results_folder);
%  average over conditions
meanFixationDurations = getMeanFixationDurationPerCondition(fixation_stats, conditions);

% Fixation Duration: Compare the average duration of fixations between ADHD 
% and non-ADHD participants. ADHD individuals often have shorter fixation durations, 
% indicating less stable attention on specific areas of interest.
plotGroupFixationDurations(fixation_stats, group_labels, conditions, color_map);% line
plotFixationDifferences(fixation_stats, group_labels, conditions, color_map); % line
plotViolinFixationStats(fixation_stats, group_labels, conditions, color_map, comparison_results_folder, safe);

% Fixation Distribution: Analyze the spatial distribution of fixations. 
% ADHD participants might have a more scattered fixation pattern, with less 
% focus on task-relevant areas compared to non-ADHD participants, 
% who might show more concentrated fixations on targets.
% this just plots mean of all the fixations per condition (useless)
%plotFixationSpatialDistribution(fixationStats, group_labels, conditions, comparison_results_folder, safe);


%% trash
% % distance of first gaze trial to targer
% trial_distance = sqrt((endCenter(1) - startCenter(1)).^2 + (endCenter(2) - startCenter(2)).^2);
% % distance of all saccades in trial 
% saccade_distance = sqrt((endCenter(1) - startCenter(1)).^2 + (endCenter(2) - startCenter(2)).^2);
% % total distance of gaze in trial
% total_distance = sqrt((endCenter(1) - startCenter(1)).^2 + (endCenter(2) - startCenter(2)).^2);




%% --- Saccade Metrics ---
% Saccade Amplitude and Velocity: ADHD participants may exhibit larger saccade 
% amplitudes and higher velocities, reflecting more erratic or less controlled eye movements.
saccadeStats = analyzeSaccadeAmplitudeAndVelocity(all_saccades, data, conditions);
plotSaccadeDifferences(saccadeStats, group_labels, conditions, color_map);
% same but as violin
plotSaccadesViolin(saccadeStats, group_labels, conditions, color_map, comparison_results_folder, safe)

% Saccade Frequency: Higher saccade frequency might indicate more frequent 
% shifts in attention, which is often characteristic of ADHD.
plotSaccadeCount(saccadeStats, group_labels, conditions, color_map);

ttestFixationDuration(group_labels, fixation_stats, saccadeStats, conditions)


%% --- Pupil Diameter ---
% check target coords/eye coords (+ tolerance)
% get 30 datapoints before target is reached with eyes + 10 datapoints after
% padded with NaNs if not enough datapoints
% these 3 are set in analyseResults.m just fyi here:
% tolerance = 100;  % Tolerance in pixels for detecting gaze near the target
num_before = 30; % Number of datapoints before target found
num_after = 20;  % Number of datapoints after target found

% per participant: normalize pupil diameter by baseline: mean and median during blank screen
% each trial preceding the stimulation
% returns: 
normalized_data = calcProportionalPupilChange(data);

% get averages 
[norm_mean_diam_around_stim_adhd, norm_mean_diam_around_stim_nonadhd] = getPupilDiamsPerGroup(normalized_data, 'mean', group_labels, conditions);
[norm_median_diam_around_stim_adhd, norm_median_diam_around_stim_nonadhd] = getPupilDiamsPerGroup(normalized_data, 'median', group_labels, conditions);

% plot 30 datapoints before target was reached with gaze and 10 after
% average per participant/group/condition
plotAvgPupilDiamsBeforeAfterStimFoundByGroupByCondition('Mean', norm_mean_diam_around_stim_adhd, norm_mean_diam_around_stim_nonadhd, conditions, color_map, num_before, num_after, comparison_results_folder);
plotAvgPupilDiamsBeforeAfterStimFoundByGroupByCondition('Median', norm_median_diam_around_stim_adhd, norm_median_diam_around_stim_nonadhd, conditions, color_map, num_before, num_after, comparison_results_folder);







%% --- Search Efficiency ---
% Path Efficiency: Measure the efficiency of the visual search path. 
% Non-ADHD participants might display more direct paths towards the target, 
% whereas ADHD participants could show more erratic paths with unnecessary movements.


% Time to First Fixation on Target: Assess the time taken for the first 
% fixation on the target. Delays in this metric could suggest differences 
% in how quickly participants orient their attention to relevant stimuli.


%% --- Search Strategy Analysis ---
% Systematic vs. Random Search Patterns: Assess whether participants use a systematic 
% search strategy (e.g., scanning in a structured way) versus a more random or 
% chaotic approach. ADHD participants might be less likely to use systematic strategies.

% Search Latency: Analyze the latency before initiating a search, and how quickly 
% participants give up on a search (indicative of persistence and focus).
searchLatencyStats = analyzeSearchLatency(all_fixations, data, conditions, comparison_results_folder);






%% --- Attention Shifts and Disengagement ---
% Task Engagement Over Time: Evaluate how task engagement changes over time. 
% ADHD participants might show a decline in performance or more frequent disengagement from the task as it progresses.

% Revisits: Track how often participants revisit previously searched areas. 
% Frequent revisits might indicate memory deficits or difficulties in sustaining attention.





























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REACTION TIME AND ACCURACY
%% Compare RT_ButtonPress and RT_Eye across the four conditions ("a", "a simple", "b", "b simple").
%% Hyp: Easier conditions ("a simple", "b simple") might have faster RTs than harder conditions ("a", "b").
%% Hyp: ADHD participants might show different reaction time patterns compared
%%      to non-ADHD participants, especially in harder conditions.

%% Compare Reaction Times Across Conditions and Groups
% H: ADHD participants might show different reaction time patterns compared
% to non-ADHD participants, especially in harder conditions


%% RT variability
% TODO Add y label einheit, vairnace, seconds squared/make more clear
% Use individual colours per person

measure = 'var';  % var/std
rt_variability = getRTvariability(data, conditions, groups, measure, color_map, comparison_results_folder);
rt_variability_2 = getAndPlotRTV_ParticipantLines(data, conditions, measure, groups, color_map, comparison_results_folder);


%% correlation analysis: RTeye, RTbuttonpress and accuracy
[R, P] = correlationRTbuttonRTeyeAccuracy(data, comparison_results_folder);


%% condition-Specific Correlation Analysis
[R_condition, P_condition] = correlationByCondition(data, conditions, comparison_results_folder);

% correlation has nans in b_simple since theres only 1 incorrect trial for
% all participants/trials
checkAccuracyDistributionByCondition(data, conditions);


%% ADHD vs. non-ADHD Group Correlation Analysis
[R_group, P_group] = correlationByGroup(data, groups, comparison_results_folder);

%% BOTH TODO look this up
[R_group_condition, P_group_condition] = correlationByGroupCondition(data, groups, conditions, comparison_results_folder);
% TODO corr of groups by condition
% TODO add Rtbutton-rteye (lapse?) 
plotGroupConditionCorrelations(R_group_condition, groups, conditions)
compareGroupCorrelations(R_group_condition, groups, data, conditions)



%% BAR and LINE plot: mean RT and error bars (~variability) per group and condition
% plotBarSEMMeanRTConditonGroup(groups, conditions, data, color_map, safe, comparison_results_folder)

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
% TODO Raw data: Ratio between a and b (lapse time between a and b; a/b)
plotMeanRTButtonPressVsEyecomparison(data, color_map, comparison_results_folder, safe)

%% Confusion/ RTa - RTb or RTa/RTb
% TODO add errorbars, add button press,, Do for eye rt, ttest
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

% Merge questionnaire data with RT/accuracy data based on ID
merged_data = innerjoin(data_table, quest_table, 'Keys', 'id');
merged_data = innerjoin(merged_data, quest_scores, 'Keys', 'id');
merged_data = innerjoin(merged_data, struct2table(rt_variability_2), 'Keys', 'id');


%% define attentional levels based on ASRS guidelines
% define thresholds for categorizing attention levels (based on guidelines)
thresholds = [0.5, 3.5, 4.5, 5.5];  % Thresholds between categories
attention_levels = discretize(merged_data.PartASymptomsCount, thresholds, 'categorical', {'Low', 'Moderate', 'High'});

% Add the attention level to your data
merged_data.InattentionLevel = attention_levels;

% Initialize empty arrays to hold mean RTs and accuracies
mean_rt_adhd = [];
mean_accuracy_adhd = [];
mean_rt_nonadhd = [];
mean_rt_a_adhd = [];
mean_rt_b_adhd = [];
mean_rt_asimple_adhd = [];
mean_rt_bsimple_adhd = [];
mean_rt_a_nonadhd = [];
mean_rt_b_nonadhd = [];
mean_rt_asimple_nonadhd = [];
mean_rt_bsimple_nonadhd = [];
mean_accuracy_nonadhd = [];
mean_6inatt_adhd = [];
count_6inatt_adhd = [];
mean_allinatt_adhd = [];
mean_inatt_nonadhd = [];
count_6inatt_nonadhd = [];
mean_allinatt_nonadhd = [];
mean_con_adhd = [];
mean_con_nonadhd = [];
rtv_a_adhd = [];
rtv_b_adhd = [];
rtv_asimple_adhd = [];
rtv_bsimple_adhd = [];
rtv_a_nonadhd = [];
rtv_b_nonadhd = [];
rtv_asimple_nonadhd = [];
rtv_bsimple_nonadhd = [];


% Loop through each participant and calculate mean RT and accuracy
for i = 1:size(merged_data,1)
    % Check if the participant belongs to the ADHD group
    if strcmp(merged_data.group_merged_data{i}, 'ADHD')
        % the mean RT
        mean_rt_a_adhd = [mean_rt_a_adhd; merged_data.nRTa(i)];
        mean_rt_b_adhd = [mean_rt_b_adhd; merged_data.nRTb(i)];
        mean_rt_asimple_adhd = [mean_rt_asimple_adhd; merged_data.nRTasimple(i)];
        mean_rt_bsimple_adhd = [mean_rt_bsimple_adhd; merged_data.nRTbsimple(i)];
        % rtv
        rtv_a_adhd = [rtv_a_adhd; merged_data.a(i).RT_ButtonPress_var];
        rtv_b_adhd = [rtv_b_adhd; merged_data.b(i).RT_ButtonPress_var];
        rtv_asimple_adhd = [rtv_asimple_adhd; merged_data.a_simple(i).RT_ButtonPress_var];
        rtv_bsimple_adhd = [rtv_bsimple_adhd; merged_data.b_simple(i).RT_ButtonPress_var];
        % mean accuracy
        mean_accuracy_adhd = [mean_accuracy_adhd; mean(merged_data.accuracy{i}, 'omitnan')];
        % confusion
        mean_con_adhd  = [mean_con_adhd; merged_data.nRTa(i)/merged_data.nRTb(i)];
        % inattention level
        mean_6inatt_adhd  = [mean_6inatt_adhd; merged_data.PartAMeanSymptoms(i)];
        count_6inatt_adhd = [count_6inatt_adhd; merged_data.PartASymptomsCount(i)];
        mean_allinatt_adhd = [mean_allinatt_adhd; merged_data.MeanSymptoms(i)];      
        
    % Check if the participant belongs to the non-ADHD group
    elseif strcmp(merged_data.group_merged_data{i}, 'nonADHD')
        % mean RT
        %rt_values = [merged_data.nRTa(i), merged_data.nRTb(i), merged_data.nRTasimple(i), merged_data.nRTbsimple(i)];
        %mean_rt_nonadhd = [mean_rt_nonadhd; rt_values]; % mean(rt_values, 'omitnan')];
        mean_rt_a_nonadhd = [mean_rt_a_nonadhd; merged_data.nRTa(i)];
        mean_rt_b_nonadhd = [mean_rt_b_nonadhd; merged_data.nRTb(i)];
        mean_rt_asimple_nonadhd = [mean_rt_asimple_nonadhd; merged_data.nRTasimple(i)];
        mean_rt_bsimple_nonadhd = [mean_rt_bsimple_nonadhd; merged_data.nRTbsimple(i)];
        % rtv
        rtv_a_nonadhd = [rtv_a_nonadhd; merged_data.a(i).RT_ButtonPress_var];
        rtv_b_nonadhd = [rtv_b_nonadhd; merged_data.b(i).RT_ButtonPress_var];
        rtv_asimple_nonadhd = [rtv_asimple_nonadhd; merged_data.a_simple(i).RT_ButtonPress_var];
        rtv_bsimple_nonadhd = [rtv_bsimple_nonadhd; merged_data.b_simple(i).RT_ButtonPress_var];
        % confusion
        mean_con_nonadhd  = [mean_con_nonadhd; merged_data.nRTa(i)/merged_data.nRTb(i)];
        % mean accuracy
        mean_accuracy_nonadhd = [mean_accuracy_nonadhd; mean(merged_data.accuracy{i}, 'omitnan')];
        % inattention level
        mean_inatt_nonadhd = [mean_inatt_nonadhd; merged_data.PartAMeanSymptoms(i)];
        count_6inatt_nonadhd = [count_6inatt_nonadhd; merged_data.PartASymptomsCount(i)];
        mean_allinatt_nonadhd = [mean_allinatt_nonadhd; merged_data.MeanSymptoms(i)];
    end
end


%% ttest 
ttestResults = runTTests_ADHD_vs_nonADHD(...
    mean_rt_a_adhd, mean_rt_b_adhd, ...
    mean_rt_asimple_adhd, mean_rt_bsimple_adhd, ...
    mean_rt_a_nonadhd, mean_rt_b_nonadhd, ...
    mean_rt_asimple_nonadhd, mean_rt_bsimple_nonadhd, ...
    rtv_a_adhd, rtv_a_nonadhd,...
    rtv_b_adhd,rtv_b_nonadhd, ... 
    rtv_asimple_adhd,rtv_asimple_nonadhd,... 
    rtv_bsimple_adhd,rtv_bsimple_nonadhd, ...
    mean_accuracy_adhd, mean_accuracy_nonadhd, ...
    mean_con_adhd, mean_con_nonadhd, ...
    mean_6inatt_adhd, mean_inatt_nonadhd, ...
    count_6inatt_adhd, count_6inatt_nonadhd, ...
    mean_allinatt_adhd, mean_allinatt_nonadhd);


%% Correlation Analysis
% Plot both correlations

% Correlation Analysis for BOTH group
% [r_adhd_rt, p_adhd_rt] = corr(mean_inatt_adhd, [mean_rt_a_adhd, mean_rt_b_adhd, mean_rt_asimple_adhd, mean_rt_bsimple_adhd]);
% [r_adhd_rtv, p_adhd_rtv] = corr(mean_inatt_adhd, [rtv_a_adhd, rtv_b_adhd, rtv_asimple_adhd, rtv_bsimple_adhd]);
% [r_adhd_con, p_adhd_con] = corr(mean_inatt_adhd, mean_con_adhd);
% [r_adhd_acc, p_adhd_acc] = corr(mean_inatt_adhd, mean_accuracy_adhd);
% 
% % Display results for BOTH group
% fprintf('Correlation between inattention levels and RT: r = %.2f, p = %.3f\n', r_adhd_rt, p_adhd_rt);
% fprintf('Correlation between inattention levels and RTV: r = %.2f, p = %.3f\n', r_adhd_rtv, p_adhd_rtv);
% fprintf('Correlation between inattention levels and Accuracy: r = %.2f, p = %.3f\n', r_adhd_acc, p_adhd_acc);
% fprintf('Correlation between inattention levels and Confusion: r = %.2f, p = %.3f\n', r_adhd_con, p_adhd_con);





