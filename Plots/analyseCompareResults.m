function analyseCompareResults(color_map, color_map_other, conditions, condition_labels, groups, ...
                              screenXpixels, screenYpixels, sr, safe, ...
                              path_to_folders, valid_participant_folders, ids, analysis_subfolder, comparison_results_folder, ...
                              subfolder_fixation, fixation_duration)


data_struct = struct();
for i = 1:length(valid_participant_folders) % Loop through each folder and load the data
    folder = strcat(path_to_folders, valid_participant_folders{i});
    data_file   = dir(fullfile(folder, '\results\trial_results*.mat')); % original data
    eyetracking_file = dir(fullfile(folder, analysis_subfolder, '\cut_trials*.mat'));
    load(strcat(folder, analysis_subfolder, '\eye_rt.mat')); % calculated eye rt from analyseResults.m
    %pupil_file = dir(strcat(folder, analysis_subfolder, '\pupil_20before_31after_finding_stim.mat'));
    %pupil_diam_norm_file = dir(strcat(folder, analysis_subfolder, '\pupil_processed_found_stim_StartOfTrial.mat'));
    
    participant_file = fullfile(folder, 'info.csv'); % participant infos
    info_data = readtable(participant_file);
    group = "";
    if strcmp(info_data.SubjectADHD, 'yes')
        group = 'ADHD';        
    elseif strcmp(info_data.SubjectADHD, 'no')
        group = 'nonADHD'; 
    end
    
    if ~isempty(data_file)
        load(fullfile(data_file.folder, data_file.name), 'trial_results');
        data_struct(i).id                = ids(i); %info_data.SubjectID;
        data_struct(i).group             = group;
        data_struct(i).age               = info_data.Subject_Age;
        data_struct(i).gender            = info_data.SubjectGender;
        data_struct(i).Condition         = trial_results.Condition;
        data_struct(i).rt                = trial_results.rt;
        data_struct(i).accuracy          = trial_results.correct;
        data_struct(i).StimulusOnsetTime = trial_results.StimulusOnsetTime;
        data_struct(i).trialEndTime      = trial_results.trialEndTime;
        % eye reaction times
        %load(fullfile(eye_rt_file.folder, eye_rt_file.name), 'eye_rt');
        data_struct(i).rt_right          = eye_rt.RightEyeRT;
        data_struct(i).rt_left           = eye_rt.LeftEyeRT;
        
        try% eyetracking data cut into trials
            load(fullfile(eyetracking_file.folder, eyetracking_file.name), 'cutData');
        catch % take the first you can find in case theres serveral
            load(fullfile(eyetracking_file(1).folder, eyetracking_file(1).name), 'cutData');
        end
        
        data_struct(i).eyeTrial           = cutData.eyeTrial;
        data_struct(i).stimulusTrial      = cutData.stimulusTrial; % cut to start stim
        data_struct(i).TargetPosition     = cutData.TargetPosition;
        data_struct(i).x_centers          = cutData.x_centers;
        data_struct(i).y_centers          = cutData.y_centers;
        data_struct(i).angleMatrix        = cutData.AngleMatrix;
        %         % pupil dilation data
        %         load(fullfile(pupil_file.folder, pupil_file.name), 'result_table');
        %         data_struct(i).pupilDiam          = result_table;
        %         load(fullfile(pupil_diam_norm_file.folder, pupil_diam_norm_file.name), 'condition_diam_avg');
        %         data_struct(i).pupilDiamNorm      = condition_diam_avg;
        
        % the following results depend on their fixation duration
        load(strcat(folder, '\',subfolder_fixation, '\diam_t0_ntse.mat'));
        load(strcat(folder, '\',subfolder_fixation, '\diam_t0_ntss.mat'));
        load(strcat(folder, '\',subfolder_fixation, '\diam_t0_tse.mat'));
        load(strcat(folder, '\',subfolder_fixation, '\diam_t0_tss.mat'));
        load(strcat(folder, '\',subfolder_fixation, '\diam_t0_tss_first.mat'));
        data_struct(i).ntse         = diam_t0_ntse;
        data_struct(i).ntss         = diam_t0_ntss;
        data_struct(i).tse          = diam_t0_tse;
        data_struct(i).tss          = diam_t0_tss;
        data_struct(i).tss_first   = diam_t0_tss_first;
        
        load(strcat(folder, '\',subfolder_fixation, '\trial_metrics.mat'));
        data_struct(i).trial_metrics= trial_metrics;
        
        load(strcat(folder, '\',subfolder_fixation, '\first_sacc.mat'));
        data_struct(i).first_sacc = first_sacc;
    else
        warning(['Data file not found in: ', folder]);
    end
end

group_labels = arrayfun(@(x) x.group, data_struct, 'UniformOutput', false);
group_labels = group_labels';
ids = arrayfun(@(x) x.id, data_struct);

% individual symbol and colour per participant
% set font size for all plots
set(gca, 'FontSize', 12); % axis and tick labels 
set(0, 'DefaultAxesFontSize', 12); % Applies to all new figures
set(0, 'DefaultTextFontSize', 14); % Applies to text objects

color_map_individual = getIndividualMarkersPerParticipant(group_labels, ids, comparison_results_folder);

% get the faster of the two eye RTs; ignore 0s where target was not looked at
data_struct = processEyeRTs(data_struct);


% normalized RT per participan, take mean or median for rt/normalisation
% median
average = 'median'; % 'mean' or 'median'
data_median = normalizeMeanRTsBySimpleConditions(data_struct, average, conditions, comparison_results_folder);

%mean
%average = 'mean'; % 'mean' or 'median'
%data_mean = normalizeMeanRTsBySimpleConditions(data_struct, average, conditions, comparison_results_folder);


% TODO "wrong" nontarget bar on target
%[data_all_median, data_correct_median] = wrongTargetBarWhichAndRatio(data_struct);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  PLOTS/ANALYSIS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% central vs peripheral
plotRTbyEccentricity(data_median, screenXpixels, screenYpixels, comparison_results_folder, conditions)

% get group labels and ids per participant into lists for easier plotting
getGroupStats(data_struct, group_labels);




%% ---- REACTION TIME AND ACCURACY ----
%% Reaction Time
% 3x3 subplots: rt button press vs rt eyes vs diff; means of groups per condition
%plotMeanRTButtonPressVsEyecomparison(data_median, color_map, conditions, condition_labels, comparison_results_folder, safe);

% rt over time per participant
%plotRTperParticipantOverTimeColouredByGroup(data_median, color_map, safe, comparison_results_folder)

%% RT variability
% values used here dont differ between data_median and data_mean so either can be used
measure = 'Standard Deviation';  % Use "Standard Deviation", "Coefficient of Variation", "Variance" or "Mean Squared Successive Difference" (or "Ex-Gaussian Analysis")
[rt_median_eye,rt_median_button, rt_mean_eye, rt_mean_button]  = getRTvariability(data_median, group_labels, ids, conditions, condition_labels, measure, color_map, color_map_individual, comparison_results_folder);
[rt_median_eye_2,rt_median_button_2, rt_mean_eye_2, rt_mean_button_2]  = getRTvariability(data_median, group_labels, ids, conditions, condition_labels, 'Mean Squared Successive Difference', color_map, color_map_individual, comparison_results_folder);
close all;
%% Confusion/ RTa - RTb or RTa/RTb for both button press/eye RT
% Mean confusion of participants for each condition plus p values in names
data_mean = plotConfusionPerGroup('Mean', data_mean, group_labels, ids, conditions, color_map, color_map_other,color_map_individual, comparison_results_folder, safe);
data_median = plotConfusionPerGroup('Median', data_median, group_labels, ids, conditions, color_map, color_map_other,color_map_individual, comparison_results_folder, safe);
close all;
%% Accuracy
% Accuracy per condition and group
checkAccuracyDistributionByCondition(data_median, conditions);

% Difference in correct trials; calcs median didnt change name
[mistakes, accuracy] = plotMeanAccuracyPerGroupCondition(data_median, group_labels, ids, conditions, condition_labels, color_map,color_map_individual, comparison_results_folder, safe);

% accuracy vs rt button press per condition and group
plotAccuracyVsButtonPressRT(data_median, mistakes, accuracy, rt_median_eye, rt_median_button, group_labels, ids, conditions, condition_labels, color_map,color_map_individual, comparison_results_folder, safe)

% Correlations: speed accuracy tradeoff
% info: corrs higher in simple conditions since the target was almost never fixated when mistakes were made!
correlationsSpeedAccuracyTradeOff(data_median, accuracy,rt_median_eye,rt_median_button, group_labels, groups, conditions, comparison_results_folder)

% Interaction Between Group, Condition (Difficulty), and RT
plotInteractionRTeyebuttonConditions(groups, conditions, data_median, color_map, safe, comparison_results_folder)
close all;
%% --- parallel advantage -2022 zhaoping paper plots ---
% Mean reaction times of participants for each condition
plotBarRTmeanParticipantsCondition(data_median, color_map, comparison_results_folder, color_map_individual, ids,conditions, condition_labels);
%same but per group
plotBarRTmeanParticipantsConditionPerGroup(data_median, color_map, comparison_results_folder, color_map_individual);

%% ---- EYETRACKING DATA ----
% Load fixations/saccades or calculate if nothing to load
dist_threshold = 100; % max distance between consecutive points (in pixels) for a fixations cluster
plot_these = []; % [2, 119]just plot two for checking
comp_results_fix = strcat(comparison_results_folder,subfolder_fixation, '\');
try
    mkdir(comp_results_fix);
catch % folders already exists
end

% -01 load/detect fixations -
all_fixations = analyseFixationAllWrapper(data_median, plot_these, dist_threshold, screenXpixels, screenYpixels, safe, path_to_folders ,valid_participant_folders, subfolder_fixation,comp_results_fix, fixation_duration);

% -02 load/reverse engineer "saccades" from fixation clusters -
all_saccades = analyseSaccadesAllWrapper(data_median, all_fixations, plot_these, screenXpixels, screenYpixels, safe, path_to_folders, valid_participant_folders, subfolder_fixation, comp_results_fix, fixation_duration);

% -03 get num saccades, distance, angle, diff, cosine from trial start to target center/saccades,... -
tolerance = 200; % num pixels away to count as target found
trial_metrics = calcGazeDistanceDirection(data_median, all_fixations, all_saccades, tolerance, screenXpixels, screenYpixels);

                    
%% --- 01_fixation: Fixation Duration and Distribution ---
% get fixation durations, counts per trial , avg per condition
fixation_stats = getFixationDuration(all_fixations, data_median,conditions, comp_results_fix);

% fixation duration
plotGroupFixationDurations(fixation_stats, ids, group_labels, condition_labels, conditions, color_map_individual, color_map, comp_results_fix, fixation_duration, safe);

% number of fixations
plotBarNumFixations(fixation_stats, group_labels, ids, conditions,condition_labels, color_map,color_map_individual, comp_results_fix, safe)

%% --- 01_saccade: Saccade Metrics ---
% Saccade Amplitude and Velocity
saccadeStats = analyzeSaccadeAmplitudeAndVelocity(all_saccades, data_median, conditions);
plotSaccadeDifferences(saccadeStats,data_median, group_labels, ids, conditions,condition_labels, color_map, color_map_individual, comp_results_fix, safe);
% -> looks basically the same, just as zp said it would be

    % Saccade Frequency -> is basically plotBarNumFixations since between each 2 fixations i assume a saccade
    %plotSaccadeCount(saccadeStats, group_labels, conditions, color_map, ids, color_map_individual, comp_results_fix, safe);

% Time to First Saccade
plotAverageStartFirstSaccadeGroup(data_median, ids, conditions, condition_labels, group_labels, color_map, color_map_individual, comp_results_fix, safe);

    % Search Latency -> same as avg first saccade start
    %searchLatencyStats = analyzeSearchLatency(all_fixations, data_median, conditions, comparison_results_folder);
    %nanmean(searchLatencyStats(1).conditions(1).searchLatencies)

% Arrive Abandon Return Trials
% 1. If a target fixation is followed by a non-target fixation -> AARTs
% 2. If the target fixation is the last fixation -> no abandon
% 3. combo of the two
is_abandon = 201; % 1 pixel more than counts as target found
participant_data = plotAbandonReturnTrials(is_abandon, data_struct, all_fixations, conditions, condition_labels, group_labels, color_map, comp_results_fix, safe);

% RThand - RTeye for AARTs
plotAARtrialsRThandeye(participant_data, data_struct, ids, conditions, condition_labels, group_labels, color_map, color_map_individual, comp_results_fix, safe);

% Cosine Similarity: search path/saliency
plotCosineSimilarityGroup(data_median, group_labels, ids, all_saccades, all_fixations, ...
    conditions, condition_labels, color_map, color_map_individual, comp_results_fix);

%% --- QUESTIONAIRE RESULTS --- 
% get results
[quest_struct, quest_table] = loadQuesionnaires(path_to_folders, valid_participant_folders, group_labels, ids);% load questionnaire results
% get stats table
quest_table = groupQuestionnaire(quest_table, group_labels);
% usual plots
plotQuestionnaires(ids, quest_struct, quest_table, group_labels, color_map, color_map_other, color_map_individual,comparison_results_folder, safe)


%% INFLUENCE OF (IN)ATTENTION
%% combine questionnaire results with RT and accuracy
    % Merge questionnaire data with RT/accuracy/saccade data based on ID
    %data_table = struct2table(data_median);
    %merged_data = innerjoin(data_table, quest_table, 'Keys', 'id');

% correlation rtbutton, accuracy and quest scores 
correlateQuestionnaire(conditions, rt_median_button, accuracy, quest_table)

%% PUPIL DIAMETER
% pupil diameter over all participants, grouped by adhd nonadhd and averages of both each
[diam_t0_ADHD, diam_t0_nonADHD] = plotGroupPupilDiams(data_struct,conditions, color_map, sr, condition_labels, group_labels, comparison_results_folder);



% % Normalized and Max Pupil Diameter 
% 
% % --- plot nRT versus maxPupilDiam after target reaches stim ---
% % Set threshold for outlier detection
% % 2 plots: per group and per participant
% outlier_threshold = 4; % remove outliers from rt/rtv, already done for pupil dilation in analyseResults
% [avg_adhd, avg_nonadhd] = extractPupilData(data_median, outlier_threshold, 'ADHD');
% plotPupilRTRelationships(avg_adhd, avg_nonadhd, conditions, color_map, comparison_results_folder);
% 
% % --- pupil diameter over time per group with stdabweichung ---
% % get data into format needed for plot
% outlier_threshold = 4;
% [avg_adhd, avg_nonadhd] = extractPupilData(data_median, outlier_threshold, 'ADHD');
% before = size(data_median(1).pupilDiamNorm.a.beforeMean,2);
% after = size(data_median(1).pupilDiamNorm.a.afterMean,2);
% [avg_adhd_struct, avg_nonadhd_struct] = changeFormatPupilRT(avg_adhd, avg_nonadhd, before, after);
% plotAvgPupilDiamsBeforeAfterStimFoundByGroupByCondition('Mean', 'Mean Change in Pupil Diameter (z-score)', ...
%     avg_adhd_struct, avg_nonadhd_struct, conditions, color_map, before, after, outlier_threshold, comparison_results_folder);
% avg_both = [avg_adhd; avg_nonadhd];
% 
% % --- max pupil diam after stim found per group ---
% plotPupilMaxStimFound(avg_adhd, avg_nonadhd, conditions, color_map, comparison_results_folder);
% % --- max min and mean pupil diam during stimulus presentation ---
% plotPupilMagnitudesTwoGroups(avg_adhd, avg_nonadhd, conditions, color_map, comparison_results_folder)

close all;
end