function analyseChangePupilDilation(id, cut_data, trial_metrics, conditions, condition_labels, analysis_folder, compare_folder, ...
                                    do_plots, time_vector, x_label_text, color_map)


%%  --- 4.1 plot saccades TOWARDS target aligned by START of saccade [raw, 0-aligned, cleaned] ---
diam_t0_tss = getPupilDiamsBasedOnt0(cut_data, trial_metrics, 'tss_diam_before','tss_diam_after', 'ts_start', analysis_folder, compare_folder, ...
                    do_plots, conditions, condition_labels, time_vector, x_label_text, color_map, ...
                    'aligned_diams', strcat('ID ', num2str(id), ': Pupil Response Around Target Saccade Starts [zero-aligned]'), ...
                    't0_s', strcat('diam_t0_tss_zero_aligned_subj',num2str(id)), id);
% Fill outliers
diam_t0_tss.cleanedData = filloutliers(diam_t0_tss.aligned_diams,'nearest', 'percentiles',[10 90]); 
if do_plots == 1
    plotCompareTwoStepsDiam(diam_t0_tss, 'aligned_diams', 'ZeroAligned',  ...
                            diam_t0_tss, 'cleanedData',   '80Percentile', ...
                            strcat('ID ', num2str(id), ': Pupil Response Around Target Saccade Starts'), ...
                            color_map, conditions, condition_labels, time_vector, x_label_text, 't0_s',analysis_folder, compare_folder, 'tss',id);
    plotHistNaNsBufferedDiamSacc(diam_t0_tss, analysis_folder, 'Target Saccade NaN Counts before/after t0start', 'diam_t0_tss_nan_counts_histogram');
else
end
save(fullfile(analysis_folder, strcat('\diam_t0_tss.mat')), 'diam_t0_tss');

%%  --- 4.2 plot saccades NOT towards target aligned by START of saccade [raw, 0-aligned, cleaned] ---
diam_t0_ntss = getPupilDiamsBasedOnt0(cut_data, trial_metrics, 'ntss_diam_before','ntss_diam_after', 'nts_start', analysis_folder, compare_folder, ...
                    do_plots, conditions, condition_labels, time_vector, x_label_text, color_map, ...
                    'aligned_diams', strcat('ID ', num2str(id), ': Pupil Response Around non-Target Saccade Starts [ZeroAligned]'), ...
                    't0_s', strcat('diam_t0_ntss_zero_aligned_subj',num2str(id)), id);
%remove outliers
diam_t0_ntss.cleanedData = filloutliers(diam_t0_ntss.aligned_diams,'nearest','percentiles',[10 90]); 
save(fullfile(analysis_folder, strcat('\diam_t0_ntss.mat')), 'diam_t0_ntss');
if do_plots == 1
    plotCompareTwoStepsDiam(diam_t0_ntss, 'aligned_diams', 'ZeroAligned', ...
                            diam_t0_ntss, 'cleanedData',   '80Percentile', ...
                            strcat('ID ', num2str(id), ': Pupil Response Around non-Target Saccade Starts'), ...
                            color_map, conditions, condition_labels, time_vector, x_label_text, 't0_s', analysis_folder, compare_folder, 'ntss',id);
    plotHistNaNsBufferedDiamSacc(diam_t0_ntss, analysis_folder, 'NonTarget Saccade NaN Counts before/after t0start', 'diam_t0_ntss_nan_counts_histogram');
else
end
%%  --- 4.3 plot ALL saccades TOWARDS target aligned by END of saccade [raw, 0-aligned, cleaned] ---
diam_t0_tse = getPupilDiamsBasedOnt0(cut_data, trial_metrics, 'tse_diam_before','tse_diam_after', 'ts_end', analysis_folder, compare_folder, ...
                    do_plots, conditions, condition_labels, time_vector, x_label_text, color_map, ...
                    'aligned_diams', strcat('ID ', num2str(id), ': Pupil Response Around Target Saccade End [ZeroAligned]'), ...
                    't0_e', strcat('diam_t0_tse_zero_aligned_subj',num2str(id)), id);
%remove outliers
diam_t0_tse.cleanedData = filloutliers(diam_t0_tse.aligned_diams,'nearest', 'percentiles',[10 90]); 
save(fullfile(analysis_folder, strcat('\diam_t0_tse.mat')), 'diam_t0_tse');
if do_plots == 1
    plotCompareTwoStepsDiam(diam_t0_tse, 'aligned_diams', 'ZeroAligned', ...
                            diam_t0_tse, 'cleanedData',   '80Percentile', ...
                            strcat('ID ', num2str(id), ': Pupil Response Around Target Saccade Ends'), ...
                            color_map, conditions, condition_labels, time_vector, x_label_text, 't0_e', analysis_folder, compare_folder, 'tse',id);
    plotHistNaNsBufferedDiamSacc(diam_t0_tse, analysis_folder, 'Target Saccade NaN Counts before/after t0end', 'diam_t0_tse_nan_counts_histogram');
else
end
%%  --- 4.4 plot saccades NOT towards target aligned by END of saccade ---
%has different values than NTSS since zero aligning only works if first value is non-NaN
diam_t0_ntse = getPupilDiamsBasedOnt0(cut_data, trial_metrics, 'ntse_diam_before','ntse_diam_after', 'nts_end', analysis_folder, compare_folder, ...
                    do_plots, conditions, condition_labels, time_vector, x_label_text, color_map, ...
                    'aligned_diams', strcat('ID ', num2str(id), ': Pupil Response Around Non-Target Saccade End [ZeroAligned]'),... 
                    't0_e', strcat('diam_t0_ntse_zero_aligned_subj',num2str(id)), id);
%remove outliers                
diam_t0_ntse.cleanedData = filloutliers(diam_t0_ntse.aligned_diams,'nearest', 'percentiles',[10 90]);
save(fullfile(analysis_folder, strcat('\diam_t0_ntse.mat')), 'diam_t0_ntse');
if do_plots == 1
    plotCompareTwoStepsDiam(diam_t0_ntse, 'aligned_diams', 'ZeroAligned', ...
                            diam_t0_ntse, 'cleanedData',   '80Percentile', ...
                            strcat('ID ', num2str(id), ': Pupil Response Around non-Target Saccade Ends'), ...
                            color_map, conditions, condition_labels,time_vector, x_label_text, 't0_e', analysis_folder, compare_folder, 'ntse',id);
    plotHistNaNsBufferedDiamSacc(diam_t0_ntse, analysis_folder, 'NonTarget Saccade NaN Counts before/after t0end', 'diam_t0_ntse_nan_counts_histogram');
else
end
%%  --- 4.5 plot FIRST START target saccade of trial ---
result_table = [];  num_trials = size(cut_data.Condition, 1);
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
%%  --- 4.6 plot FIRST saccade TOWARDS target aligned by END of saccade ---
result_table = [];
num_trials = size(cut_data.Condition, 1);
for ts = 1:num_trials
    current_condition = cut_data.Condition(ts);
    num_target_saccades = height(trial_metrics.tse_diam_before{ts});
    if ~isempty(num_target_saccades)
        try
            diam_before = trial_metrics.tse_diam_before{ts}(1,:);
            diam_after = trial_metrics.tse_diam_after{ts}(1,:);
            diam_combined = [diam_before, diam_after];
            saccade_idx = trial_metrics.ts_end{ts}(1); % t0
            result_table = [result_table; {ts, current_condition, diam_before, diam_after, diam_combined, saccade_idx}];
        catch
            continue; % no target saccades in this trial
        end
    else
    end
end
diam_t0_tse_first = cell2table(result_table, 'VariableNames', {'TrialNumber', 'Condition', 'DataPointsBefore', 'DataPointsAfter', 'DiamCombined', 'SaccadeTrialIndex'});
for ts = 1:size(diam_t0_tse_first, 1) % zero align
    if ~isnan(diam_t0_tse_first.SaccadeTrialIndex(ts))
        try current_diam = diam_t0_tse_first.DiamCombined{ts,:};
        catch, current_diam = diam_t0_tse_first.DiamCombined(ts,:); 
        end
        try diam_t0_tse_first.aligned_diams(ts,:) = current_diam - current_diam(1);  % Aligning each trial to start at zero
        catch, diam_t0_tse_first.aligned_diams(ts,:) = NaN(1,51);
        end
    else
        diam_t0_tse_first.aligned_diams(ts,:) = NaN(1,51);
    end
end

diam_t0_tse_first.cleanedData = filloutliers(diam_t0_tse_first.aligned_diams,'nearest','percentiles',[10 90]);% remove outliers
save(fullfile(analysis_folder, strcat('\diam_t0_tse_first.mat')), 'diam_t0_tse_first');
if do_plots == 1
    plotCompareTwoStepsDiam(diam_t0_tse_first, 'aligned_diams', 'ZeroAligned', ...
                            diam_t0_tse_first, 'cleanedData',   '80Percentile', ...
                            strcat('ID ', num2str(id), ': Pupil Response Around First Target Saccade Ends'), ...
                            color_map, conditions, condition_labels, time_vector, x_label_text, 't0_e', analysis_folder, compare_folder, 'tse_first',id);
    
    plotDiamsPerConditionAndAverage(diam_t0_tse_first, 'aligned_diams', conditions, condition_labels, time_vector, x_label_text, color_map, ...
            strcat('ID ', num2str(id), ': Pupil Response Around First Target Saccade End [ZeroAligned]'), 't0_e',...
            analysis_folder, compare_folder, 'tse_first',id);
    
else
end
%%  --- 4.7 only take start of last target saccades per trial ---
result_table = []; num_trials = size(cut_data.Condition, 1);
for i=1:num_trials
    current_saccades = diam_t0_tss(ismember(diam_t0_tss.TrialNumber, i), :);
    try
        last_ts_index = current_saccades.SaccadeTrialIndex(end);
        last_nts_index = trial_metrics.nts_start{i,:}(end);
        if ~isnan(current_saccades.SaccadeTrialIndex(end)) && last_ts_index > last_nts_index
            result_table = [result_table; {i, current_saccades.Condition(end), current_saccades.SaccadeTrialIndex(end), [current_saccades.aligned_diams(end,:)], [current_saccades.cleanedData(end,:)]}];
        else % if no target saccades in trial skip it
            %result_table = [result_table; {i, current_saccades.Condition(end), current_saccades.SaccadeTrialIndex(end), [NaN(1,51)], [NaN(1,51)]}];
        end
    catch % if theres no target saccades in this trial ignore
        %result_table = [result_table; {i, current_saccades.Condition(end), NaN, [NaN(1,51)], [NaN(1,51)]}];
    end
end
diam_t0_tss_last = cell2table(result_table, 'VariableNames', {'TrialNumber', 'Condition', 'SaccadeTrialIndex', 'aligned_diams', 'cleanedData'});
diam_t0_tss_last.cleanedData = filloutliers(diam_t0_tss_last.aligned_diams,'nearest','percentiles',[10 90]);% remove outliers
if do_plots == 1
    plotDiamsPerConditionAndAverage(diam_t0_tss_last, 'cleanedData', conditions, condition_labels, time_vector, x_label_text, color_map, ...
            strcat('ID ', num2str(id), ': Pupil Response if Last Saccade is to Target [80Percentile]'), 't0_s',...
            analysis_folder, compare_folder, 'tss_last_cleaned',id);
    plotDiamsPerConditionAndAverage(diam_t0_tss_last, 'aligned_diams', conditions, condition_labels, time_vector, x_label_text, color_map, ...
            strcat('ID ', num2str(id), ': Pupil Response if Last Saccade is to Target [ZeroAligned]'), 't0_s', ...
            analysis_folder, compare_folder, 'tss_last',id);
else
end
%%  --- 4.8 combo of zero aligned data means per condition ---
y_axis_labels = [];
plotPupilDiameterResponses4subplotsAveragePerCondition( ...
        diam_t0_tss,  'aligned_diams', 'START of Target Saccades', 't0_s',...
        diam_t0_ntss, 'aligned_diams', 'START of non-TargetSaccades', 't0_s',...
        diam_t0_tse,  'aligned_diams', 'END of Target Saccades', 't0_e',...
        diam_t0_ntse, 'aligned_diams', 'END of non-TargetSaccades', 't0_e',...
        strcat('ID ', num2str(id), ': Average Pupil Response per Condition [ZeroAligned]'), ...
        conditions, condition_labels, time_vector, color_map,...
        0, [min(time_vector), max(time_vector)], y_axis_labels, [0 0 1.2 0.7]);
print(gcf, fullfile(analysis_folder,strcat('diam_t0_averages_tss_ntss_tse_ntse.svg')), '-dsvg');
print(gcf, fullfile(compare_folder, strcat('diam_t0_averages_tss_ntss_tse_ntse_id',num2str(id),'.svg')), '-dsvg');
print(gcf, fullfile(compare_folder, strcat('diam_t0_averages_tss_ntss_tse_ntse_id',num2str(id),'.svg')), '-dsvg');

%%  --- 4.9 TSE TSS FIRST START: combo plot of aligned & cleaned data means per condition ---
plotPupilDiameterResponses4subplotsAveragePerCondition( ...
        diam_t0_tss_first,  'aligned_diams', 'START of First Target Saccades', 't0_s', ...
        diam_t0_ntss,       'aligned_diams', 'START of non-TargetSaccades',    't0_s', ...
        diam_t0_tse_first,  'aligned_diams', 'END of First Target Saccades',   't0_e', ...
        diam_t0_ntse,       'aligned_diams', 'END of non-TargetSaccades',      't0_e',  ...
        strcat('ID ', num2str(id), ': Average Pupil Response per Condition [ZeroAligned]'), ...
        conditions, condition_labels, time_vector, color_map, ...
        0, [min(time_vector), max(time_vector)],y_axis_labels, [0 0 1.2 0.7]);
print(gcf, fullfile(analysis_folder, strcat('diam_t0_average_tssFIRST_ntss_tseFIRST_ntse.svg')), '-dsvg');
print(gcf, fullfile(compare_folder,  strcat('diam_t0_average_tssFIRST_ntss_tseFIRST_ntse_id',num2str(id),'.svg')), '-dsvg');
saveas(gcf, fullfile(compare_folder, strcat('diam_t0_average_tssFIRST_ntss_tseFIRST_ntse_id',num2str(id),'.png')));

plotPupilDiameterResponses4subplotsAveragePerCondition( ...
        diam_t0_tss_first,  'cleanedData', 'START of First Target Saccades', 't0_s', ...
        diam_t0_ntss,       'cleanedData', 'START of non-TargetSaccades',    't0_s', ...
        diam_t0_tse_first,  'cleanedData', 'END of First Target Saccades',   't0_e', ...
        diam_t0_ntse,       'cleanedData', 'END of non-TargetSaccades',      't0_e', ...
        strcat('ID ', num2str(id), ': Average Pupil Response per Condition [80Percentile]'), ...
        conditions, condition_labels, time_vector, color_map, ...
        0, [min(time_vector), max(time_vector)], y_axis_labels, [0 0 1.2 0.7]);
print(gcf, fullfile(analysis_folder, strcat('diam_t0_tssFIRST_ntss_tseFIRST_ntse_cleaned.svg')), '-dsvg');
print(gcf, fullfile(compare_folder,  strcat('diam_t0_tssFIRST_ntss_tseFIRST_ntse_cleaned_id',num2str(id),'.svg')), '-dsvg');
