function diam_t0 = getPupilDiamsBasedOnt0(cut_data, trial_metrics, which_diam_before,which_diam_after, s_start, analysis_folder, compare_folder, ...
                                          do_plots, conditions, condition_labels, time_vector, x_label_text, color_map, ... 
                                          alignement, plot_title, t0, safe_name, id)

result_table = [];
num_trials = size(cut_data.Condition, 1);
for ts = 1:num_trials
    current_condition = cut_data.Condition(ts);
    num_target_saccades = height(trial_metrics.(which_diam_before){ts});
    for s = 1:num_target_saccades
        diam_before = trial_metrics.(which_diam_before){ts}(s,:);
        diam_after = trial_metrics.(which_diam_after){ts}(s,:);
        diam_combined = [diam_before, diam_after];
        saccade_idx = trial_metrics.(s_start){ts}(s,:);
        result_table = [result_table; {ts, current_condition, diam_before, diam_after, diam_combined, saccade_idx}];
    end
end
diam_t0 = cell2table(result_table, 'VariableNames', {'TrialNumber', 'Condition', 'DataPointsBefore', 'DataPointsAfter', 'DiamCombined', 'SaccadeTrialIndex'});

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

% plot
if do_plots == 1
    plotDiamsPerConditionAndAverage(diam_t0, 'aligned_diams', conditions, condition_labels, ...
                                    time_vector, x_label_text, color_map, plot_title, ...
                                    t0, analysis_folder, compare_folder, t0, id);
                                                                
    saveas(gcf, fullfile(analysis_folder, strcat(safe_name,'.png')));
    print(gcf, fullfile(compare_folder, strcat(safe_name,'.svg')), '-dsvg');
else
end