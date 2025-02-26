function [diam_t0, diam_t0_avg] = calcDiamNormalizePlot(diam_t0, baseline, outlier_tolerance, baseline_min, num_before, num_after, conditions, base, id, analysis_folder, color_map, sr)

time_vector = (-num_before:num_after - 1) * sr; % Time in seconds
x_label_text = 'Time from t0 (s)';

%% 0. plot the raw cut data
num_rows = 2;
num_col = 2;
figure; sgtitle(strcat('Subject ', num2str(id), ': Pupil Response Before and After Gaze Reached Target [Raw]'));
subplot(num_rows, num_col, 1); hold on; title('Individual Trials (Raw)');
for i = 1:4
    condition = conditions{i};
    condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
    for j = 1:height(condition_trials)
        current_diam = [condition_trials.DataPointsBefore(j, :), condition_trials.DataPointsAfter(j, :)];
        plot(time_vector, current_diam, 'Color', color_map(condition));
    end
end
xline(0, '--', 'Gaze at Target', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
xlabel(x_label_text); ylabel('Pupil Diameter (mm)'); hold off;

subplot(num_rows, num_col, 2); hold on; title('Average per Condition (Raw)');
for i = 1:4
    condition = conditions{i};
    condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
    current_diams = [condition_trials.DataPointsBefore, condition_trials.DataPointsAfter];
    plotConditionDataWithShading(current_diams, time_vector, color_map(condition), condition);
end
xline(0, '--', 'Gaze at Target', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
xlabel(x_label_text); ylabel('Pupil Diameter (mm)'); legend('show'); legend('boxoff');hold off;

% set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
% saveas(gcf, fullfile(analysis_folder, 'diam_t0_0_raw.png'));


%%  Second plot: Zero-aligned pupil diameter (pi(t) - pi(t1))
% Concatenate DataPointsBefore and DataPointsAfter for each trial
% Calculate pi(t) - pi(t1) by subtracting the first value of DataPointsBefore
for trial = 1:size(diam_t0, 1)
    current_diam = [diam_t0.DataPointsBefore(trial, :), diam_t0.DataPointsAfter(trial, :)];
    diam_t0.aligned_diams(trial,:) = current_diam - current_diam(1);  % Aligning each trial to start at zero
end
        subplot(num_rows, num_col, 3); hold on; title('Individual Trials (Zero-Aligned)');
        for i = 1:4
            condition = conditions{i};
            condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
            for j = 1:height(condition_trials)
                plot(time_vector, condition_trials.aligned_diams(j,:), 'Color', color_map(condition));
            end
        end
        xline(0, '--', 'Gaze at Target', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
        xlabel(x_label_text); ylabel('Pupil Diameter Change (mm)'); hold off;

        subplot(num_rows, num_col, 4); hold on; title('Average per Condition (Zero-Aligned)');
        for i = 1:4
            condition = conditions{i};
            condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
            plotConditionDataWithShading(condition_trials.aligned_diams, time_vector, color_map(condition), condition);
        end
        xline(0, '--', 'Gaze at Target', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
        xlabel(x_label_text); ylabel('Pupil Diameter (mm)'); legend('show'); legend('boxoff');hold off;

        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf, fullfile(analysis_folder, 'diam_t0_0_raw.png'));





% Set figure position and size
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

% Save the figure if needed
% saveas(gcf, fullfile(analysis_folder, 'diam_t0_0_raw_zero_aligned.png'));






%% 1. Normalize the diameter data relative to the baseline first [which baseline?] directly before or the mean one?
for trial = 1:size(diam_t0, 1)
    current_diam = [diam_t0.DataPointsBefore(trial, :), diam_t0.DataPointsAfter(trial, :)];
    diam_t0.t0_diam_mean_baseline(trial,:)    = current_diam - mean(baseline);
    diam_t0.t0_diam_current_baseline(trial,:) = current_diam - mean(diam_t0.DiamBothStart(trial, :));
end
        num_rows = 2; num_col = 2; figure; sgtitle(strcat('Subject ', num2str(id), ':Pupil Response Before and After Gaze Reached Target (', base, ')'));
        % t0_diam_mean_baseline
        subplot(num_rows, num_col, 1); hold on; title('Individual Trials [Mean Baseline]');
        for i = 1:4
            condition = conditions{i}; condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
            for j = 1:height(condition_trials)
                plot(time_vector, condition_trials.t0_diam_mean_baseline(j,:), 'Color', color_map(condition));
            end
        end
        xline(0, '--', 'Gaze at Target', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
        xlabel(x_label_text); ylabel('Pupil Diameter (mm)'); hold off;

        subplot(num_rows, num_col, 2); hold on; title('Average per Condition [Mean Baseline]');
        for i = 1:4
            condition = conditions{i};
            condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
            plotConditionDataWithShading(condition_trials.t0_diam_mean_baseline, time_vector, color_map(condition), condition);
        end
        xline(0, '--', 'Gaze at Target', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
        xlabel(x_label_text); ylabel('Pupil Diameter (mm)'); legend('show'); legend('boxoff'); hold off;

% t0_diam_current_baseline
subplot(num_rows, num_col, 3); hold on; title('Individual Trials [Mean(CurrentBaseline)]');
for i = 1:4
    condition = conditions{i};
    condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
    for j = 1:height(condition_trials)
        plot(time_vector, condition_trials.t0_diam_current_baseline(j,:), 'Color', color_map(condition));
    end
end
xline(0, '--', 'Gaze at Target', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
xlabel(x_label_text); ylabel('Pupil Diameter (mm)'); hold off;

subplot(num_rows, num_col, 4); hold on; title('Average per Condition [Mean(CurrentBaseline)]');
for i = 1:4
    condition = conditions{i};
    condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
    plotConditionDataWithShading(condition_trials.t0_diam_current_baseline, time_vector, color_map(condition), condition);
end
xline(0, '--', 'Gaze at Target', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
xlabel(x_label_text); ylabel('Pupil Diameter (mm)'); legend('show'); legend('boxoff'); hold off;

set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf, fullfile(analysis_folder, strcat('diam_t0_', base, '_1_baseline.png')));

%% 2. Interpolate only internal NaNs
for trial = 1:size(diam_t0.t0_diam_mean_baseline, 1)
    diam_int = interpolateInternalNaNs(diam_t0.t0_diam_mean_baseline(trial,:));
    diam_t0.diamInt(trial,:) = diam_int;
end
num_rows = 2; num_col = 2; figure; sgtitle(strcat('Subject ', num2str(id), ': Pupil Response Before and After Gaze Reached Target [Interpolate only internal NaNs & Remove Trials with too many NaNs]'));

subplot(num_rows, num_col, 1); hold on; title('Individual Trials [Interpolate NaNs]');
for i = 1:4
    condition = conditions{i}; condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
    for j = 1:height(condition_trials)
        plot(time_vector, condition_trials.diamInt(j,:), 'Color', color_map(condition));
    end
end
xline(0, '--', 'Gaze at Target', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
xlabel(x_label_text); ylabel('Pupil Diameter (mm)'); hold off;

subplot(num_rows, num_col, 2); hold on; title('Average per Condition [Interpolate NaNs]');
for i = 1:4
    condition = conditions{i};
    condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
    plotConditionDataWithShading(condition_trials.diamInt, time_vector, color_map(condition), condition);
end
xline(0, '--', 'Gaze at Target', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
xlabel(x_label_text); ylabel('Pupil Diameter (mm)'); legend('show'); legend('boxoff'); hold off;


%% 3. Check if list has enough valid data points (e.g., more than 50%)
for trial = 1:size(diam_t0.diamInt, 1)
    diam_int = diam_t0.diamInt(trial,:);
    if numel(diam_int(~isnan(diam_int))) >= abs((size(diam_int,2))/2)
        diam_t0.diamIntValid(trial,:) =  diam_int;
    else % if too many nans diregard trial
        diam_t0.diamIntValid(trial,:) = NaN;
    end
end
subplot(num_rows, num_col, 3); hold on; title('Individual Trials [>12 datapoints]');
for i = 1:4
    condition = conditions{i}; condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
    for j = 1:height(condition_trials)
        plot(time_vector, condition_trials.diamIntValid(j,:), 'Color', color_map(condition));
    end
end
xline(0, '--', 'Gaze at Target', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
xlabel(x_label_text); ylabel('Pupil Diameter (mm)');hold off;

subplot(num_rows, num_col, 4); hold on; title('Average per Condition [>12 datapoints]');
for i = 1:4
    condition = conditions{i};
    condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
    plotConditionDataWithShading(condition_trials.diamIntValid, time_vector, color_map(condition), condition);
end
xline(0, '--', 'Gaze at Target', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
xlabel(x_label_text); ylabel('Pupil Diameter (mm)'); legend('show'); legend('boxoff'); hold off;

set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf, fullfile(analysis_folder, strcat('diam_t0_', base,'_23_interpolate_50valid.png')));





%% Step 4.1: Outlier Detection Across All Trials and timepoints
pupil_derivative = diff(diam_t0.diamIntValid, 1, 1);

mean_derivative = mean(pupil_derivative(:), 'omitnan');
std_derivative = std(pupil_derivative(:), 'omitnan');
outliers = abs(pupil_derivative) > (mean_derivative + outlier_tolerance * std_derivative);
outlier_mask = [false(1, size(outliers, 2)); outliers];
diam_t0.diamOutliersRemoved = diam_t0.diamIntValid;
diam_t0.diamOutliersRemoved(outlier_mask) = NaN;

% Step 5.1: Interpolate Data After Outlier Removal Across All Trials/timepoints
diam_t0.diamIntOutliers = diam_t0.diamOutliersRemoved;
for trial = 1:size(diam_t0.diamOutliersRemoved, 1)
    diam_t0.diamIntOutliers(trial, :) = interpolateInternalNaNs(diam_t0.diamOutliersRemoved(trial, :));
end

%% Step 4.2: Outlier Detection Per Time Point
mean_derivative_per_time = mean(pupil_derivative, 1, 'omitnan');
std_derivative_per_time = std(pupil_derivative, 0, 1, 'omitnan');
outliers_per_time = abs(pupil_derivative) > (mean_derivative_per_time + outlier_tolerance * std_derivative_per_time);
outlier_mask_per_time = [false(1, size(outliers_per_time, 2)); outliers_per_time];
diam_t0.diamOutliersRemovedPerTime = diam_t0.diamIntValid;
diam_t0.diamOutliersRemovedPerTime(outlier_mask_per_time) = NaN;

% Step 5.2: Interpolate Data After Outlier Removal Per Time Point
diam_t0.diamIntOutliersPerTime = diam_t0.diamOutliersRemovedPerTime;
for trial = 1:size(diam_t0.diamOutliersRemovedPerTime, 1)
    diam_t0.diamIntOutliersPerTime(trial, :) = interpolateInternalNaNs(diam_t0.diamOutliersRemovedPerTime(trial, :));
end

% plot all 4
plotOutlierAndInterpolationComparison(diam_t0, base, conditions, time_vector, color_map, x_label_text, outlier_tolerance, id, analysis_folder)




%% Step 6: Z-score normalization

%% Global Time-Point Mean Z-Scoring: centers each time point around the mean across all trials 
%% for that time, emphasizing condition-wide changes that are time-locked.

% z-score normalization per trial - cant be time locked-> mean and std per trial

% z-score normalization per condition
for i = 1:4
    condition = conditions{i};
    condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
    condition_idx = find(strcmp(diam_t0.Condition, condition));
    % Calculate the mean and standard deviation, ignoring NaNs
    mean_trial = mean(condition_trials.diamIntValid,1,'omitnan');
    std_trial = std(condition_trials.diamIntValid, 1,'omitnan');
    for trial = 1:size(condition_trials,1)
        trial_data = condition_trials.diamIntValid(trial, :);
        % Perform z-score normalization on non-NaN values
        trial_index_in_diam_t0 = condition_idx(trial);
        diam_t0.diamZtimeCondition(trial_index_in_diam_t0,:) = (trial_data - mean_trial) ./ std_trial;
    end
end

%% % z-score all trials
% With group-level z-scoring, all trials share the same baseline, 
% which can help in comparing values across trials directly.
group_mean = mean(diam_t0.diamIntValid, 1, 'omitnan'); % Calculate mean across all trials, ignoring NaNs
group_std = std(diam_t0.diamIntValid, 1, 'omitnan');   % Calculate std across all trials, ignoring NaNs
% Z-score each trial based on the group-level mean and std
for trial = 1:size(diam_t0.diamIntValid, 1)
    trial_data = diam_t0.diamIntValid(trial, :);
    % Perform z-score normalization on non-NaN values
    diam_t0.diamZtimeGroup(trial, :) = (trial_data - group_mean) ./ group_std;
end



%% Global Mean Z-Scoring: This single mean and standard deviation for the whole dataset treats every 
%% data point as comparable, showing overall deviations from the grand baseline.

% z-score normalization per trial
% In trial-level z-scoring, each trial has a different mean and standard deviation, 
% so the z-scores are relative to that trial.
for trial = 1:size(diam_t0.diamIntValid, 1)
    % Extract the data for the current trial
    trial_data = diam_t0.diamIntValid(trial, :);
    % Calculate the mean and standard deviation, ignoring NaNs
    mean_trial = mean(trial_data, 'omitnan');
    std_trial = std(trial_data, 'omitnan');
    % Perform z-score normalization on non-NaN values
    diam_t0.diamZ(trial,:) = (trial_data - mean_trial) / std_trial;
end

%% z-score normalization per condition
for i = 1:4
    condition = conditions{i};
    condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
    condition_idx = find(strcmp(diam_t0.Condition, condition));
    % Calculate the mean and standard deviation, ignoring NaNs
    mean_trial = mean(condition_trials.diamIntValid(:),1,'omitnan');
    std_trial = std(condition_trials.diamIntValid(:),1,'omitnan');
    for trial = 1:size(condition_trials,1)
        trial_data = condition_trials.diamIntValid(trial, :);
        % Perform z-score normalization on non-NaN values
        trial_index_in_diam_t0 = condition_idx(trial);
        diam_t0.diamZcondition(trial_index_in_diam_t0,:) = (trial_data - mean_trial) / std_trial;
    end
end

% z-score all trials
% With group-level z-scoring, all trials share the same baseline, 
% which can help in comparing values across trials directly.
all_data = diam_t0.diamIntValid(:); % Flatten to a single vector
group_mean = mean(all_data, 'omitnan'); % Calculate mean across all trials, ignoring NaNs
group_std = std(all_data, 'omitnan');   % Calculate std across all trials, ignoring NaNs
% Z-score each trial based on the group-level mean and std
for trial = 1:size(diam_t0.diamIntValid, 1)
    trial_data = diam_t0.diamIntValid(trial, :);
    % Perform z-score normalization on non-NaN values
    diam_t0.diamZgroup(trial, :) = (trial_data - group_mean) / group_std;
end


%% each z-scoring method
% nontime locked
plozZscoring(diam_t0, 'diamZ', conditions, time_vector, color_map, 'Time (s)', ...
    strcat('Subject ', num2str(id), ': Pupil Response Before and After Gaze Reached Target [z-scored per trial]'));
saveas(gcf, fullfile(analysis_folder, strcat('diam_t0_', base, '_6_zscore_per_trial.png')));

plozZscoring(diam_t0, 'diamZgroup', conditions, time_vector, color_map, 'Time (s)', ...
    strcat('Subject ', num2str(id), ': Pupil Response Before and After Gaze Reached Target [z-scored over all trials]'));
saveas(gcf, fullfile(analysis_folder, strcat('diam_t0_', base, '_6_zscore_all.png')));

plozZscoring(diam_t0, 'diamZcondition', conditions, time_vector, color_map, 'Time (s)', ...
    strcat('Subject ', num2str(id), ': Pupil Response Before and After Gaze Reached Target [z-scored per condition]'));
saveas(gcf, fullfile(analysis_folder, strcat('diam_t0_', base, '_6_zscore_per_condition.png')));

% time locked
plozZscoring(diam_t0, 'diamZtimeGroup', conditions, time_vector, color_map, 'Time (s)', ...
    strcat('Subject ', num2str(id), ': Pupil Response Before and After Gaze Reached Target [time-locked z-scored over all trials]'));
saveas(gcf, fullfile(analysis_folder, strcat('diam_t0_', base, '_6_zscore_timelocked_all.png')));

plozZscoring(diam_t0, 'diamZtimeCondition', conditions, time_vector, color_map, 'Time (s)', ...
    strcat('Subject ', num2str(id), ': Pupil Response Before and After Gaze Reached Target [time-locked z-scored per condition]'));
saveas(gcf, fullfile(analysis_folder, strcat('diam_t0_', base, '_6_zscore_timelocked_per_condition.png')));




%%
%plotAllZScoringMethods(diam_t0.diamIntValid, conditions, time_vector, color_map, 'Time (s)');


save(fullfile(analysis_folder, strcat('\diam_t0_normalized_', base,'.mat')), 'diam_t0');










%% part two min max mean
diam_t0_avg = struct();
% for c = 1:length(conditions)
%     condition = conditions{c};
%     condition_mask = strcmp(diam_t0.Condition, condition);
%     condition_diam = diam_t0.diamZgroup(condition_mask,:);
%     condition_idx = find(strcmp(diam_t0.Condition, condition));
%     
%     for trial = 1:size(condition_diam, 1)
%         current_diam = condition_diam(trial,:);
%         % Compute the maximum pupil dilation after shortly before and after target was found onset
%         trial_idx_diam_t0 = condition_idx(trial);
%         diam_t0_avg.(condition).maxDiam(trial_idx_diam_t0)  = max(current_diam);
%         diam_t0_avg.(condition).minDiam(trial_idx_diam_t0)  = min(current_diam);
%         diam_t0_avg.(condition).meanDiam(trial_idx_diam_t0) = mean(current_diam,'omitnan');
%     end
%     diam_t0_avg.(condition).maxMean  = max(diam_t0_avg.(condition).maxDiam);
%     diam_t0_avg.(condition).minMean  = min(diam_t0_avg.(condition).minDiam);
%     diam_t0_avg.(condition).meanMean = mean(diam_t0_avg.(condition).meanDiam, 'omitnan');
% end
% save(fullfile(analysis_folder, strcat(strcat('diam_t0_', base, '_normalized_maxmeanmin.mat')), 'diam_t0_avg'));


% Plot the normalized pupil responses before and after the target is found
%plotPupilResponses(condition_diam_avg, id, conditions, analysis_folder, color_map, sr, base)
end


%%
function plotConditionDataWithShading(data, x_values, color, condition)
% Plot mean ± SEM with shading for each condition
avg = mean(data, 1, 'omitnan');
std_error = std(data, 0, 1, 'omitnan') / sqrt(size(data, 1));
nan_mask = ~isnan(avg) & ~isnan(std_error);
x_shaded = x_values(nan_mask);
y_upper = avg(nan_mask) + std_error(nan_mask);
y_lower = avg(nan_mask) - std_error(nan_mask);
fill([x_shaded, fliplr(x_shaded)], [y_upper, fliplr(y_lower)], color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(x_values, avg, 'Color', color, 'LineWidth', 2, 'DisplayName', condition);
end


%%
% Helper function to interpolate only internal NaNs (between two valid points)
function interpolated_data = interpolateInternalNaNs(pupil_data)
% Find the first and last non-NaN values
first_valid_idx = find(~isnan(pupil_data), 1, 'first');
last_valid_idx = find(~isnan(pupil_data), 1, 'last');

% Only interpolate NaNs that are between two valid data points
interpolated_data = pupil_data;  % Start with the original data

% Use interpolation only on the range between first_valid_idx and last_valid_idx
interpolated_range = first_valid_idx:last_valid_idx;

% Interpolate internal NaNs within the valid range
interpolated_data(interpolated_range) = fillmissing(pupil_data(interpolated_range), 'spline');
end


%%
function plotOutlierAndInterpolationComparison(diam_t0, base, conditions, time_vector, color_map, x_label_text, outlier_tolerance, id, analysis_folder)
    % diam_t0: Table containing pupil data
    % conditions: Cell array of condition names
    % time_vector: Time axis for x-axis labeling
    % color_map: Map of colors for each condition
    % x_label_text: Label for x-axis
    % outlier_tolerance: Number of standard deviations to consider an outlier
    % id: Subject ID
    % analysis_folder: Folder to save the figures


    % Create the figure
    figure;
    sgtitle(['Subject ', num2str(id), ': Pupil Response Outlier Detection & Interpolation Comparison']);

    % Titles for each subplot pair
    method_titles = {'Outliers Removed (Overall)', 'Outliers Removed (Per Time)', ...
                     'Outliers Interpolated (Overall)', 'Outliers Interpolated (Per Time)'};

    % Data to plot for each method
    plot_data = {diam_t0.diamOutliersRemoved, diam_t0.diamOutliersRemovedPerTime, ...
                 diam_t0.diamIntOutliers, diam_t0.diamIntOutliersPerTime};

    % Plot each method with individual trials and condition averages
    for method = 1:4
        % Plot individual trials
        subplot(4, 2, (method - 1) * 2 + 1); hold on; title([method_titles{method}, ' - Individual Trials']);
        for i = 1:size(conditions,2)
            condition = conditions{i};
            condition_trials = plot_data{method}(strcmp(diam_t0.Condition, condition), :);
            for j = 1:size(condition_trials, 1)
                plot(time_vector, condition_trials(j, :), 'Color', color_map(condition));
            end
        end
        xline(0, '--', 'Gaze at Target', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3]);
        xlabel(x_label_text); ylabel('Pupil Diameter (mm)'); hold off;

        % Plot condition averages with shading
        subplot(4, 2, method * 2); hold on; title([method_titles{method}, ' - Average per Condition']);
        for i = 1:size(conditions,2)
            condition = conditions{i};
            condition_trials = plot_data{method}(strcmp(diam_t0.Condition, condition), :);
            plotConditionDataWithShading(condition_trials, time_vector, color_map(condition), condition)
        end
        xline(0, '--', 'Gaze at Target', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3]);
        xlabel(x_label_text); ylabel('Pupil Diameter (mm)'); hold off;
    end
    legend('show'); legend('boxoff'); 
    % Full figure adjustments
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf, fullfile(analysis_folder, strcat('diam_t0_', base, '_outlier_interpolation_comparison.png')));
end



%% 
function plozZscoring(diam_data, z_method, conditions, time_vector, color_map, x_label_text, title_text)
% diam_data: Matrix of z-scored pupil data to plot
% conditions: Cell array of condition names
% time_vector: Time axis for x-axis labeling
% color_map: Map of colors for each condition
% x_label_text: Label for x-axis
% title_text: Title for the plot
num_conditions = length(conditions);

figure;
sgtitle(title_text); % Set overall title for the figure

% Plot individual trials for each condition
subplot(1, 2, 1); hold on; title('Individual Trials]');
    for i = 1:num_conditions
        condition = conditions{i};
        condition_trials = diam_data(strcmp(diam_data.Condition, condition), :);
        for j = 1:height(condition_trials)
            plot(time_vector, condition_trials.(z_method)(j, :), 'Color', color_map(condition));
        end
    end
    xline(0, '--', 'Gaze at Target', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3]);
    xlabel(x_label_text); ylabel('Pupil Diameter (z-scored)'); hold off;
   
    % Plot condition averages with shading for each condition
    subplot(1, 2, 2); hold on; title('Average per Condition');
    for i = 1:num_conditions
        condition = conditions{i};
        condition_trials = diam_data(strcmp(diam_data.Condition, condition), :);
        plotConditionDataWithShading(condition_trials.(z_method), time_vector, color_map(condition), condition);

        % Calculate mean and standard error across trials in the condition
%         mean_vals = nanmean(condition_trials.(z_method), 1);
%         std_err = nanstd(condition_trials.(z_method), 0, 1) / sqrt(size(condition_trials.(z_method), 1));
%         % Plot mean with shading for standard error
%         fill([time_vector, fliplr(time_vector)], ...
%             [mean_vals + std_err, fliplr(mean_vals - std_err)], ...
%             color_map(condition), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%         plot(time_vector, mean_vals, 'Color', color_map(condition), 'LineWidth', 2);
    end
     xline(0, '--', 'Gaze at Target', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3]);
     xlabel(x_label_text); ylabel('Pupil Diameter (z-scored)'); hold off;

legend(conditions, 'Location', 'bestoutside'); legend('boxoff');
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Fullscreen figure
end




%% 
function plotAllZScoringMethods(diam_data, conditions, time_vector, color_map, x_label_text)
% diam_data: Original pupil data (pre-interpolated and outlier-removed)
% conditions: Cell array of condition names
% time_vector: Time axis for x-axis labeling
% color_map: Map of colors for each condition
% x_label_text: Label for x-axis

num_conditions = length(conditions);
num_z_methods = 4; % Number of z-scoring methods to compare

% Pre-allocate for different z-scoring approaches
diamZTrial = NaN(size(diam_data));
diamZCondition = NaN(size(diam_data));
diamZTimePoint = NaN(size(diam_data));
diamZGlobal = NaN(size(diam_data));

% 1. Trial-Level Z-Scoring (Per Trial)
for trial = 1:size(diam_data, 2)
    trial_data = diam_data(:, trial);
    mean_trial = mean(trial_data, 'omitnan');
    std_trial = std(trial_data, 'omitnan');
    diamZTrial(:, trial) = (trial_data - mean_trial) / std_trial;
end

% 2. Condition-Level Z-Scoring (Per Condition)
for i = 1:num_conditions
    condition = conditions{i};
    condition_trials = diam_data(strcmp(diam_data.Condition, condition), :);
    mean_condition = nanmean(condition_trials, 1);
    std_condition = nanstd(condition_trials, 0, 1);
    diamZCondition(strcmp(diam_data.Condition, condition), :) = ...
        (condition_trials - mean_condition) ./ std_condition;
end

% 3. Time-Point Mean Z-Scoring (Global Time-Point Mean)
time_point_mean = nanmean(diam_data, 2);
time_point_std = nanstd(diam_data, 0, 2);
diamZTimePoint = (diam_data - time_point_mean) ./ time_point_std;

% 4. Global Mean Z-Scoring (Across All Trials and Time Points)
global_mean = nanmean(diam_data(:));
global_std = nanstd(diam_data(:));
diamZGlobal = (diam_data - global_mean) / global_std;

% Create a figure with subplots for each z-scoring method and condition
figure;
sgtitle('Comparison of Z-Scoring Methods for Pupil Response');

% Z-scoring method titles
method_titles = {'Trial-Level Z-Scoring', 'Condition-Level Z-Scoring', ...
    'Time-Point Z-Scoring', 'Global Z-Scoring'};

% Plotting each z-scoring method with a row per method
for z_method = 1:num_z_methods
    % Select the appropriate data for the current z-scoring method
    switch z_method
        case 1
            z_data = diamZTrial;
        case 2
            z_data = diamZCondition;
        case 3
            z_data = diamZTimePoint;
        case 4
            z_data = diamZGlobal;
    end
    
    % Loop over each condition for the current z-scoring method
    for i = 1:num_conditions
        condition = conditions{i};
        condition_trials = z_data(strcmp(diam_data.Condition, condition), :);
        %condition_idx = find(strcmp(diam_data.Condition, condition));

        % Calculate mean and standard error for shaded plot
        mean_vals = nanmean(condition_trials, 1);
        std_err = nanstd(condition_trials, 0, 1) / sqrt(size(condition_trials, 1));
        
        % Plot average with shading in the appropriate subplot
        subplot(num_z_methods, num_conditions, (z_method - 1) * num_conditions + i); hold on;
        title([method_titles{z_method}, ' - ', condition]);
        
        fill([time_vector, fliplr(time_vector)], ...
            [mean_vals + std_err, fliplr(mean_vals - std_err)], ...
            color_map(condition), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(time_vector, mean_vals, 'Color', color_map(condition), 'LineWidth', 2);
        
        xline(0, '--', 'Gaze at Target', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3]);
        xlabel(x_label_text); ylabel('Pupil Diameter (z-scored)');
        hold off;
    end
end

% Set full figure size and layout
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Fullscreen figure
legend(conditions, 'Location', 'bestoutside'); legend('boxoff');
end
