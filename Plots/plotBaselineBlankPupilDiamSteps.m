function diam_t0 = plotBaselineBlankPupilDiamSteps(diam_t0, baseline_blank_diam_mean, outlier_tolerance, baseline_min, conditions, color_map, analysis_folder)

%% individual baselines
num_rows = 4;
num_col = 2;
figure
sgtitle('Pupil Response During Baseline (Blank Screen)')
subplot(num_rows, num_col, 1); hold on;
%y_limits = [];
for i = 1:4
    condition = conditions{i};
    condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
    for j = 1:height(condition_trials)
        plot(condition_trials.DiamBothBlank(j,:), 'Color', color_map(condition)); % Plot each trial in light grey        y_limits = [y_limits; ylim];  % Store y-axis limits
    end % y_limits = [y_limits; ylim];  % Store y-axis limits
end
title('Individual Trials (Raw)');
xlabel('Datapoints During BlankScreen'); ylabel('Pupil Diameter (mm)');
hold off;

subplot(num_rows, num_col, 2); hold on;
y_limit_avg = []; 
for i = 1:4
    condition = conditions{i};
    condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
    avg_response = nanmean(condition_trials.DiamBothBlank, 1);
    plot(avg_response, 'Color', color_map(condition), 'DisplayName', condition, 'LineWidth', 1.5);
end
hold off;
title('Average per Condition (Raw)');
xlabel('Datapoints During BlankScreen'); ylabel('Pupil Diameter (mm)'); 
legend('show');
legend('boxoff');
y_limit_avg = [y_limit_avg; ylim];


%% individual baselines - baseline corrected
subplot(num_rows, num_col, 3); hold on;
diam_t0.DiamBothBlankBaselineAdjusted = diam_t0.DiamBothBlank - baseline_blank_diam_mean;
for i = 1:4
    condition = conditions{i};
    condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
    for j = 1:height(condition_trials)
        plot(condition_trials.DiamBothBlankBaselineAdjusted(j,:), 'Color', color_map(condition)); % Plot each trial in light grey        y_limits = [y_limits; ylim];  % Store y-axis limits
    end
end
title('Individual Trials (-avg over all baselines)');
xlabel('Datapoints During BlankScreen'); ylabel('Pupil Diameter (mm)');
hold off;

subplot(num_rows, num_col, 4); hold on;
for i = 1:4
    condition = conditions{i};
    condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
    avg_response = nanmean(condition_trials.DiamBothBlankBaselineAdjusted, 1);
    plot(avg_response, 'Color', color_map(condition), 'DisplayName', condition, 'LineWidth', 1.5);
end
hold off;
title('Average Pupil Response per Condition (-avg over all baselines)');
xlabel('Datapoints During BlankScreen'); ylabel('Pupil Diameter (mm)');


%% individual baselines - remove outliers
subplot(num_rows, num_col, 5); hold on;
diam_t0.DiamBothBlankOutliersRemoved = diam_t0.DiamBothBlankBaselineAdjusted;
% Calculate the derivative across all data at once
pupil_derivative = diff(diam_t0.DiamBothBlankBaselineAdjusted, 1,1);
mean_derivative  = mean(pupil_derivative(:),'omitnan');
std_derivative   = std(pupil_derivative(:),'omitnan');
% Identify outliers in the derivative (across all trials)
outliers = abs(pupil_derivative) > (mean_derivative + outlier_tolerance * std_derivative);
% Create a mask for the outlier values in `pupil_data`
outlier_mask = [false(1, size(outliers, 2)); outliers]; % Add a row of `false` to align with `pupil_data`
% Copy the data and set outlier values to NaN in `pupil_data`
diam_t0.DiamBothBlankOutliersRemoved = diam_t0.DiamBothBlankBaselineAdjusted;
diam_t0.DiamBothBlankOutliersRemoved(outlier_mask) = NaN;
for i = 1:4
    condition = conditions{i};
    condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
    condition_idx = find(strcmp(diam_t0.Condition, condition));
    outlier_count = 0;
    for j = 1:height(condition_trials)
        current_data = condition_trials.DiamBothBlank(j, :);
        if numel(current_data(~isnan(current_data))) > baseline_min % valid data points
            plot(condition_trials.DiamBothBlankBaselineAdjusted(j,:), 'Color', color_map(condition)); % Plot each trial in light grey        y_limits = [y_limits; ylim];  % Store y-axis limits
        else % fill completely with nans to mark as unusable
            outlier_count = outlier_count+1;
            trial_index_in_result_table = condition_idx(j);
            diam_t0.DiamBothBlankOutliersRemoved(trial_index_in_result_table,:) = NaN;
            continue;  % Skip processing and plotting this trial
        end
    end
    fprintf('Removed %d ouliers in condition %s.\n', outlier_count, condition);
end
title('Individual Trials (Outliers Removed)');
xlabel('Datapoints During BlankScreen'); ylabel('Pupil Diameter (mm)');
hold off;

subplot(num_rows, num_col, 6); hold on;
for i = 1:4
    condition = conditions{i};
    condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
    avg_response = nanmean(condition_trials.DiamBothBlankOutliersRemoved, 1);
    plot(avg_response, 'Color', color_map(condition), 'DisplayName', condition, 'LineWidth', 1.5);
end
hold off;
title('Average per Condition (Outliers Removed)');
xlabel('Datapoints During BlankScreen'); ylabel('Pupil Diameter (Baseline Adjusted, mm)');


%% individual baselines - zscore data
% Loop over each trial and perform z-score normalization individually
for trial = 1:size(diam_t0.DiamBothBlankOutliersRemoved, 2)
    % Extract the data for the current trial
    trial_data = diam_t0.DiamBothBlankOutliersRemoved(:, trial);
    % Calculate the mean and standard deviation, ignoring NaNs
    mean_trial = mean(trial_data, 'omitnan');
    std_trial = std(trial_data, 'omitnan');
    % Perform z-score normalization on non-NaN values
    diam_t0.DiamBothBlankZ(:, trial) = (trial_data - mean_trial) / std_trial;
end

subplot(num_rows, num_col,7); hold on;
for i = 1:4
    condition = conditions{i};
    condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
    for j = 1:height(condition_trials)
        plot(condition_trials.DiamBothBlankZ(j,:), 'Color', color_map(condition)); % Plot each trial in light grey        y_limits = [y_limits; ylim];  % Store y-axis limits
    end
end
title('Individual Trials (z-scored)');
xlabel('Datapoints During BlankScreen'); ylabel('Pupil Diameter (z-score)');
hold off;

subplot(num_rows, num_col, 8); hold on;
for i = 1:4
    condition = conditions{i};
    condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
    avg_response = nanmean(condition_trials.DiamBothBlankZ, 1);
    plot(avg_response, 'Color', color_map(condition), 'DisplayName', condition, 'LineWidth', 1.5);
end
hold off;
title('Average per Condition (z-scored)');
xlabel('Datapoints During BlankScreen'); ylabel('Pupil Diameter (z-score)');





set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf, fullfile(analysis_folder, 'baseline_blank_steps.png'));