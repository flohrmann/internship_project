function result_table = plotPupilDiamPerStartBaseline(result_table, outlier_threshold, baseline_min_length, color_map, conditions, analysis_folder)
% Extract unique conditions from result_table
num_conditions = length(conditions);

% Add a column for storing mean filtered response per condition
result_table.DiamBothStart3StdFilter = nan(height(result_table), size(result_table.DiamBothStart,2));
result_table.DiamBothStart3StdFilterMean = nan(height(result_table), 1);
result_table.DiamBothStart3StdFilterMedian = nan(height(result_table), 1);

figure;
%% Plot 1: Individual trials per condition before outlier removal
subplot(3, 2, 1);
hold on;
y_limits = [];  % Track y-axis limits for alignment
for i = 1:num_conditions
    % Get the current condition's data
    condition = conditions{i};
    condition_trials = result_table(strcmp(result_table.Condition, condition), :);
    % Plot all trials for the current condition
    for j = 1:height(condition_trials)
        plot(condition_trials.DiamBothStart(j,:), 'Color', color_map(condition)); % Plot each trial in light grey        y_limits = [y_limits; ylim];  % Store y-axis limits
    end
    % Plot average response for the current condition in black
%     avg_response = nanmean(condition_trials.DiamBothBlank, 1);
%     plot(avg_response, 'Color', color_map(condition), 'LineWidth', 1.5);
     y_limits = [y_limits; ylim];  % Store y-axis limits
end
title('Individual Trials per Condition (Raw)');
xlabel('Time (s)');
ylabel('Pupil Diameter');
hold off;



%% Average response per condition before outlier removal
subplot(3, 2, 2);
hold on;
y_limit_avg = [];  % Track y-axis limits for avg data
for i = 1:num_conditions
    condition = conditions{i};
    condition_trials = result_table(strcmp(result_table.Condition, condition), :);
    avg_response = nanmean(condition_trials.DiamBothStart, 1);
    plot(avg_response, 'Color', color_map(condition), 'DisplayName', condition, 'LineWidth', 1.5);
end
hold off;
title('Average Pupil Response per Condition (Raw)');
xlabel('Time (s)');
ylabel('Pupil Diameter');
legend('show');
y_limit_avg = [y_limit_avg; ylim];  % Store y-axis limits
%y_limits = [y_limits; ylim];



%% 2. Remove outliers and process data for individual trials
subplot(3, 2, 3);
hold on;
for i = 1:num_conditions
    condition = conditions{i};
    condition_idx = find(strcmp(result_table.Condition, condition));
    condition_trials = result_table(strcmp(result_table.Condition, condition), :);
    outlier_count = 0;
    for j = 1:size(condition_trials,1)
        current_data = condition_trials.DiamBothStart(j,:);
        if numel(current_data(~isnan(current_data))) > baseline_min_length % valid data points
            current_data = condition_trials.DiamBothStart(j, :);
            
            %Outlier removal based on derivative
            derivative_data = diff(current_data);
            mean_derivative = mean(derivative_data, 'omitnan');
            std_derivative = std(derivative_data, 'omitnan');
            outliers = abs(derivative_data) > (mean_derivative + outlier_threshold * std_derivative);
            current_data(2:end) = current_data(2:end) .* ~outliers;
%             z_scores = abs((current_data - mean(current_data, 'omitnan')) ./ std(current_data, 'omitnan'));
%             outliers = z_scores > outlier_threshold;  % outlier_threshold = 3 for 3 standard deviations
%                         current_data = current_data .* ~outliers;

            if any(outliers)% If any outliers are found, skip this trial
                outlier_count = outlier_count+1;
                continue;  % Skip processing and plotting this trial
            end
            
            % Interpolate and normalize
            %processed_data_no_outliers = fillmissing(current_data, 'spline');
            %normalized_data = zscore(processed_data_no_outliers);
            
            % Store processed data back in the current condition trial
            %                 processed_data = nan(size(current_data));
            %                 processed_data(~nan_mask) = normalized_data;
            
            
            % Store the processed data back in the corresponding row of result_table
            trial_index_in_result_table = condition_idx(j);
            result_table.DiamBothStart3StdFilter(trial_index_in_result_table, :) = current_data;
            result_table.DiamBothStart3StdFilterMean(trial_index_in_result_table, :) = mean(current_data, 'omitnan');
            result_table.DiamBothStart3StdFilterMedian(trial_index_in_result_table, :) = median(current_data, 'omitnan');
            %condition_trials.DiamBothStart3StdFilter(trial_index_in_result_table,:) = current_data;
            plot(current_data, 'Color', color_map(condition)); % Plot each trial in light grey
        else
            outlier_count = outlier_count+1;
        end
    end
    fprintf('Removed %d ouliers in condition %s.\n', outlier_count, condition);
%     condition_trials = result_table(strcmp(result_table.Condition, condition), :);
%     avg_response = nanmean(condition_trials.DiamBothBlank3StdFilter, 1);
%     plot(avg_response, 'Color', color_map(condition), 'LineWidth', 1.5);
    y_limits = [y_limits; ylim];  % Store y-axis limits
end
title(strcat('Individual Trials per Condition (Without Outliers ',num2str(outlier_threshold),' SE))'));
xlabel('Time (s)');
ylabel('Pupil Diameter');
hold off;



%% average response per condition after outlier removal
subplot(3, 2, 4);
hold on;
for i = 1:num_conditions
    condition = conditions{i};
    condition_trials = result_table(strcmp(result_table.Condition, condition), :);
    avg_response = nanmean(condition_trials.DiamBothStart3StdFilter, 1);
    plot(avg_response, 'Color', color_map(condition), 'DisplayName', condition, 'LineWidth', 1.5);
end

hold off;
title(strcat('Average Pupil Response per Condition (Without Outliers ',num2str(outlier_threshold),' SE))'));
xlabel('Time (s)');
ylabel('Pupil Diameter');
legend('show');
y_limit_avg = [y_limit_avg; ylim];  % Store y-axis limits
%y_limits = [y_limits; ylim];







%% 3. z score trials and average
subplot(3, 2, 5);
hold on;
for i = 1:num_conditions
    % Get the current condition's data
    condition = conditions{i};    
    condition_idx = find(strcmp(result_table.Condition, condition));
    condition_trials = result_table(strcmp(result_table.Condition, condition), :);
    
    % Plot all trials for the current condition
    for j = 1:size(condition_trials,1)

        current_data = condition_trials.DiamBothStart3StdFilter(j,:);
            % Interpolate and normalize
            %if ~isnan(current_data(1))
                %processed_data_no_nans = fillmissing(current_data, 'spline');
            %    normalized_data = zscore(current_data);
            %else
            %    continue;
            %end
            normalized_data = zscore(current_data);
            %Store processed data back in the current condition trial
            trial_index_in_result_table = condition_idx(j);
            result_table.DiamBothStart3StdFilterZScore(trial_index_in_result_table, :) = normalized_data;
            %condition_trials.DiamBothStart3StdFilterZScore(trial_index_in_result_table,:) = processed_data_no_nans;
            plot(normalized_data, 'Color', color_map(condition)); % Plot each trial in light grey
    end   
%     condition_trials = result_table(strcmp(result_table.Condition, condition), :);
%     avg_response = nanmean(condition_trials.DiamBothBlank3StdFilterZScore, 1);
%     plot(avg_response, 'Color', color_map(condition), 'LineWidth', 1.5);
    %y_limits = [y_limits; ylim];  % Store y-axis limits
end
title('Individual Trials per Condition (Without Outliers, z-scored)');
xlabel('Time (s)');
ylabel('Pupil Diameter (z-score)');
hold off;



%% Normalized average response per condition after outlier removal
subplot(3, 2, 6);
hold on;
for i = 1:num_conditions
    condition = conditions{i};
    condition_trials = result_table(strcmp(result_table.Condition, condition), :);
    avg_response = nanmean(condition_trials.DiamBothStart3StdFilterZScore, 1);
    plot(avg_response, 'Color', color_map(condition), 'DisplayName', condition, 'LineWidth', 1.5);
end

hold off;
title('Average Pupil Response per Condition (Without Outliers, z-scored)');
xlabel('Time (s)');
ylabel('Pupil Diameter (z-score)');
legend('show');
%y_limit_avg = [y_limit_avg; ylim];  % Store y-axis limits








sgtitle('Pupil Response During Baseline (Start of Stimulation Screen)')



%% Apply consistent y-axis limits for aligned plots
% For raw data plots (subplot 1 and 2)
y_limit_common = [min(y_limits(:)), max(y_limits(:))];
subplot(3, 2, 1); ylim(y_limit_common);
subplot(3, 2, 3); ylim(y_limit_common);
%subplot(3, 3, 5); ylim(y_limit_common);

% For averaged plots data plots (subplot 3 and 4)
y_limit_common_filtered = [min(y_limit_avg(:)), max(y_limit_avg(:))];
subplot(3, 2, 2); ylim(y_limit_common_filtered);
subplot(3, 2, 4); ylim(y_limit_common_filtered);
%subplot(3, 3, 6); ylim(y_limit_common_filtered);


set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf, fullfile(analysis_folder, sprintf('start_baseline_pupil_diam_filter.png')));


end
