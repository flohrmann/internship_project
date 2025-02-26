function result_table = plotPupilDiamPerTrial(result_table, outlier_threshold, color_map, conditions)
% Extract unique conditions from result_table
num_conditions = length(conditions);

% Add a column for storing mean filtered response per condition
result_table.DiamBothStart3StdFilter = nan(height(result_table), size(result_table.DiamBothBlank,2));

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
        plot(condition_trials.DiamBothBlank(j,:), 'Color', [0.8, 0.8, 0.8]); % Plot each trial in light grey        y_limits = [y_limits; ylim];  % Store y-axis limits
    end
    % Plot average response for the current condition in black
    avg_response = nanmean(condition_trials.DiamBothStart, 1);
    plot(avg_response, 'Color', color_map(condition), 'LineWidth', 1.5);
    y_limits = [y_limits; ylim];  % Store y-axis limits
end
title('Individual Trials per Condition (Raw)');
xlabel('Time (s)');
ylabel('Pupil Diameter');
hold off;



%% Plot 2: Average response per condition before outlier removal
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


%% Remove outliers and process data for individual trials
subplot(3, 2, 3);
hold on;
for i = 1:num_conditions
    % Get the current condition's data
    condition = conditions{i};
    condition_idx = find(strcmp(result_table.Condition, condition));
    condition_trials = result_table(strcmp(result_table.Condition, condition), :);
    outlier_count = 0;
    % Plot all trials for the current condition
    for j = 1:size(condition_trials,1)
        current_data = condition_trials.DiamBothBlank(j,:);
        
        if numel(current_data(~isnan(current_data))) > 8
            % Process data: handle NaNs and filter
            current_data = condition_trials.DiamBothBlank(j, :);
            %                 nan_mask = isnan(current_data);
            %                 data_without_nans = current_data(~nan_mask);
            %
            %                 % Apply bandpass filter to clean noise
            %                 [b, a] = butter(2, [0.025, 4] / (60 / 2), 'bandpass');  % 60Hz sampling rate
            %                 filtered_data = filtfilt(b, a, data_without_nans);
            % Outlier removal based on derivative
            derivative_data = diff(current_data);
            mean_derivative = mean(derivative_data, 'omitnan');
            std_derivative = std(derivative_data, 'omitnan');
            outliers = abs(derivative_data) > (mean_derivative + outlier_threshold * std_derivative);
            current_data(2:end) = current_data(2:end) .* ~outliers;
            
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
            %condition_trials.DiamBothStart3StdFilter(trial_index_in_result_table,:) = current_data;
            plot(current_data, 'Color', [0.8, 0.8, 0.8]); % Plot each trial in light grey
        else
            outlier_count = outlier_count+1;
        end
    end
    fprintf('Removed %d ouliers in condition %s.\n', outlier_count, condition);
    
    condition_trials = result_table(strcmp(result_table.Condition, condition), :);
    avg_response = nanmean(condition_trials.DiamBothStart3StdFilter, 1);
    plot(avg_response, 'Color', color_map(condition), 'LineWidth', 1.5);
    y_limits = [y_limits; ylim];  % Store y-axis limits
end
title('Individual Trials per Condition (Without Outliers)');
xlabel('Time (s)');
ylabel('Pupil Diameter');
hold off;




%% Plot 4:average response per condition after outlier removal
subplot(3, 2, 4);
hold on;
for i = 1:num_conditions
    condition = conditions{i};
    condition_trials = result_table(strcmp(result_table.Condition, condition), :);
    avg_response = nanmean(condition_trials.DiamBothStart3StdFilter, 1);
    plot(avg_response, 'Color', color_map(condition), 'DisplayName', condition, 'LineWidth', 1.5);
end

hold off;
title('Average Pupil Response per Condition (Without Outliers)');
xlabel('Time (s)');
ylabel('Pupil Diameter');
legend('show');
y_limit_avg = [y_limit_avg; ylim];  % Store y-axis limits














%% z score trials and average
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
            if ~isnan(current_data(1))
                processed_data_no_nans = fillmissing(current_data, 'spline');
                %normalized_data = zscore(processed_data_no_nans);
            else
                continue;
            end
            %Store processed data back in the current condition trial
            trial_index_in_result_table = condition_idx(j);
            result_table.DiamBothStart3StdFilterZScore(trial_index_in_result_table, :) = processed_data_no_nans;
            %condition_trials.DiamBothStart3StdFilterZScore(trial_index_in_result_table,:) = processed_data_no_nans;
            plot(processed_data_no_nans, 'Color', [0.8, 0.8, 0.8]); % Plot each trial in light grey
    end   
    condition_trials = result_table(strcmp(result_table.Condition, condition), :);
    avg_response = nanmean(condition_trials.DiamBothStart3StdFilterZScore, 1);
    plot(avg_response, 'Color', color_map(condition), 'LineWidth', 1.5);
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

end
