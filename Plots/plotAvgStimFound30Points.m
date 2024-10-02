function plotAvgStimFound30Points(avg_stim_30, conditions, id, analysis_folder, color_map)

% Define the custom x-axis: 30 points, with the last point being 0 and each earlier point spaced by -0.016 seconds
num_points = 30;
x_values = linspace(-0.016 * (num_points - 1), 0, num_points);  % Time points

%% Plot 1: Mean with Standard Error
figure;
hold on;

% Loop through each condition and compute the average of the 30 data points
for i = 1:length(conditions)
    condition = conditions{i};  % Get the current condition
    condition_mask = strcmp(avg_stim_30.Condition, condition);  % Find rows for this condition
    data_points = avg_stim_30.DataPoints(condition_mask,:);  % Extract the data points for this condition
    
    num_trials = length(data_points);  % Number of trials
    % Convert cell array into a matrix of size [num_trials x num_points]
    data_matrix = nan(num_trials, num_points);  % Initialize with NaNs to handle missing values
    for j = 1:num_trials
        data_matrix(j, :) = data_points{j};  % Fill each row with the 30 data points from each trial
    end
    
    % Calculate the mean and standard error across trials for the current condition
    mean_values = mean(data_matrix, 1, 'omitnan');  % Mean across trials
    std_error = std(data_matrix, 0, 1, 'omitnan') / sqrt(sum(condition_mask));  % Standard error
    
    % Plot the mean line for this condition
    color = color_map(condition);  % Get the color for the condition
    plot(x_values, mean_values, 'Color', color, 'LineWidth', 2);  % Plot the mean

    % Plot the shaded error region (mean +/- standard error)
    fill([x_values fliplr(x_values)], [mean_values + std_error fliplr(mean_values - std_error)], ...
         color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Shaded region with transparency
end

% X-axis adjustments: remove white space on the left
xlim([-0.016 * (num_points - 1), 0]);
xlabel('Time Before Button Press (s)');
ylabel('Mean Pupil Diameter and Standard Error (mm)');
title(['Mean Pupil Diameter Before Button Press - Participant ', num2str(id)]);
% in legend only show the lines (not the shaded areas)
%legend('a', 'a SE', 'a simple', 'a simple SE','b','b SE', 'b simple', 'b simple SE', 'Location', 'northeastoutside', 'Orientation', 'vertical');
legend_objects = findobj(gca, 'Type', 'Line');  % Find all line objects in the plot
legend(legend_objects(end:-1:1), {'a', 'a simple', 'b', 'b simple'}, ...
       'Location', 'northwest', 'Orientation', 'vertical');
legend('boxoff')      
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1])
saveas(gcf, fullfile(analysis_folder, 'button_press_mean_30_points_by_condition_se.png'));
hold off;

%% Plot 2: Median with Standard Error
figure;
hold on;

% Loop through each condition and compute the median of the 30 data points
for i = 1:length(conditions)
    condition = conditions{i};  % Get the current condition
    condition_mask = strcmp(avg_stim_30.Condition, condition);  % Find rows for this condition
    data_points = avg_stim_30.DataPoints(condition_mask,:);  % Extract the data points for this condition
    
    num_trials = length(data_points);  % Number of trials
    % Convert cell array into a matrix of size [num_trials x num_points]
    data_matrix = nan(num_trials, num_points);  % Initialize with NaNs to handle missing values
    for j = 1:num_trials
        data_matrix(j, :) = data_points{j};  % Fill each row with the 30 data points from each trial
    end
    
    % Calculate the median and standard error across trials for the current condition
    median_values = median(data_matrix, 1, 'omitnan');  % Median across trials
    std_error = std(data_matrix, 0, 1, 'omitnan') / sqrt(sum(condition_mask));  % Standard error
    
    % Plot the median line for this condition
    color = color_map(condition);  % Get the color for the condition
    plot(x_values, median_values, 'Color', color, 'LineWidth', 2);  % Plot the median

    % Plot the shaded error region (median +/- standard error)
    fill([x_values fliplr(x_values)], [median_values + std_error fliplr(median_values - std_error)], ...
         color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Shaded region with transparency
end

% X-axis adjustments: remove white space on the left
xlim([-0.016 * (num_points - 1), 0]);
xlabel('Time Before Button Press (s)');
ylabel('Median Pupil Diameter and Standard Error (mm)');
title(['Median Pupil Diameter Before Button Press - Participant ', num2str(id)]);
% in legend only show the lines (not the shaded areas)
%legend('a', 'a SE', 'a simple', 'a simple SE','b','b SE', 'b simple', 'b simple SE', 'Location', 'northeastoutside', 'Orientation', 'vertical');
legend_objects = findobj(gca, 'Type', 'Line');  % Find all line objects in the plot
legend(legend_objects(end:-1:1), {'a', 'a simple', 'b', 'b simple'}, ...
       'Location', 'northwest', 'Orientation', 'vertical');
legend('boxoff')  
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1])
saveas(gcf, fullfile(analysis_folder, 'button_press_median_30_points_by_condition_se.png'));
hold off;




%% Plot 3: Median with IQR (Interquartile Range)
% shaded region now represents the interquartile range (IQR), which captures 
% the variability between the 25th percentile (Q1) and the 75th percentile (Q3).

figure;
hold on;

% Loop through each condition and compute the median of the 30 data points
for i = 1:length(conditions)
    condition = conditions{i};  % Get the current condition
    condition_mask = strcmp(avg_stim_30.Condition, condition);  % Find rows for this condition
    data_points = avg_stim_30.DataPoints(condition_mask,:);  % Extract the data points for this condition
    
    num_trials = length(data_points);  % Number of trials
    % Convert cell array into a matrix of size [num_trials x num_points]
    data_matrix = nan(num_trials, num_points);  % Initialize with NaNs to handle missing values
    for j = 1:num_trials
        data_matrix(j, :) = data_points{j};  % Fill each row with the 30 data points from each trial
    end
    
    % Calculate the median and IQR across trials for the current condition
    median_values = median(data_matrix, 1, 'omitnan');  % Median across trials
    lower_quartile = prctile(data_matrix, 25, 1);  % 25th percentile (Q1)
    upper_quartile = prctile(data_matrix, 75, 1);  % 75th percentile (Q3)
    
    % Plot the median line for this condition
    color = color_map(condition);  % Get the color for the condition
    plot(x_values, median_values, 'Color', color, 'LineWidth', 2);  % Plot the median
    
    % Plot the shaded region for the interquartile range (IQR)
    fill([x_values fliplr(x_values)], [upper_quartile fliplr(lower_quartile)], ...
         color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Shaded region with transparency
end

% X-axis adjustments: remove white space on the left
xlim([-0.016 * (num_points - 1), 0]);
xlabel('Time Before Button Press (s)');
ylabel('Median Pupil Diameter and Interquartile Range Q1-Q3 (mm)');
title(['Median Pupil Diameter Before Button Press - Participant ', num2str(id)]);
% Adjust the legend to only show lines
legend_objects = findobj(gca, 'Type', 'Line');  % Find all line objects in the plot
legend(legend_objects(end:-1:1), {'a', 'a simple', 'b', 'b simple'}, ...
       'Location', 'northwest', 'Orientation', 'vertical');
legend('boxoff')    
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1])
saveas(gcf, fullfile(analysis_folder, 'button_press_median_30_points_by_condition_iqr.png'));
hold off;

end
