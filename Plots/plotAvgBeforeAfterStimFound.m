function plotAvgBeforeAfterStimFound(result_table, conditions, id, analysis_folder, color_map, num_before, num_after)

total_points = num_before + num_after;

% Create x_values with 30 points before (-0.016 * (num_before - 1) to 0) and 10 points after (0 to 0.016 * num_after)
x_values = [-0.016 * (num_before - 1):0.016:0, 0.016 * (1:num_after)];

% Loop through each condition and compute the mean and median for the 30 data points before and 10 data points after
for method = {'Mean', 'Median'}  % Create separate plots for mean and median
    method_name = method{1};
    
    figure; % Create a new figure for each method (mean/median)
    hold on;
    
    for i = 1:length(conditions)
        condition = conditions{i};  % Get the current condition
        condition_mask = strcmp(result_table.Condition, condition);  % Find rows for this condition
        data_before = result_table.DataPointsBefore(condition_mask,:);  % Extract the 30 data points before
        data_after = result_table.DataPointsAfter(condition_mask,:);  % Extract the 10 data points after
        
        % Combine the before and after data
        data_combined = [data_before, data_after];
        
        num_trials = size(data_combined, 1);  % Number of trials
        
        % Convert the table into a matrix for easier averaging
        data_matrix = nan(num_trials, total_points);  % Initialize with NaNs to handle missing values
        for j = 1:num_trials
            data_matrix(j, :) = data_combined(j, :);  % Fill each row with the 30 + 10 data points from each trial
        end
        
        % Calculate the mean or median across trials for the current condition
        if strcmp(method_name, 'Mean')
            values = mean(data_matrix, 1, 'omitnan');  % Average across trials
            std_error = std(data_matrix, 0, 1, 'omitnan') / sqrt(num_trials);  % Standard error for mean
        else
            values = median(data_matrix, 1, 'omitnan');  % Median across trials
            std_error = std(data_matrix, 0, 1, 'omitnan') / sqrt(num_trials);  % Standard error for median
        end
        
        % Plot the average/median line for this condition
        color = color_map(condition);  % Get the color for the condition
        plot(x_values, values, 'Color', color, 'LineWidth', 2);  % Plot the average/median
        
        % Plot the shaded error region (mean/median +/- standard error)
        fill([x_values, fliplr(x_values)], [values + std_error, fliplr(values - std_error)], ...
            color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Shaded region with transparency
    end
    
    % Plot the grey dashed line at 0 (gaze reaches target)
    xline(0, '--', 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5]);
    % hidden plot just for legend
    h = plot([0 0], ylim, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5, 'Visible', 'off');
    % X-axis adjustments: remove white space on the left
    xlim([-0.016 * (num_before - 1), 0.016 * num_after]);
    
    xlabel('Time (s) before and after target found');
    ylabel([method_name,' Pupil Diameter (mm)']);
    used_trials = height(result_table);
    
%     if id == 0 % for groups (adhd/nonadhd)
%         title('Average Pupil Diameter Once Gaze Reaches Target - ADHD');
%     elseif id == 100
%         title('Average Pupil Diameter Once Gaze Reaches Target - nonADHD');
%     else % for individual participants
        title(['Average Pupil Diameter Once Gaze Reaches Target - Participant ', num2str(id), ...
            ' (', num2str(used_trials), ' trials)']);
%     end
    
    % Manually set the legend with the line and the shaded error regions
    legend_objects = findobj(gca, 'Type', 'Line');  % Find all line objects in the plot
    legend(legend_objects(end:-1:1), 'a', 'a simple', 'b', 'b simple', 'Gaze Reached Target', ...
        'Location', 'northwest', 'Orientation', 'vertical');
    legend('boxoff')
    
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf, fullfile(analysis_folder, [method_name, 'diam_before_after_stim_found_by_condition.png']));
    hold off;
end
end
