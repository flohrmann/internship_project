function plotDataWithoutNaNs(condition_data_avg_blank, conditions, color_map, analysis_folder, base)
    % Create figure to plot data without NaNs for each condition
    figure('Name', sprintf('Data Without NaNs - %s', base), 'NumberTitle', 'off');
    hold on;
    
    % Define the maximum number of points across conditions for consistent x-axis
    max_points = 0;
    for i = 1:length(conditions)
        condition = conditions{i};
        steps_data = condition_data_avg_blank.(condition).steps_both_median;
        
        % Find the max length of `data_without_nans` across trials in this condition
        max_points = max(max_points, max(cellfun(@(x) length(x.data_without_nans), steps_data)));
    end
    x = 1:max_points;  % Set x-axis range based on the longest trial
    
    % Plot each condition
    for i = 1:length(conditions)
        condition = conditions{i};
        color = color_map(condition);
        steps_data = condition_data_avg_blank.(condition).steps_both_median;
        
        % Collect data without NaNs for all trials
        step_data = nan(length(steps_data), max_points);  % Initialize with NaNs for consistent lengths
        for j = 1:length(steps_data)
            trial_data = steps_data{j}.data_without_nans;
            step_data(j, 1:length(trial_data)) = trial_data;  % Fill trial data up to its length
        end
        
        % Plot individual trials for the condition
        for j = 1:size(step_data, 1)
            plot(x, step_data(j, :), 'Color', [color, 0.3]);  % Semi-transparent individual trials
        end
        
        % Plot mean line for this condition
        mean_data = nanmean(step_data, 1);
        plot(x, mean_data, 'Color', color, 'LineWidth', 2, 'DisplayName', condition);
    end
    
    % Add labels and legend
    xlabel('Time (arbitrary units)');
    ylabel('Pupil Diameter (data without NaNs)');
    title(['Data Without NaNs - ', base]);
    legend('show', 'Location', 'northeast');
    hold off;
    
    % Save the figure
    saveas(gcf, fullfile(analysis_folder, sprintf('data_without_nans_%s.png', base)));
end
