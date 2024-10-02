function slope_table = analyzeSlopesByCondition(avg_stim_30, conditions, id, analysis_folder, color_map)
    % Initialize a table to store slopes for each trial and condition
    slope_table = table('Size', [0, 4], ...
                        'VariableTypes', {'string', 'double', 'double', 'double'}, ...
                        'VariableNames', {'Condition', 'Trial', 'Mean_Slope', 'Median_Slope'});

    % Prepare two figures: one for mean slopes and one for median slopes
    figure_mean = figure;
    hold on;
    title(['Change in Pupil Diameter Before Button Press - Participant ', num2str(id)]);
    xlabel('Time Before Button Press [s]');
    ylabel('Mean Pupil Diameter (mm)');

    figure_median = figure;
    hold on;
    title(['Change in Pupil Diameter Before Button Press - Participant ', num2str(id)]);
    xlabel('Time Before Button Press (s)');
    ylabel('Median Pupil Diameter (mm)');

    % Loop through each condition
    for i = 1:length(conditions)
        condition = conditions{i};  % Get the current condition
        condition_mask = strcmp(avg_stim_30.Condition, condition);  % Find rows for this condition
        data_points = avg_stim_30.DataPoints(condition_mask,:);  % Extract the data points for this condition
        
        % Find the original trial numbers (row numbers in avg_stim_30)
        trial_numbers = find(condition_mask);

        num_trials = length(data_points);  % Number of trials
        % Convert cell array into a matrix of size [num_trials x num_points]
        data_matrix = nan(num_trials, 30);  % Initialize with NaNs to handle missing values
        for j = 1:num_trials
            data_matrix(j, :) = data_points{j};  % Fill each row with the 30 data points from each trial
        end

        % Define the x-values for the 30 time points (time before stimulus)
        x_values = -0.016 * (29:-1:0);  % Time from -0.16 * 29 to 0

        % Loop through each trial to calculate the individual mean and median slopes
        for trial_idx = 1:num_trials
            % Fit a linear regression model to the current trial data (mean)
            p_mean_trial = polyfit(x_values, data_matrix(trial_idx, :), 1);  % Linear fit (1st-degree polynomial)

            % Extract the slope (p_mean_trial(1) is the slope, p_mean_trial(2) is the intercept)
            mean_slope_trial = p_mean_trial(1);

            % Fit a linear regression model to the current trial data (median)
            median_values_trial = median(data_matrix, 1, 'omitnan');  % Median across time points
            p_median_trial = polyfit(x_values, median_values_trial, 1);  % Linear fit for median

            % Extract the median slope (p_median_trial(1) is the slope, p_median_trial(2) is the intercept)
            median_slope_trial = p_median_trial(1);

            % Add the condition, actual trial number, and slopes to the slope_table
            slope_table = [slope_table; {condition, trial_numbers(trial_idx), mean_slope_trial, median_slope_trial}];
        end

        % Calculate the mean and median across trials for the current condition
        mean_values = mean(data_matrix, 1, 'omitnan');  % Mean across trials
        median_values = median(data_matrix, 1, 'omitnan');  % Median across trials
        
        % Fit a linear regression model to the mean and median values for plotting
        p_mean = polyfit(x_values, mean_values, 1);  % Linear fit for mean
        p_median = polyfit(x_values, median_values, 1);  % Linear fit for median

        % Get the color for the condition
        color = color_map(condition);

        % Plot the mean values and linear fit on the mean figure
        figure(figure_mean);
        plot(x_values, mean_values, 'Color', color, 'LineWidth', 2);  % Plot the mean
        plot(x_values, polyval(p_mean, x_values), '--', 'Color', color, 'LineWidth', 1);  % Plot the mean linear fit

        % Plot the median values and linear fit on the median figure
        figure(figure_median);
        plot(x_values, median_values, 'Color', color, 'LineWidth', 2);  % Plot the median
        plot(x_values, polyval(p_median, x_values), '--', 'Color', color, 'LineWidth', 1);  % Plot the median linear fit
    end

    % Add legends and save the figures
    figure(figure_mean);
    legend('a', 'a slope', 'a simple', 'a simple slope', ...
           'b', 'b slope', 'b simple', 'b simple slope', ...
           'Location', 'northeast');
    legend('boxoff');
    xlim([x_values(1), x_values(end)]);
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1])
    saveas(figure_mean, fullfile(analysis_folder, 'button_press_pupil_diameter_mean_slopes_by_condition.png'));

    figure(figure_median);
    legend('a', 'a slope', 'a simple', 'a simple slope', ...
           'b', 'b slope', 'b simple', 'b simple slope', ...
           'Location', 'northeast');
    legend('boxoff');
    xlim([x_values(1), x_values(end)]);
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure_median, fullfile(analysis_folder, 'button_press_pupil_diameter_median_slopes_by_condition.png'));

    % Save the slope data to a .mat file
    save(fullfile(analysis_folder, 'slope_table.mat'), 'slope_table');
    hold off;
end
