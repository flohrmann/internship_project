function plotAvgPupilDiamsBeforeAfterStimFoundByGroupByCondition(mean_diam_around_stim_adhd, mean_diam_around_stim_nonadhd, conditions, analysis_folder)
    % Create a tiled layout for the 4 plots and 1 legend tile
    figure;
    t = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); % 4 tiles for plots, 1 for the legend

    % Loop through each condition
    for condIdx = 1:length(conditions)
        condition = conditions{condIdx};  % Current condition

        % Extract ADHD group data
        adhd_data_before = mean_diam_around_stim_adhd.(condition).before;  % [n_participants x 30]
        adhd_data_after = mean_diam_around_stim_adhd.(condition).after;    % [n_participants x 10]
        
        % Extract non-ADHD group data
        nonadhd_data_before = mean_diam_around_stim_nonadhd.(condition).before;  % [n_participants x 30]
        nonadhd_data_after = mean_diam_around_stim_nonadhd.(condition).after;    % [n_participants x 10]

        % Compute the average for each group and condition
        avg_adhd_before = mean(adhd_data_before, 1, 'omitnan');
        avg_adhd_after = mean(adhd_data_after, 1, 'omitnan');
        avg_nonadhd_before = mean(nonadhd_data_before, 1, 'omitnan');
        avg_nonadhd_after = mean(nonadhd_data_after, 1, 'omitnan');

        % Compute the standard error for the shaded error regions
        std_error_adhd_before = std(adhd_data_before, 0, 1, 'omitnan') / sqrt(size(adhd_data_before, 1));
        std_error_adhd_after = std(adhd_data_after, 0, 1, 'omitnan') / sqrt(size(adhd_data_after, 1));
        std_error_nonadhd_before = std(nonadhd_data_before, 0, 1, 'omitnan') / sqrt(size(nonadhd_data_before, 1));
        std_error_nonadhd_after = std(nonadhd_data_after, 0, 1, 'omitnan') / sqrt(size(nonadhd_data_after, 1));

        % Plot the results for both groups per condition in the current tile
        nexttile;
        hold on;
        % X-axis: time points (before and after the stimulus)
        x_before = linspace(-0.016 * 29, 0, 30);  % 30 data points before
        x_after = linspace(0.016, 0.016 * 10, 10);  % 10 data points after
        x_values = [x_before, x_after];  % Combine x-values before and after target found

        % Plot ADHD group data
        plot(x_values, [avg_adhd_before, avg_adhd_after], 'r', 'LineWidth', 2);
        % Plot the shaded error region (mean +/- standard error) for ADHD
        fill([x_values, fliplr(x_values)], ...
             [avg_adhd_before, avg_adhd_after, fliplr([avg_adhd_before + std_error_adhd_before, avg_adhd_after + std_error_adhd_after])], ...
             'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Upper bound
        fill([x_values, fliplr(x_values)], ...
             [avg_adhd_before - std_error_adhd_before, avg_adhd_after - std_error_adhd_after, fliplr([avg_adhd_before, avg_adhd_after])], ...
             'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Lower bound

        % Plot Non-ADHD group data
        plot(x_values, [avg_nonadhd_before, avg_nonadhd_after], 'b', 'LineWidth', 2);
        % Plot the shaded error region (mean +/- standard error) for Non-ADHD
        fill([x_values, fliplr(x_values)], ...
             [avg_nonadhd_before, avg_nonadhd_after, fliplr([avg_nonadhd_before + std_error_nonadhd_before,  avg_nonadhd_after + std_error_nonadhd_after])], ...
             'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Upper bound
        fill([x_values, fliplr(x_values)], ...
             [avg_nonadhd_before - std_error_nonadhd_before, avg_nonadhd_after - std_error_nonadhd_after, fliplr([avg_nonadhd_before, avg_nonadhd_after])], ...
             'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Lower bound

        % Add a vertical line at x = 0 (time point where the target is found)
        plot([0, 0], ylim, '--k', 'LineWidth', 1.5, 'DisplayName', 'Gaze Reached Target');

        % Add labels and title for each condition
        xlabel('Time (s)');
        ylabel('Average Pupil Diameter (mm)');
        title(['Condition: ', condition]);

        hold off;
    end

    % Add a tile for the legend
    nexttile;
    hold on;
    % Create the legend
    legend({'ADHD', 'Non-ADHD', 'Gaze Reached Target'}, 'Location', 'best', 'Orientation', 'vertical');

    % Save the plot
    saveas(gcf, fullfile(analysis_folder, 'pupil_diameter_comparison_by_condition_group.png'));
end
