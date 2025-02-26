function plotAvgPupilDiamsBeforeAfterStimFoundByGroupByCondition(average, y_title, mean_diam_around_stim_adhd, mean_diam_around_stim_nonadhd, conditions, color_map, num_before, num_after, outlier_threshold, comparison_results_folder)
    % Create a tiled layout for the 4 plots and 1 legend tile
    figure;
    %t = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); % 4 tiles for plots, 1 for the legend
    t = tiledlayout('flow','TileSpacing','compact');
    title(t, [average, ' Pupil Diameter Per Condition Before and After Looking at the Target'], 'FontSize', 14, 'FontWeight', 'bold');
    condition_names = {'a', 'a simple', 'b', 'b simple'};
    
    for condIdx = 1:length(conditions)
        condition = conditions{condIdx};  
        % Extract ADHD group data
        adhd_data_before = mean_diam_around_stim_adhd.(condition).before;  
        adhd_data_after = mean_diam_around_stim_adhd.(condition).after;   
        
        % Extract non-ADHD group data
        nonadhd_data_before = mean_diam_around_stim_nonadhd.(condition).before; 
        nonadhd_data_after = mean_diam_around_stim_nonadhd.(condition).after;    

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
        x_before = linspace(-0.016 * num_before, 0, num_before);  % 30 data points before
        x_after = linspace(0.016, 0.016 * num_after, num_after);  % 10/20 data points after
        x_values = [x_before, x_after];  % Combine x-values before and after target found

        % Plot ADHD group data
        plot(x_values, [avg_adhd_before, avg_adhd_after], 'Color', color_map('ADHD'), 'MarkerFaceColor', color_map('ADHD'), 'LineWidth', 2);
        % Plot the shaded error region (mean +/- standard error) for ADHD
        fill([x_values, fliplr(x_values)], ...
             [avg_adhd_before, avg_adhd_after, fliplr([avg_adhd_before + std_error_adhd_before, avg_adhd_after + std_error_adhd_after])], ...
             color_map('ADHD'), 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Upper bound
        fill([x_values, fliplr(x_values)], ...
             [avg_adhd_before - std_error_adhd_before, avg_adhd_after - std_error_adhd_after, fliplr([avg_adhd_before, avg_adhd_after])], ...
             color_map('ADHD'), 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Lower bound

        % Plot Non-ADHD group data
        plot(x_values, [avg_nonadhd_before, avg_nonadhd_after], 'Color', color_map('nonADHD'), 'MarkerFaceColor', color_map('nonADHD'), 'LineWidth', 2);
        % Plot the shaded error region (mean +/- standard error) for Non-ADHD
        fill([x_values, fliplr(x_values)], ...
             [avg_nonadhd_before, avg_nonadhd_after, fliplr([avg_nonadhd_before + std_error_nonadhd_before,  avg_nonadhd_after + std_error_nonadhd_after])], ...
             color_map('nonADHD'), 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Upper bound
        fill([x_values, fliplr(x_values)], ...
             [avg_nonadhd_before - std_error_nonadhd_before, avg_nonadhd_after - std_error_nonadhd_after, fliplr([avg_nonadhd_before, avg_nonadhd_after])], ...
             color_map('nonADHD'), 'FaceAlpha', 0.2, 'EdgeColor', 'none');  % Lower bound

        % Add a vertical line at x = 0 (time point where the target is found)
        plot([0, 0], ylim, '--k', 'LineWidth', 1.5, 'DisplayName', 'Gaze Reached Target');
        % X-axis adjustments: remove white space on the left
        xlim([-0.016 * num_before, 0.016 * num_after]);


        % Add labels and title for each condition
        xlabel('Time (s)');
        ylabel(y_title);
        title(condition_names(condIdx));

        hold off;
    end

    %% Add the legend and global title outside the subplots
    legend_objects = findobj(gca, 'Type', 'Line');  % Find all line objects in the plot

    lgd = legend(legend_objects(end:-1:1), {'ADHD', 'nonADHD', 'Gaze Reaches Target',}, ...
        'Orientation', 'horizontal');
    lgd.NumColumns = 1;
    lgd.Layout.Tile = 5;
    lgd.Location = 'bestoutside';

    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1])
    saveas(gcf, fullfile(comparison_results_folder, strcat('norm_pupil_diameter_comparison_by_condition_group_outlier_threshold_', outlier_threshold, '.png')));
end
