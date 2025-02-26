function plotPupilMaxStimFound(avg_adhd, avg_nonadhd, conditions, color_map, comparison_results_folder)
    % Create a new figure for the subplots
    figure;
    t = tiledlayout(1, length(conditions), 'TileSpacing', 'compact', 'Padding', 'compact');
    
    condition_names = {'a', 'a simple', 'b', 'b simple'};
    pupil_fields = {'a_max', 'as_max', 'b_max', 'bs_max'};
    
    % Loop through each condition and plot the max pupil dilation for each group
    for i = 1:length(conditions)
        condition = condition_names{i};  % Get the current condition
        max_diam_adhd = avg_adhd.(pupil_fields{i});
        max_diam_nonadhd = avg_nonadhd.(pupil_fields{i});
        
        % Calculate the means for each group
        mean_adhd = mean(max_diam_adhd, 'omitnan');
        mean_nonadhd = mean(max_diam_nonadhd, 'omitnan');
        
        % Plot the scatter plot and mean lines
        nexttile;
        hold on;
        sz = 70;
        
        % Scatter plots for individual participants; x value fake/just for
        % plotting groups seperately
        scatter(ones(size(max_diam_adhd)), max_diam_adhd, sz, color_map('ADHD'), 'filled');
        scatter(2 * ones(size(max_diam_nonadhd)), max_diam_nonadhd, sz, color_map('nonADHD'), 'filled');
        
        % Plot mean lines
        plot([0.8, 1.2], [mean_adhd, mean_adhd], 'Color', color_map('ADHD'), 'LineWidth', 4);  % Mean ADHD
        plot([1.8, 2.2], [mean_nonadhd, mean_nonadhd], 'Color', color_map('nonADHD'), 'LineWidth', 4);  % Mean non-ADHD

        xticks([1 2]);
        xticklabels({'ADHD', 'Non-ADHD'});
        xlim([0.5 2.5]);
        title(condition);
        ylabel('Pupil Diameter (z-score)');
        grid off;
        hold off;
    end
    title(t, 'Max Pupil Diameter After Gaze Reached Target', 'FontWeight', 'bold');
    saveas(gcf, fullfile(comparison_results_folder, 'pupil_max_target_found_comparison_small.png'));
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf, fullfile(comparison_results_folder, 'pupil_max_target_found_comparison.png'));
end
