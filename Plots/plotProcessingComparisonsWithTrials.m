function plotProcessingComparisonsWithTrials(condition_data_avg_blank, conditions, analysis_folder, color_map, sr, base)
    % Create figure with 6 rows for each processing stage: without NaNs, corrected, filtered, interpolated, normalized, final
    figure('Name', sprintf('Processing Steps Comparison with Individual Trials [%s]', base), 'NumberTitle', 'off');
    
    % Define number of time points for the x-axis
    num_points = length(condition_data_avg_blank.a.steps_both_median{1}.data_without_nans);
    x_values = linspace(-sr * (num_points - 1), 0, num_points);

    % Define processing steps and manually extract them
    titles = {'Data Without NaNs', 'Baseline Corrected Data', 'Filtered Data', ...
              'Interpolated Data', 'Normalized Data', 'Final Processed Data'};
    y_labels = {'Raw Pupil', 'Baseline Corrected Pupil', 'Filtered Pupil', ...
                'Interpolated Pupil', 'Normalized Pupil', 'Final Pupil'};

    % Loop over each processing step and generate plots
    for row_idx = 1:length(titles)
        switch row_idx
            case 1
                step_field = 'data_without_nans';
            case 2
                step_field = 'data_corrected';
            case 3
                step_field = 'filtered_data';
            case 4
                step_field = 'interpolated_data';
            case 5
                step_field = 'normalized_data';
            case 6
                step_field = 'final_data';
        end
        
        % Column 1: Mean ± SEM of the processing step
        subplot(length(titles), 2, (row_idx - 1) * 2 + 1);
        hold on;
        for i = 1:length(conditions)
            condition = conditions{i};
            step_data = [];
            % Manually access data for the current step across trials
            for j = 1:length(condition_data_avg_blank.(condition).steps_both_median)
                step_data = [step_data; condition_data_avg_blank.(condition).steps_both_mean{j}.(step_field)];
            end
            plotConditionDataWithShading(step_data, x_values, color_map(condition), condition);
        end
        title(['Mean ± SEM - ', titles{row_idx}]);
        xlabel('Time (s)');
        ylabel(y_labels{row_idx});
        legend(conditions, 'Location', 'northeast');
        hold off;

        % Column 2: Individual Trials plot for the processing step
        subplot(length(titles), 2, (row_idx - 1) * 2 + 2);
        hold on;
        for i = 1:length(conditions)
            condition = conditions{i};
            step_data = [];
            % Manually access data for the current step across trials
            for j = 1:length(condition_data_avg_blank.(condition).steps_both_median)
                step_data = [step_data; condition_data_avg_blank.(condition).steps_both_median{j}.(step_field)];
            end
            plotIndividualTrials(step_data, x_values, color_map(condition), condition);
        end
        title(['Individual Trials - ', titles{row_idx}]);
        xlabel('Time (s)');
        ylabel(y_labels{row_idx});
        hold off;
    end

    % Set overall figure title and save
    sgtitle(['Processing Steps Comparison with Individual Trials - ', base]);
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf, fullfile(analysis_folder, sprintf('processing_steps_comparisons_%s.png', base)));
end

function plotConditionDataWithShading(data, x_values, color, condition)
    % Plot mean ± SEM with shading for each condition
    avg = mean(data, 1, 'omitnan');
    std_error = std(data, 0, 1, 'omitnan') / sqrt(size(data, 1));
    nan_mask = ~isnan(avg) & ~isnan(std_error);
    x_shaded = x_values(nan_mask);
    y_upper = avg(nan_mask) + std_error(nan_mask);
    y_lower = avg(nan_mask) - std_error(nan_mask);
    fill([x_shaded, fliplr(x_shaded)], [y_upper, fliplr(y_lower)], color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(x_values, avg, 'Color', color, 'LineWidth', 2, 'DisplayName', condition);
end

function plotIndividualTrials(data, x_values, color, condition)
    % Plot individual trials for each condition as semi-transparent lines
    for trial = 1:size(data, 1)
        plot(x_values, data(trial, :), 'Color', [color, 0.3]);  % Semi-transparent lines for individual trials
    end
    % Mean line for reference
    plot(x_values, mean(data, 1, 'omitnan'), 'Color', color, 'LineWidth', 2, 'DisplayName', [condition, ' Mean']);
    xline(0, '--', 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5]);  % Reference line for target found
end
