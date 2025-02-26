function plotProcessingComparisonsWithTrials(condition_data_avg_blank, conditions, analysis_folder, color_map, sr, base)
    % Create figure with 4 rows of plots, each representing a processing step
    % Each row has two columns for: Mean ± SEM and individual trials
    figure('Name', sprintf('Processing Comparisons with Individual Trials [%s]', base), 'NumberTitle', 'off');
    
    % Define number of time points before and after based on the first condition
    num_before = size(condition_data_avg_blank.a.beforeMean, 2);
    num_after = size(condition_data_avg_blank.a.afterMean, 2);
    x_before = linspace(-sr * (num_before - 1), 0, num_before);
    x_after = linspace(sr, sr * num_after, num_after);
    x_values = [x_before, x_after];

    % Define plot position counters
    plot_row_count = 4;  % Rows for each processing step
    plot_column_count = 2;  % Columns: Mean ± SEM and individual trials

    % 1st Row: Mean ± SEM and individual trials for concatenated beforeMean and afterMean
    plotDataComparison(condition_data_avg_blank, conditions, x_values, 'beforeMean', 'afterMean', ...
                       color_map, 1, plot_row_count, plot_column_count, 'Before & After (Mean)', ...
                       'Pupil Dilation (Mean)', sr, base, analysis_folder);

    % 2nd Row: Mean ± SEM and individual trials for concatenated beforeMedian and afterMedian
    plotDataComparison(condition_data_avg_blank, conditions, x_values, 'beforeMedian', 'afterMedian', ...
                       color_map, 2, plot_row_count, plot_column_count, 'Before & After (Median)', ...
                       'Pupil Dilation (Median)', sr, base, analysis_folder);

    % 3rd Row: processed_both_mean
    plotSingleData(condition_data_avg_blank, conditions, x_values, 'processed_both_mean', color_map, 3, ...
                   plot_row_count, plot_column_count, 'Processed Both (Mean)', 'Pupil Dilation (Mean)', sr, base);

    % 4th Row: processed_both_median
    plotSingleData(condition_data_avg_blank, conditions, x_values, 'processed_both_median', color_map, 4, ...
                   plot_row_count, plot_column_count, 'Processed Both (Median)', 'Pupil Dilation (Median)', sr, base);

    % Set overall figure title and save
    sgtitle(['Processing Steps Comparison with Individual Trials - ', base]);
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf, fullfile(analysis_folder, sprintf('processing_comparisons_trials_%s.png', base)));
end

function plotDataComparison(data, conditions, x_values, field_before, field_after, color_map, ...
                            row_idx, row_count, col_count, plot_title, ylabel_text, sr, base, analysis_folder)
    % Plot average (Mean ± SEM) and individual trials for each condition

    % Mean or Median plot (Column 1)
    subplot(row_count, col_count, (row_idx - 1) * col_count + 1);
    hold on;
    for i = 1:length(conditions)
        condition = conditions{i};
        data_before = data.(condition).(field_before);
        data_after = data.(condition).(field_after);
        concatenated_data = [data_before, data_after];
        plotConditionDataWithShading(concatenated_data, x_values, color_map(condition), condition);
    end
    title(['Mean ± SEM - ', plot_title]);
    xlabel('Time (s)');
    ylabel(ylabel_text);
    legend(conditions, 'Location', 'northeast');
    hold off;

    % Individual Trials plot (Column 2)
    subplot(row_count, col_count, (row_idx - 1) * col_count + 2);
    hold on;
    for i = 1:length(conditions)
        condition = conditions{i};
        data_before = data.(condition).(field_before);
        data_after = data.(condition).(field_after);
        concatenated_data = [data_before, data_after];
        plotIndividualTrials(concatenated_data, x_values, color_map(condition), condition);
    end
    title(['Individual Trials - ', plot_title]);
    xlabel('Time (s)');
    ylabel(ylabel_text);
    hold off;
end

function plotSingleData(data, conditions, x_values, field_name, color_map, row_idx, row_count, col_count, ...
                        plot_title, ylabel_text, sr, base)
    % Plot individual trials and mean for a single data field

    % Mean ± SEM plot (Column 1)
    subplot(row_count, col_count, (row_idx - 1) * col_count + 1);
    hold on;
    for i = 1:length(conditions)
        condition = conditions{i};
        field_data = data.(condition).(field_name);
        %x_values = linspace(-sr * (size(field_data, 2) - 1) / 2, sr * (size(field_data, 2) / 2), size(field_data, 2));
        plotConditionDataWithShading(field_data, x_values, color_map(condition), condition);
    end
    title(['Mean ± SEM - ', plot_title]);
    xlabel('Time (s)');
    ylabel(ylabel_text);
    legend(conditions, 'Location', 'northeast');
    hold off;

    % Individual Trials plot (Column 2)
    subplot(row_count, col_count, (row_idx - 1) * col_count + 2);
    hold on;
    for i = 1:length(conditions)
        condition = conditions{i};
        field_data = data.(condition).(field_name);
        plotIndividualTrials(field_data, x_values, color_map(condition), condition);
    end
    title(['Individual Trials - ', plot_title]);
    xlabel('Time (s)');
    ylabel(ylabel_text);
    hold off;
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
        plot(x_values, data(trial, :), 'Color', [color, 0.3]);  % Make lines semi-transparent
    end
    plot(x_values, mean(data, 1, 'omitnan'), 'Color', color, 'LineWidth', 2, 'DisplayName', [condition, ' Mean']);  % Add mean line for visibility
    xline(0, '--', 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5]);  % Reference line for target found
end
