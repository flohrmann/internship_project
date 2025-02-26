function plotPupilMagnitudesTwoGroups(avg_adhd, avg_nonadhd, conditions, color_map, comparison_results_folder)
    % Define the groups for the plot
    groups = {'ADHD', 'Non-ADHD'};

    % Conditions and pupil measures
    condition_names = {'a', 'as', 'b', 'bs'};
    xlabels = {'a', 'a simple', 'b', 'b simple'};
    measures = {'min', 'max', 'mean'};  % min, max, mean pupil diameter

    % Initialize figure with 2 rows (groups) and 3 columns (measures)
    figure;
    t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');  % 2 rows, 3 columns layout

    % Group data for easier access
    group_data = {avg_adhd, avg_nonadhd};

    % Initialize separate arrays for axes handles for linking later
    ax_min = [];
    ax_max = [];
    ax_mean = [];

    % Loop through each group (ADHD, Non-ADHD)
    for group_idx = 1:length(groups)
        group_table = group_data{group_idx};  % Get the table for the current group
        num_subjects = size(group_table, 1);

        % Loop through each measure (min, max, mean pupil diameter)
        for measure_idx = 1:length(measures)
            measure = measures{measure_idx};

            % Create a subplot for each group-measure combination
            ax = nexttile;  
            hold on;

            % Initialize arrays to store data for correlation calculation
            measure_data_across_conditions = [];

            % Loop through each condition ('a', 'a_simple', etc.)
            for condition_idx = 1:length(conditions)
                condition = condition_names{condition_idx};

                % Extract the relevant data for this condition and measure
                measure_field = [condition, '_', measure];  % e.g., 'a_min', 'a_max', 'a_mean'
                measure_data = group_table.(measure_field);  % Column for the current condition and measure

                % Stack data for each subject
                measure_data_across_conditions(:, condition_idx) = measure_data;  % Columns for each condition

                % Scatter plot the averages of individual subjects for this condition
                scatter(repmat(condition_idx, num_subjects, 1), measure_data, 50, ...
                    color_map(conditions{condition_idx}), 'filled');

                % Compute group mean and 95% confidence interval
                mean_values = mean(measure_data, 'omitnan');
                ci_95 = 1.96 * std(measure_data, 0, 'omitnan') / sqrt(num_subjects);

                % Plot the group mean and confidence interval as horizontal lines
                plot([condition_idx - 0.2, condition_idx + 0.2], [mean_values, mean_values], ...
                    'Color', color_map(conditions{condition_idx}), 'LineWidth', 2);  % Mean line
                plot([condition_idx - 0.2, condition_idx + 0.2], [mean_values + ci_95, mean_values + ci_95], ...
                    '--', 'Color', color_map(conditions{condition_idx}));  % Upper CI
                plot([condition_idx - 0.2, condition_idx + 0.2], [mean_values - ci_95, mean_values - ci_95], ...
                    '--', 'Color', color_map(conditions{condition_idx}));  % Lower CI
            end

%             % Calculate the Spearman correlation between the conditions
%             rho_values = [];
%             p_values = [];
%             for subj_idx = 1:num_subjects
%                 [rho, p_value] = corr((1:length(conditions))', measure_data_across_conditions(subj_idx, :)', ...
%                                       'Type', 'Spearman', 'Rows', 'complete');
%                 rho_values = [rho_values; rho];
%                 p_values = [p_values; p_value];
%             end
%             % Average correlation values across subjects
%             avg_rho = mean(rho_values, 'omitnan');
%             avg_p = mean(p_values, 'omitnan');

            % X-axis settings
            xticks(1:length(conditions));
            xticklabels(xlabels);

            % Title and labels with Spearman rho and p-value
            ylabel([measure, ' Pupil Diameter']);
%             title([groups{group_idx}, ' - ', measure, sprintf(' (rho = %.2f, p = %.2f)', avg_rho, avg_p)]);
            title([groups{group_idx}, ' - ', measure]);

            % Store axes handle for linking later based on measure type
            if strcmp(measure, 'min')
                ax_min = [ax_min, ax];  % Store axes handle for min
            elseif strcmp(measure, 'max')
                ax_max = [ax_max, ax];  % Store axes handle for max
            elseif strcmp(measure, 'mean')
                ax_mean = [ax_mean, ax];  % Store axes handle for mean
            end

            % Grid on for clarity
            grid on;
            hold off;
        end
    end

    % Link the y-axes separately for min, max, and mean measures
    linkaxes(ax_min, 'y');  % Link axes for "min"
    linkaxes(ax_max, 'y');  % Link axes for "max"
    linkaxes(ax_mean, 'y');  % Link axes for "mean"

    % Add a general title to the entire figure
    title(t, 'Pupil Diameter Magnitudes by Group and Measure with Spearman Correlations', 'FontWeight', 'bold');

    % Save the figure
    saveas(gcf, fullfile(comparison_results_folder, 'pupil_magnitudes_comparison_twogroups_table.png'));
end
