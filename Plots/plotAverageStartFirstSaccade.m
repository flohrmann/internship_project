function first_sacc = plotAverageStartFirstSaccade(id, trial_metrics, trial_results, conditions, ...
                                                     condition_labels, color_map, compare_folder, analysis_folder)

% Extract data
first_saccade_start_time = trial_metrics.first_saccade_start_time;
conditions_data = trial_results.Condition; % Cell array of condition labels

% Initialize data for plotting
condition_colors = cellfun(@(cond) color_map(cond), conditions, 'UniformOutput', false);

% Group data by condition
first_sacc = cell(size(conditions));
for i = 1:numel(conditions)
    first_sacc{i} = first_saccade_start_time(strcmp(conditions_data, conditions{i}));
end


figure;
hold on;
% Plot individual scatter points and boxplots
for i = 1:numel(conditions)
    % Scatter points
    jittered_x = repmat(i, size(first_sacc{i})) + 0.1 * (rand(size(first_sacc{i})) - 0.5);
    scatter(jittered_x, first_sacc{i}, 20, 'filled', 'MarkerFaceColor', condition_colors{i}, 'MarkerEdgeColor', 'none');
    %plot(jittered_x, grouped_data{i}, '.', 'MarkerSize', 20, 'Color', condition_colors{i});
    
    % Boxplot for the condition
    %         boxchart(repmat(i, size(grouped_data{i})), grouped_data{i}, ...
    %             'BoxFaceColor', condition_colors{i}, 'BoxWidth', 0.5);
    % Boxplot for the condition (disable outliers)
    boxchart(repmat(i, size(first_sacc{i})), first_sacc{i}, ...
        'BoxFaceColor', condition_colors{i}, 'BoxWidth', 0.5, 'MarkerStyle', 'none');
end

xticks(1:numel(conditions));
xticklabels(condition_labels);
%xlabel('Condition');
ylabel('First Saccade Start Time (s)');
title(strcat('ID', num2str(id) ,': Start of first Saccade after Search Task Starts'));
grid on; hold off;

saveas(gcf, fullfile(analysis_folder, strcat('first_saccade_rt_by_condition_colored_subj_',num2str(id),'.svg')));
saveas(gcf, fullfile(compare_folder,  strcat('first_saccade_rt_by_condition_colored_subj_',num2str(id),'.svg')));


% get mean per condition
for i = 1:numel(conditions)
    current = first_sacc{1,i};
    first_sacc{2,i} = mean(current, 'omitnan');
end
first_sacc = cell2table(first_sacc, 'VariableNames', {'a', 'a_simple', 'b', 'b_simple'});

save(fullfile(analysis_folder, strcat('\first_sacc.mat')), 'first_sacc');



