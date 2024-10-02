function plotQPieChart(data_numeric, group_name, color_map, number_questions, safe, comparison_results_folder)
% Function to plot a pie chart of response proportions for a group
%
% Inputs:
%   data_numeric: Numeric responses for the group
%   group_name: Name of the group ('ADHD' or 'non-ADHD')
%   color_map: Color mapping for groups

% Count occurrences of each response category
counts_total = histcounts(data_numeric(:), 0.5:1:5.5);

% Plot the pie chart
figure;
pie(counts_total);
legend({'Never', 'Rarely', 'Sometimes', 'Often', 'Very Often'}, 'Location', 'bestoutside');
title([number_questions, ' Response Proportions for ', group_name, ' Group']);

if safe == 1
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf, fullfile(comparison_results_folder, strcat('quest_pie_answer_distr_group_', group_name, '_', number_questions, '.png')));
end
end
