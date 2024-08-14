function plotMeanAccuracyPerGroupCondition(data, color_map, safe, comparison_results_folder)
accuracy_all = [];
group_all = {};
condition_all = {};

for i = 1:length(data)
    accuracy_all = [accuracy_all; data(i).accuracy];
    group_all = [group_all; repmat({data(i).group}, length(data(i).accuracy), 1)];
    condition_all = [condition_all; data(i).Condition];
end

% Convert group and condition labels to categorical for easier comparison
group_all = categorical(group_all);
condition_all = categorical(condition_all);

% Calculate mean accuracy per condition and group
[unique_conditions, ~, cond_idx] = unique(condition_all);
[unique_groups, ~, group_idx] = unique(group_all);

mean_accuracy = zeros(length(unique_conditions), length(unique_groups));

for i = 1:length(unique_conditions)
    for j = 1:length(unique_groups)
        idx = cond_idx == i & group_idx == j;
        mean_accuracy(i, j) = mean(accuracy_all(idx));
    end
end

% Create figure
figure;
hold on;

% Plot mean points and connect them with lines for each group
for j = 1:length(unique_groups)
    group = char(unique_groups(j));
    x_data = 1:length(unique_conditions);
    y_data = mean_accuracy(:, j);
    plot(x_data, y_data, '-o', 'Color', color_map(group), 'MarkerFaceColor', color_map(group), ...
        'LineWidth', 2, 'MarkerSize', 8);
end

% Customize the plot
xticks(1:length(unique_conditions));
xticklabels(unique_conditions);
xlabel('Condition');
ylabel('Mean Accuracy');
title('Mean Accuracy Comparison Between ADHD and nonADHD Groups per Condition');
legend(unique_groups, 'Location', 'northeastoutside');

% Save the figure
if safe == 1
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf, fullfile(comparison_results_folder, 'mean_accuracy_group_condition.png'));
else
end
hold off;
end
