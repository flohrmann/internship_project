function plotBarAccuracyGroupCondition(data, color_map, comparison_results_folder, safe)

accuracy_all = [];
group_all = {};
condition_all = {};

for i = 1:length(data)
    accuracy_all = [accuracy_all; data(i).accuracy];
    group_all = [group_all; repmat({data(i).group}, length(data(i).accuracy), 1)];
    condition_all = [condition_all; data(i).Condition];
end

group_all = categorical(group_all);
condition_all = categorical(condition_all);

[unique_conditions, ~, cond_idx] = unique(condition_all);
[unique_groups, ~, group_idx] = unique(group_all);

mean_accuracy = zeros(length(unique_conditions), length(unique_groups));

for i = 1:length(unique_conditions)
    for j = 1:length(unique_groups)
        idx = cond_idx == i & group_idx == j;
        mean_accuracy(i, j) = mean(accuracy_all(idx));
    end
end

figure;
b = bar(categorical(unique_conditions), mean_accuracy);

for k = 1:length(b)
    b(k).FaceColor = color_map(char(unique_groups(k)));
    hold on;
    % Calculate the x-data positions for the mean points
    x_data = b(k).XEndPoints;
    y_data = mean_accuracy(:, k);
    % Convert x_data to categorical for plotting with plot()
    x_data_categorical = categorical(unique_conditions);
    plot(x_data_categorical, y_data, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);  % Add mean points
end

xlabel('Condition');
ylabel('Accuracy');
title('Accuracy Comparison Between ADHD and non-ADHD Groups per Condition');
legend(unique_groups, 'Location', 'northeastoutside');
if safe == 1
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf, fullfile(comparison_results_folder, 'accuracy_comparison_per_condition.png'));
end