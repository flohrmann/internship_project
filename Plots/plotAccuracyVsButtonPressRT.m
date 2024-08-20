function plotAccuracyVsButtonPressRT(data, color_map, comparison_results_folder, safe)
    % Manually set the x-axis labels
    x_labels = {'a', 'a simple', 'b', 'b simple'};

    rt_button_press = [];
    accuracy_all = [];
    condition_all = [];
    group_all = {};

    % Loop through the struct to gather data for RT comparison and accuracy
    for i = 1:length(data)
        rt_button_press = [rt_button_press; data(i).rt];
        accuracy_all = [accuracy_all; data(i).accuracy];
        condition_all = [condition_all; data(i).Condition];
        group_all = [group_all; repmat({data(i).group}, length(data(i).Condition), 1)];
    end

    % Convert group and condition labels to categorical
    condition_all = categorical(condition_all);
    group_all = categorical(group_all);

    % Calculate mean RT and accuracy per condition and group
    [unique_conditions, ~, cond_idx] = unique(condition_all);
    [unique_groups, ~, group_idx] = unique(group_all);

    mean_rt_button_press = zeros(length(unique_conditions), length(unique_groups));
    mean_accuracy = zeros(length(unique_conditions), length(unique_groups));

    for i = 1:length(unique_conditions)
        for j = 1:length(unique_groups)
            idx = cond_idx == i & group_idx == j;
            mean_rt_button_press(i, j) = mean(rt_button_press(idx));
            mean_accuracy(i, j) = mean(accuracy_all(idx));
        end
    end

    % Create figure
    figure;

    % Plot accuracy vs. button press RT for each group and condition
    hold on;
    for j = 1:length(unique_groups)
        group = char(unique_groups(j));
        x_data = mean_rt_button_press(:, j);
        y_data = mean_accuracy(:, j);
        scatter(x_data, y_data, 100, 'MarkerFaceColor', color_map(group), ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        
        % Label each point with the corresponding label from x_labels at the top right with a small offset
        x_offset = 0.003;  % Adjust this value to control the horizontal distance
        y_offset = 0.003;  % Adjust this value to control the vertical distance
        for k = 1:length(x_labels)
            text(x_data(k) + x_offset, y_data(k) + y_offset, x_labels{k}, ...
                'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
                'FontSize', 10, 'Color', color_map(group), 'FontWeight', 'bold');
        end
    end
    hold off;

    % Customize the plot
    xlabel('Button Press RT (s)');
    ylabel('Accuracy');
    title('Accuracy vs. Button Press Reaction Time per Condition/Group');
    legend(unique_groups, 'Location', 'northeastoutside');
    grid on;

    % Save the figure if safe is set to 1
    if safe == 1
        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf, fullfile(comparison_results_folder, 'accuracy_vs_button_press_rt.png'));
    end
end
