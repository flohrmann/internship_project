function plotMeanRTButtonPressVsEyecomparison(data, color_map, comparison_results_folder, safe)
    % Manually set the x-axis labels
    x_labels = {'a', 'a simple', 'b', 'b simple'};

    rt_button_press = [];
    norm_rt_button_press = [];
    rt_eye = [];
    norm_rt_eye = [];
    condition_all = [];
    group_all = {};

    % Loop through the struct to gather data for RT comparison
    for i = 1:length(data)
        % Filter out trials where rt_eye is NaN
        valid_trials = ~isnan(data(i).rt_eye);
        
        % Only use data from valid trials
        rt_button_press = [rt_button_press; data(i).rt(valid_trials)];
        norm_rt_button_press = [norm_rt_button_press; data(i).normalized_rt(valid_trials)];
        
        rt_eye = [rt_eye; data(i).rt_eye(valid_trials)];
        norm_rt_eye = [norm_rt_eye; data(i).normalized_rt_eye(valid_trials)];

        condition_all = [condition_all; data(i).Condition(valid_trials)];
        group_all = [group_all; repmat({data(i).group}, sum(valid_trials), 1)];
    end

    % Convert group and condition labels to categorical
    condition_all = categorical(condition_all);
    group_all = categorical(group_all);

    % Calculate mean RT per condition and group
    [unique_conditions, ~, cond_idx] = unique(condition_all);
    [unique_groups, ~, group_idx] = unique(group_all);

    mean_rt_button_press = zeros(length(unique_conditions), length(unique_groups));
    mean_norm_rt_button_press = zeros(length(unique_conditions), length(unique_groups));
    mean_rt_eye = zeros(length(unique_conditions), length(unique_groups));
    mean_norm_rt_eye = zeros(length(unique_conditions), length(unique_groups));
    mean_rt_diff = zeros(length(unique_conditions), length(unique_groups));
    mean_norm_rt_diff = zeros(length(unique_conditions), length(unique_groups));

    for i = 1:length(unique_conditions)
        for j = 1:length(unique_groups)
            idx = cond_idx == i & group_idx == j;
            mean_rt_button_press(i, j) = mean(rt_button_press(idx));  % Use mean since NaNs are removed
            mean_norm_rt_button_press(i, j) = mean(norm_rt_button_press(idx));  % Use mean since NaNs are removed
            
            mean_rt_eye(i, j) = mean(rt_eye(idx));  % Use mean since NaNs are removed
            mean_norm_rt_eye(i, j) = mean(norm_rt_eye(idx));  % Use mean since NaNs are removed
            
            mean_rt_diff(i, j) = mean_rt_button_press(i, j) - mean_rt_eye(i, j);
            mean_norm_rt_diff(i, j) = mean_norm_rt_button_press(i, j) - mean_norm_rt_eye(i, j);
        end
    end

    % Create figure with 3x2 subplots
    figure;

    % Subplot 1: Eye movement RTs
    subplot(3, 2, 1);
    hold on;
    for j = 1:length(unique_groups)
        group = char(unique_groups(j));
        x_data = 1:length(unique_conditions);
        y_data_eye = mean_rt_eye(:, j);
        plot(x_data, y_data_eye, '-o', 'Color', color_map(group), ...
            'MarkerFaceColor', color_map(group), 'LineWidth', 2, 'MarkerSize', 8);
    end
    hold off;
    xticks(1:length(x_labels));
    xticklabels(x_labels);
    ylabel('Eye RT (s)');
    title('Eye Movement Reaction Time');
    legend(unique_groups, 'Orientation', 'horizontal');
    grid on;

    % Subplot 2: Normalized Eye movement RTs
    subplot(3, 2, 2);
    hold on;
    for j = 1:length(unique_groups)
        group = char(unique_groups(j));
        x_data = 1:length(unique_conditions);
        y_data_norm_eye = mean_norm_rt_eye(:, j);
        plot(x_data, y_data_norm_eye, '-o', 'Color', color_map(group), ...
            'MarkerFaceColor', color_map(group), 'LineWidth', 2, 'MarkerSize', 8);
    end
    hold off;
    xticks(1:length(x_labels));
    xticklabels(x_labels);
    ylabel('Normalized Eye RT (s)');
    title('Normalized Eye Movement Reaction Time');
    grid on;

    % Subplot 3: Button press RTs
    subplot(3, 2, 3);
    hold on;
    for j = 1:length(unique_groups)
        group = char(unique_groups(j));
        x_data = 1:length(unique_conditions);
        y_data_button_press = mean_rt_button_press(:, j);
        plot(x_data, y_data_button_press, '-s', 'Color', color_map(group), ...
            'MarkerFaceColor', color_map(group), 'LineWidth', 2, 'MarkerSize', 8);
    end
    hold off;
    xticks(1:length(x_labels));
    xticklabels(x_labels);
    ylabel('Button Press RT (s)');
    title('Button Press Reaction Time');
    grid on;

    % Subplot 4: Normalized Button press RTs
    subplot(3, 2, 4);
    hold on;
    for j = 1:length(unique_groups)
        group = char(unique_groups(j));
        x_data = 1:length(unique_conditions);
        y_data_norm_button_press = mean_norm_rt_button_press(:, j);
        plot(x_data, y_data_norm_button_press, '-s', 'Color', color_map(group), ...
            'MarkerFaceColor', color_map(group), 'LineWidth', 2, 'MarkerSize', 8);
    end
    hold off;
    xticks(1:length(x_labels));
    xticklabels(x_labels);
    ylabel('Normalized Button Press RT (s)');
    title('Normalized Button Press Reaction Time');
    grid on;

    % Subplot 5: Difference between button press RT and eye RT
    subplot(3, 2, 5);
    hold on;
    for j = 1:length(unique_groups)
        group = char(unique_groups(j));
        x_data = 1:length(unique_conditions);
        y_data_diff = mean_rt_diff(:, j);
        plot(x_data, y_data_diff, '-^', 'Color', color_map(group), ...
            'MarkerFaceColor', color_map(group), 'LineWidth', 2, 'MarkerSize', 8);
    end
    hold off;
    xticks(1:length(x_labels));
    xticklabels(x_labels);
    xlabel('Condition');
    ylabel('Button Press RT - Eye RT (s)');
    title('Difference Between Button Press RT and Eye RT');
    grid on;

    % Subplot 6: Difference between normalized button press RT and normalized eye RT
    subplot(3, 2, 6);
    hold on;
    for j = 1:length(unique_groups)
        group = char(unique_groups(j));
        x_data = 1:length(unique_conditions);
        y_data_norm_diff = mean_norm_rt_diff(:, j);
        plot(x_data, y_data_norm_diff, '-^', 'Color', color_map(group), ...
            'MarkerFaceColor', color_map(group), 'LineWidth', 2, 'MarkerSize', 8);
    end
    hold off;
    xticks(1:length(x_labels));
    xticklabels(x_labels);
    xlabel('Condition');
    ylabel('Norm Button RT - Norm Eye RT (s)');
    title('Difference Normalized Button Press RT and Normalized Eye RT');
    grid on;

    % Save the figure if safe is set to 1
    if safe == 1
        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf, fullfile(comparison_results_folder, 'mean_rt_comparison_per_condition_group.png'));
    end
end
