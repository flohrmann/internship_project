function plotRTOverTimeColouredByCondition(data, color_map, comparison_results_folder, condition_labels)
    % Initialize figure with two subplots
    figure;

    %% First subplot: Reaction times over trials
    subplot(2, 1, 1);
    hold on;

    % Loop through each trial and plot reaction times without adding to legend
    for i = 1:size(data, 1)
        % Get the color directly from color_map based on Condition
        condition = char(data.Condition{i});
        plot_color = color_map(condition);

        % Plot each point with a single marker, hiding it from the legend
        plot(i, data.rt(i), 'o', 'Color', plot_color, 'MarkerSize', 6, 'HandleVisibility', 'off');
    end

    % Add legend entries for each unique condition without showing individual points
    condition_keys = keys(color_map);
    for j = 1:length(condition_keys)
        % Get the color from the color map for the legend entry
        condition = condition_keys{j};
        legend_color = color_map(condition);

        % Create a dummy plot for each unique condition for the legend
        plot(NaN, NaN, 'o', 'Color', legend_color, 'DisplayName', condition_labels{j});
    end

    % Finalize the first subplot
    title('Reaction Times Over Trials');
    xlabel('Trial Number');
    ylabel('Reaction Time');
    set(gca, 'YScale', 'log');
    legend('Location', 'northeastoutside');
    hold off;

    %% Second subplot: Violin plot with error bars for each condition
    subplot(2, 1, 2);
    hold on;

    % Prepare data for violin plot
    all_rt = [];
    all_conditions = [];

    % Collect reaction times and conditions for plot
    for j = 1:length(condition_keys)
        % Filter reaction times based on the current condition
        condition = condition_keys{j};
        rt_condition = data.rt(strcmp(data.Condition, condition));  % Use strcmp for cell array comparison

        % Append to arrays for violin plot
        all_rt = [all_rt; rt_condition];
        all_conditions = [all_conditions; repmat({condition}, size(rt_condition))];
    end

    % Create the violin plot
    %errorbar(all_rt, all_conditions, 'ShowData', false, 'ShowMean', true, 'ShowBox', false);

    % Calculate and add error bars for each condition
    for j = 1:length(condition_keys)
        condition = condition_keys{j};
        rt_condition = data.rt(strcmp(data.Condition, condition));
        
        % Calculate mean and standard error for each condition
        avg_rt = median(rt_condition);
        stderr_rt = std(rt_condition) / sqrt(length(rt_condition));  % Standard error of the mean

        % Plot error bars on top of the violin plot
        x_loc = j;  % x-location for error bars
        bar(j,avg_rt, 'FaceColor', color_map(condition_keys{j}), 'EdgeColor', 'none', ...
            'FaceAlpha', 0.3, 'DisplayName', condition_labels{j}, 'HandleVisibility', 'off');
        errorbar(x_loc, avg_rt, stderr_rt, 'k', 'LineWidth', 1.5);
    end

    % Finalize the second subplot
    title('Reaction Time Distribution per Condition with Median and STD');
    xlabel('Condition');
    ylabel('Reaction Time');
    %set(gca, 'YScale', 'log');
    xticks([1,2,3,4]);
    xticklabels(condition_labels);
    hold off;

    % Save the figure
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf, fullfile(comparison_results_folder, 'RTbutton_violin_allTrials.png'));

end
