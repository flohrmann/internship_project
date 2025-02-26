function plotRTDistributionPerCondition(trial_results, analysis_folder, condition_labels, color_map)
    % Extract conditions and reaction times
    conditions = categorical(trial_results.Condition);
    rt = trial_results.rt;
    uniqueConditions = categories(conditions);
    numConditions = length(uniqueConditions);

    % Initialize figure
    figure;
    hold on;

    % Loop through each unique condition
    for i = 1:numConditions
        % Get the reaction times for the current condition
        conditionIndices = conditions == uniqueConditions{i};
        rt_condition = rt(conditionIndices);

        % Get color from color_map based on condition label
        conditionName = char(uniqueConditions{i});
        if isKey(color_map, conditionName)
            color = color_map(conditionName);
        else
            color = [0 0 0]; % Default to black if condition not in color_map
        end

        % Plot histogram with specified color
        histogram(rt_condition, 'BinWidth', 1, 'DisplayName', conditionName, ...
            'FaceColor', color);
    end
    hold off;

    % Label and title adjustments
    xlabel('Reaction Time (s)');
    ylabel('Density');
    title('Reaction Time Distribution per Condition');
    legend(condition_labels);

    % Save the figure
    saveas(gcf, fullfile(analysis_folder, 'rt_distr_condition.png'));
end
