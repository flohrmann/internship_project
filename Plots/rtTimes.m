function rt_table = rtTimes(cutData, analysis_folder)
    conditions = unique(cutData.Condition);
    rt_table = cell(length(conditions), 5); % Initialize as a cell array with 5 columns
    all_rts = [];
    all_conditions = [];

    for i = 1:length(conditions)
        condition = conditions{i};
        condition_data = cutData(strcmp(cutData.Condition, condition), :);
        reaction_times_condition = condition_data.rt;

        min_rt = min(reaction_times_condition);
        max_rt = max(reaction_times_condition);
        mean_rt = mean(reaction_times_condition);
        median_rt = median(reaction_times_condition);

        fprintf('Condition: %s\n', condition);
        fprintf('Shortest Reaction Time: %.3f s\n', min_rt);
        fprintf('Longest Reaction Time: %.3f s\n', max_rt);
        fprintf('Mean Reaction Time: %.3f s\n', mean_rt);
        fprintf('Median Reaction Time: %.3f s\n', median_rt);

        rt_table{i, 1} = condition;
        rt_table{i, 2} = min_rt;
        rt_table{i, 3} = max_rt;
        rt_table{i, 4} = mean_rt;
        rt_table{i, 5} = median_rt;

        % Collect all reaction times and conditions for boxplot
        all_rts = [all_rts; reaction_times_condition];
        all_conditions = [all_conditions; repmat({condition}, length(reaction_times_condition), 1)];
    end

    % Convert the results to a table and name the variables
    rt_table = cell2table(rt_table, 'VariableNames', {'Condition', 'MinRT', 'MaxRT', 'MeanRT', 'MedianRT'});

    % Create the boxplot
    figure;
    boxplot(all_rts, all_conditions);
    title('Reaction Times by Condition');
    xlabel('Condition');
    ylabel('Reaction Time (s)');
    saveas(gcf,strcat(analysis_folder, '\rt_per_condition.png'));
end
