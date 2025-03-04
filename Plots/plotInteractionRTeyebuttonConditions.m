function plotInteractionRTeyebuttonConditions(groups, conditions, data, color_map, safe, comparison_results_folder)


rt_button_press_conditions = cell(length(groups), length(conditions));
rt_eye_conditions = cell(length(groups), length(conditions));
accuracy_conditions = cell(length(groups), length(conditions));

for i = 1:length(data)
    group_idx = find(strcmp(groups, data(i).group));
    if isempty(group_idx)
        warning('Group label "%s" not found in the predefined groups list.', data(i).group);
        continue;
    end
    for t = 1:length(data(i).Condition)
        condition = data(i).Condition{t};
        condition_idx = find(strcmp(conditions, condition));
        if isempty(condition_idx)
            warning('Condition label "%s" not found in the predefined conditions list.', condition);
            continue;
        end
        rt_button_press_conditions{group_idx, condition_idx} = [rt_button_press_conditions{group_idx, condition_idx}; data(i).rt(t)];
        rt_eye_conditions{group_idx, condition_idx} = [rt_eye_conditions{group_idx, condition_idx}; data(i).rt_eye(t)];
        accuracy_conditions{group_idx, condition_idx} = [accuracy_conditions{group_idx, condition_idx}; data(i).accuracy(t)];
    end
end

% Visualize Reaction Times Across Conditions and Groups
% todo change condition names
% First Figure: Each group and condition in separate subplots
plotRTScatterGroupCondition(groups, conditions, rt_eye_conditions, rt_button_press_conditions, color_map, safe, comparison_results_folder);
plotRTscatterConditionAcrossGroups(groups, conditions, rt_eye_conditions, rt_button_press_conditions, color_map, safe, comparison_results_folder)
