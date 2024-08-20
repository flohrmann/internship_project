function plotBarSEMMeanRTConditonGroup(groups, conditions, data_struct_norm_mean, color_map, safe, comparison_results_folder)

% Initialize arrays for means and SEMs
means = zeros(length(groups), length(conditions));
sems = zeros(length(groups), length(conditions)); 

% Loop through groups and conditions to retrieve normalized means and SEMs
for g = 1:length(groups)
    for c = 1:length(conditions)
        % Initialize arrays to store RTs and SEMs for each participant in the group
        rt_values = [];
        sem_values = [];
        
        % Loop through each participant in the data_struct_norm_mean
        for i = 1:length(data_struct_norm_mean)
            if strcmp(data_struct_norm_mean(i).group, groups{g})
                % Retrieve normalized RT and SEM for the current condition
                switch conditions{c}
                    case 'a'
                        rt_values = [rt_values; data_struct_norm_mean(i).nRTa];
                        sem_values = [sem_values; data_struct_norm_mean(i).nSEMa];
                    case 'b'
                        rt_values = [rt_values; data_struct_norm_mean(i).nRTb];
                        sem_values = [sem_values; data_struct_norm_mean(i).nSEMb];
                    case 'a_simple'
                        rt_values = [rt_values; data_struct_norm_mean(i).nRTasimple];
                        sem_values = [sem_values; data_struct_norm_mean(i).nSEMasimple];
                    case 'b_simple'
                        rt_values = [rt_values; data_struct_norm_mean(i).nRTbsimple];
                        sem_values = [sem_values; data_struct_norm_mean(i).nSEMbsimple];
                end
            end
        end
        
        % Calculate the mean and SEM for the group and condition
        means(g, c) = mean(rt_values);
        sems(g, c) = mean(sem_values); % Use mean of SEMs across participants
    end
end

%% Bar plot with error bars
figure;
hold on;
bar_width = 0.4;

for g = 1:length(groups)
    % Offset x-values for groups to prevent overlap
    x_values = (1:length(conditions)) + (g-1) * bar_width;
    y_values = means(g, :);
    err_values = sems(g, :);
    
    % Plot bars with the respective color
    bar(x_values, y_values, bar_width, 'FaceColor', color_map(groups{g}));
    % Plot error bars
    errorbar(x_values, y_values, err_values, 'k', 'LineStyle', 'none');
end

% Adjust x-axis and labels
xticks(1:length(conditions) + bar_width * (length(groups) - 1) / 2);
xticklabels(conditions);
xlabel('Condition');
ylabel('Normalized Mean Reaction Time');
legend(groups, 'Location', 'northeastoutside');
title('Normalized Mean Reaction Time with Standard Error');
hold off;

if safe == 1
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf, fullfile(comparison_results_folder, 'bar_norm_mean_SEM_RTbutton_group_condition.png'));
end

%% Lines mean RT per condition and group with SEM error bars
figure;
hold on;

for g = 1:length(groups)
    x_values = 1:length(conditions);
    y_values = means(g, :);
    err_values = sems(g, :);
    
    % Plot lines with error bars
    errorbar(x_values, y_values, err_values, '-o', 'Color', color_map(groups{g}), ...
        'DisplayName', groups{g}, 'LineWidth', 2, 'MarkerSize', 8);
end

% Customize the plot
xticks(1:length(conditions));
xticklabels(conditions);
xlabel('Condition');
ylabel('Normalized Mean Reaction Time');
legend('Location', 'northeastoutside');
title('Interaction Plot of Normalized Reaction Times by Condition and Group');
grid on;
hold off;

if safe == 1
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf, fullfile(comparison_results_folder, 'lines_norm_mean_SEM_RTbutton_group_condition.png'));
end

end
