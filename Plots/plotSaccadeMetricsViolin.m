function plotSaccadeMetricsViolin(trial_metrics, cutData, conditions, color_map, condition_labels, safe_name)
figure;
t = tiledlayout(2, 2, 'TileSpacing', 'Compact');

%% 1. Average number of saccades until target found per condition
data_combined = [];
group_combined = [];
for i = 1:length(conditions)
    condition = conditions{i};
    condition_idx = strcmp(cutData.Condition, condition);
    % Convert cells to numeric and compute sum
    n_sac_numeric = trial_metrics.n_sac(condition_idx,1);
    saccades_after_numeric = cellfun(@(x) x(1), trial_metrics.saccades_after_target(condition_idx));
    % Sum n_sac and saccades_after_target
    condition_data = n_sac_numeric + saccades_after_numeric;
    data_combined = [data_combined; condition_data];
    group_combined = [group_combined; repmat({condition}, length(condition_data), 1)];
end
valid_idx = ~isnan(data_combined);
data_combined = data_combined(valid_idx);
group_combined = group_combined(valid_idx);
nexttile;
vp = violinplot(data_combined, group_combined);
% Apply consistent colors from color_map
for i = 1:length(vp)
    condition = conditions{i};  % Match the condition by order
    vp(1,i).ViolinColor{1,1} = color_map(condition);
    %vp(1,i).ShowData = false;  % Disable scatter points
end
title('Saccades Until Target Found');
ylabel('Number of Saccades'); set(gca, 'YScale', 'log');
%% 2. Average number of saccades after target found per condition
data_combined = [];
group_combined = [];
for i = 1:length(conditions)
    condition = conditions{i};
    condition_idx = strcmp(cutData.Condition, condition);
    condition_data = cellfun(@(x) x(1), trial_metrics.saccades_after_target(condition_idx));
    data_combined = [data_combined; condition_data];
    group_combined = [group_combined; repmat({condition}, length(condition_data), 1)];
end
valid_idx = ~isnan(data_combined);
data_combined = data_combined(valid_idx);
group_combined = group_combined(valid_idx);
nexttile;
vp = violinplot(data_combined, group_combined);
for i = 1:length(vp)
    condition = conditions{i};  % Match the condition by order
    vp(1,i).ViolinColor{1,1} = color_map(condition);
    %vp(1,i).ShowData = false;  % Disable scatter points
end
title('Saccades After Target Found');
ylabel('Number of Saccades'); set(gca, 'YScale', 'log');
%% 3. Average distance of saccade before and towards target per condition
%     data_combined = [];
%     group_combined = [];
%     for i = 1:length(conditions)
%         condition = conditions{i};
%         condition_idx = strcmp(cutData.Condition, condition);
%         % Convert cells to numeric
%         before_data = cell2mat(trial_metrics.saccade_distance_before_target(condition_idx));
%         towards_data = cell2mat(trial_metrics.saccade_distance_towards_target(condition_idx));
%
%         % Append to combined data
%         data_combined = [data_combined; before_data; towards_data];
%         group_combined = [group_combined; ...
%                           repmat({[condition, ' (Before)']}, length(before_data), 1); ...
%                           repmat({[condition, ' (Towards)']}, length(towards_data), 1)];
%    end
%
%     % Violin plot for distance before and towards target
%     nexttile;
%     violinplot(data_combined, group_combined);
%     title('Saccade Distances');
%     ylabel('Distance (mm)');
%     legend('Before Target', 'Towards Target', 'Location', 'Best');

%% 3. Number of saccades (n_sacs) per condition
% data_combined = [];
% group_combined = [];
% for i = 1:length(conditions)
%     condition = conditions{i};
%     condition_idx = strcmp(cutData.Condition, condition);
%     n_sacs_data = trial_metrics.n_sac(condition_idx);
%     data_combined = [data_combined; n_sacs_data];
%     group_combined = [group_combined; repmat({condition}, length(n_sacs_data), 1)];
% end
%
% % Plot the violin plot
% nexttile;
% violinplot(data_combined, group_combined);
% title('Number of Saccades (n_sacs)');
% ylabel('Number of Saccades');set(gca, 'YScale', 'log');
%% 3. Distance of saccades per condition
data_combined = [];
group_combined = [];
for i = 1:length(conditions)
    condition = conditions{i};
    condition_idx = strcmp(cutData.Condition, condition);
    saccade_distances_data = trial_metrics.mean_dist_sac(condition_idx);
    data_combined = [data_combined; saccade_distances_data];
    group_combined = [group_combined; repmat({condition}, length(saccade_distances_data), 1)];
end
valid_idx = ~isnan(data_combined);
data_combined = data_combined(valid_idx);
group_combined = group_combined(valid_idx);
% Plot the violin plot
nexttile;
vp = violinplot(data_combined, group_combined);
for i = 1:length(vp)
    condition = conditions{i};  % Match the condition by order
    vp(1,i).ViolinColor{1,1} = color_map(condition);
    %vp(1,i).ShowData = false;  % Disable scatter points
end
title('Mean Saccade Distance per Trial');
ylabel('Distance (mm)');set(gca, 'YScale', 'log');


%% 4. Average difference between optimal angle and first saccade angle
data_combined = [];
group_combined = [];
for i = 1:length(conditions)
    condition = conditions{i};
    condition_idx = strcmp(cutData.Condition, condition);
    % Convert angles to numeric and compute difference
    optimal_angles = rad2deg(trial_metrics.optimal_angle(condition_idx));
    first_saccade_angles = rad2deg(trial_metrics.first_saccade_angle(condition_idx));
    angle_diff = abs(optimal_angles - first_saccade_angles);
    
    % Append to combined data
    data_combined = [data_combined; angle_diff];
    group_combined = [group_combined; repmat({condition}, length(angle_diff), 1)];
end
valid_idx = ~isnan(data_combined);
data_combined = data_combined(valid_idx);
group_combined = group_combined(valid_idx);
nexttile;
vp = violinplot(data_combined, group_combined);
for i = 1:length(vp)
    condition = conditions{i};  % Match the condition by order
    vp(1,i).ViolinColor{1,1} = color_map(condition);
    %vp(1,i).ShowData = false;  % Disable scatter points
end
title('Diff: Optimal and First Saccade Angle');
ylabel('Angle Difference (degrees)');set(gca, 'YScale', 'log');


print(gcf, safe_name, '-dsvg');
end
