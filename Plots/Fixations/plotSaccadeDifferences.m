function plotSaccadeDifferences(saccadeStats, data, group_labels, ids, unique_conditions, condition_labels, color_map, color_map_individual, comparison_results_folder, safe)
numConditions = length(unique_conditions);

% Initialize matrices to store group data
adhdAmplitudeData = zeros(numConditions, 0);
nonAdhdAmplitudeData = zeros(numConditions, 0);
adhdVelocityData = zeros(numConditions, 0);
nonAdhdVelocityData = zeros(numConditions, 0);
adhd_subject_ids = {};
nonadhd_subject_ids = {};
adhdRTdistData = zeros(numConditions, 0);
nonadhdRTdistData = zeros(numConditions, 0);

% Loop over participants to group data by ADHD status
for participant = 1:length(saccadeStats)
    conditionAmplitudes = zeros(numConditions, 1);
    conditionVelocities = zeros(numConditions, 1);
    conditionRTdistance = zeros(numConditions, 1);
    
    for c = 1:numConditions
        conditionAmplitudes(c) = saccadeStats(participant).conditions(c).medianSaccadeAmplitude;
        conditionVelocities(c) = saccadeStats(participant).conditions(c).medianSaccadeVelocity;
        conditionRTdistance(c) = saccadeStats(participant).conditions(c).medianDistancePerSecond;
    end
    
    if strcmp(group_labels{participant}, 'ADHD')
        adhdAmplitudeData = [adhdAmplitudeData, conditionAmplitudes];
        adhdVelocityData = [adhdVelocityData, conditionVelocities];
        adhdRTdistData   = [adhdRTdistData, conditionRTdistance];
        adhd_subject_ids = [adhd_subject_ids; saccadeStats(participant).id]; % Store subject IDs
    else
        nonAdhdAmplitudeData = [nonAdhdAmplitudeData, conditionAmplitudes];
        nonAdhdVelocityData = [nonAdhdVelocityData, conditionVelocities];
        nonadhdRTdistData   = [nonadhdRTdistData, conditionRTdistance];
        nonadhd_subject_ids = [nonadhd_subject_ids; saccadeStats(participant).id]; % Store subject IDs
    end
end

%% amplitude
adhdAmplitudeMedian = median(adhdAmplitudeData', 1, 'omitnan');
adhdAmplitudeSem = std(adhdAmplitudeData', 0, 1, 'omitnan')./ sqrt(size(adhdAmplitudeData', 1));
nonAdhdAmplitudeMedian = median(nonAdhdAmplitudeData', 1, 'omitnan');
nonAdhdAmplitudeSem = std(nonAdhdAmplitudeData', 0, 1, 'omitnan')./ sqrt(size(nonAdhdAmplitudeData', 1));

plotADHDnonADHDandDiff('Average Search Path Distance',...
    adhdAmplitudeData', adhdAmplitudeMedian, adhdAmplitudeSem, 'ADHD', 'northwest', ...
    nonAdhdAmplitudeData', nonAdhdAmplitudeMedian, nonAdhdAmplitudeSem, 'nonADHD', 'northeast', ...
    ids, condition_labels, 'Gaze Distance (pixels)', group_labels, unique_conditions, color_map, color_map_individual, ...
    fullfile(comparison_results_folder, '02_saccade_amplitude_allinone_median.png'));

%% velocity
adhdVelocityMedian = median(adhdVelocityData', 1, 'omitnan');
adhdVelocitySem = std(adhdVelocityData', 0, 1, 'omitnan')./ sqrt(size(adhdVelocityData', 1));
nonAdhdVelocityMedian = median(nonAdhdVelocityData', 1, 'omitnan');
nonAdhdVelocitySem = std(nonAdhdVelocityData', 0, 1, 'omitnan')./ sqrt(size(nonAdhdVelocityData', 1));

plotADHDnonADHDandDiff('Average Saccade Velocity',...
    adhdVelocityData', adhdVelocityMedian, adhdVelocitySem, 'ADHD', 'northwest', ...
    nonAdhdVelocityData', nonAdhdVelocityMedian, nonAdhdVelocitySem, 'nonADHD', 'northwest', ...
    ids, condition_labels, 'Velocity (pixels/s)', group_labels, unique_conditions, color_map, color_map_individual, ...
    fullfile(comparison_results_folder, '02_saccade_velocity_allinone_median.png'));

%% amplitude/distance per second
adhdDistMedian = median(adhdRTdistData', 1, 'omitnan');
adhdDistSem = std(adhdRTdistData', 0, 1, 'omitnan')./ sqrt(size(adhdRTdistData', 1));
nonAdhdDistMedian = median(nonadhdRTdistData', 1, 'omitnan');
nonAdhdDistSem = std(nonadhdRTdistData', 0, 1, 'omitnan')./ sqrt(size(nonadhdRTdistData', 1));
plotADHDnonADHDandDiff('Distance Gaze Travels per Second',...
    adhdRTdistData', adhdDistMedian, adhdDistSem, 'ADHD', 'northwest', ...
    nonadhdRTdistData', nonAdhdDistMedian, nonAdhdDistSem, 'nonADHD', 'northeast', ...
    ids, condition_labels, 'Median Distance (pixels/s)', group_labels, unique_conditions, color_map, color_map_individual, ...
    fullfile(comparison_results_folder, '02_distance_per_second_allinone_median.png'));

%% start first saccade




%% amplitude first saccade



end



% %% Top Left: Scatter + Bar Plot (Amplitude)
% figure;
% subplot(1, 2, 1);
% hold on;
% 
% % Generate color shades for ADHD and non-ADHD subjects
% num_adhd_subjects = length(adhd_subject_ids);
% num_nonadhd_subjects = length(nonadhd_subject_ids);
% adhd_colors = [linspace(0.3, 0, num_adhd_subjects)', linspace(0.8, 0.4, num_adhd_subjects)', linspace(0.3, 0.2, num_adhd_subjects)'];
% nonadhd_colors = [linspace(0.8, 0.4, num_nonadhd_subjects)', linspace(0.3, 0, num_nonadhd_subjects)', linspace(0.3, 0.2, num_nonadhd_subjects)'];
% 
% % Scatter individual data points for ADHD
% adhd_indices = strcmp(group_labels, 'ADHD');
% adhd_ids = ids(adhd_indices);
% 
% adhd_handles = gobjects(num_adhd_subjects, 1); % Store handles for the legend
% for s = 1:num_adhd_subjects
%     adhd_id = adhd_ids(s);
%     adhd_handles(s) = scatter(1:numConditions, adhdAmplitudeData(:, s), 30, ...
%         'MarkerEdgeColor', color_map_individual(adhd_id).color, 'MarkerFaceColor', color_map_individual(adhd_id).color, ...
%         'Marker', color_map_individual(adhd_id).marker,...% 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.3,
%         'DisplayName', strcat('id',num2str(adhd_subject_ids{s}))); % Use subject ID
% end
% 
% % Scatter individual data points for non-ADHD
% nonadhd_indices = strcmp(group_labels, 'nonADHD');
% nonadhd_ids = ids(nonadhd_indices);
% nonadhd_handles = gobjects(num_nonadhd_subjects, 1); % Store handles for the legend
% for s = 1:num_nonadhd_subjects
%     nonadhd_id = nonadhd_ids(s);
%     nonadhd_handles(s) = scatter(1:numConditions, nonAdhdAmplitudeData(:, s), 30, ...
%         'MarkerEdgeColor', color_map_individual(nonadhd_id).color, 'MarkerFaceColor', color_map_individual(nonadhd_id).color, ...
%         'Marker', color_map_individual(nonadhd_id).marker,... % 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.3,
%         'DisplayName', strcat('id',num2str(nonadhd_subject_ids{s}))); % Use subject ID
% end
% 
% % Error bars for ADHD
% adhd_errorbar_handle = errorbar(1:numConditions, adhdAmplitudeMedian, adhdAmplitudeSem, '-o', ...
%     'Color', color_map('ADHD'), 'MarkerFaceColor', color_map('ADHD'), ...
%     'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'ADHD');
% 
% % Error bars for non-ADHD
% nonadhd_errorbar_handle = errorbar(1:numConditions, nonAdhdAmplitudeMedian, nonAdhdAmplitudeSem, '-o', ...
%     'Color', color_map('nonADHD'), 'MarkerFaceColor', color_map('nonADHD'), ...
%     'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'nonADHD');
% 
% 
% % Customize the plot
% xticks(1:numConditions);
% xticklabels(condition_labels);
% %xlabel('Condition');
% ylabel('Average Saccade Amplitude (pixels)');
% %legend([adhd_errorbar_handle, nonadhd_errorbar_handle], 'Location', 'northeast'); % Legend for error bars only
% %legend([adhd_handles; nonadhd_handles], 'Location', 'northeast'); % Add only scatter handles to the legend
% title('Amplitude: Scatter + Error Bar');
% grid off; hold off;
% %% Top Right: Scatter + Bar Plot (Velocity)
% subplot(1, 2, 2); % hold on;
% % Scatter individual data points for ADHD
% adhd_handles = gobjects(num_adhd_subjects, 1); % Store handles for the legend
% for s = 1:num_adhd_subjects
%     adhd_id = adhd_ids(s);
%     adhd_handles(s) = scatter(1:numConditions, adhdVelocityData(:, s), 30, ...
%         'MarkerEdgeColor', color_map_individual(adhd_id).color, 'MarkerFaceColor', color_map_individual(adhd_id).color, ...
%         'Marker', color_map_individual(adhd_id).marker,...
%         'DisplayName', strcat('id',num2str(adhd_subject_ids{s}))); % Use subject ID
% end
% 
% % Scatter individual data points for non-ADHD
% nonadhd_handles = gobjects(num_nonadhd_subjects, 1); % Store handles for the legend
% for s = 1:num_nonadhd_subjects
%     nonadhd_id = nonadhd_ids(s);
%     nonadhd_handles(s) = scatter(1:numConditions, nonAdhdVelocityData(:, s), 30, ...
%         'MarkerEdgeColor', color_map_individual(nonadhd_id).color, 'MarkerFaceColor', color_map_individual(nonadhd_id).color, ...
%         'Marker', color_map_individual(nonadhd_id).marker,...
%         'DisplayName', strcat('id',num2str(nonadhd_subject_ids{s}))); % Use subject ID
% end
% 
% % Error bars for ADHD
% adhd_errorbar_handle = errorbar(1:numConditions, adhdVelocityMedian, adhdVelocitySem, '-o', ...
%     'Color', color_map('ADHD'), 'MarkerFaceColor', color_map('ADHD'), ...
%     'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'ADHD');
% 
% % Error bars for non-ADHD
% nonadhd_errorbar_handle = errorbar(1:numConditions, nonAdhdVelocityMedian, nonAdhdVelocitySem, '-o', ...
%     'Color', color_map('nonADHD'), 'MarkerFaceColor', color_map('nonADHD'), ...
%     'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'nonADHD');
% 
% legend2 = legend([adhd_errorbar_handle; adhd_handles; nonadhd_errorbar_handle; nonadhd_handles], 'Location', 'northeastoutside');
% %title(legend2, 'Subjects');
% legend2.Box = 'off';
% legend2.Position = [0.92, 0.58, 0.1, 0.2]; % [left, bottom, width, height]
% % Customize the plot
% xticks(1:numConditions);
% xticklabels(condition_labels);
% %xlabel('Condition');
% ylabel('Average Saccade Velocity (pixels/ms)');
% %legend_handles = legend([adhd_handles; nonadhd_handles], 'Location', 'BestOutside'); % Add only scatter handles to the legend
% %title(legend_handles, 'Subjects');
% %legend_handles.Box = 'off';
% title('Velocity: Scatter + Error Bar');
% grid off; hold off;
% 
% % Get the current positions of the subplots
% subplot_left = subplot(1, 2, 1); % Left plot
% subplot_right = subplot(1, 2, 2); % Right plot
% 
% % Get positions of each subplot
% pos_left = get(subplot_left, 'Position');
% pos_right = get(subplot_right, 'Position');
% 
% % Adjust the width of the right subplot to match the left subplot
% % Shrink the width of the right plot if its legend pushes it out
% new_width = pos_left(3); % Use the width of the left plot
% pos_right(3) = new_width; % Set the right plot's width to match
% set(subplot_right, 'Position', pos_right);
% 
% if safe == 1
%     set(gcf, 'Position', [100, 100, 1200, 600]); % Resize the figure window (x, y, width, height)
%     saveas(gcf, strcat(comparison_results_folder, '\01_saccade_amplitude_velocity_comparison.png')); % Save the figure
% end

