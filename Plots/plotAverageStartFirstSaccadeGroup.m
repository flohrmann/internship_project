function plotAverageStartFirstSaccadeGroup(data_struct, ids, unique_conditions, condition_labels, group_labels, color_map,color_map_individual, comparison_results_folder, safe)

% get mean first saccade starrs of participants per condition 
% plot them plus group averages

num_conditions = numel(unique_conditions);
adhd_indices = find(strcmp(group_labels, 'ADHD'));
nonadhd_indices = find(strcmp(group_labels, 'nonADHD'));
adhd_condition_means = zeros(numel(adhd_indices), num_conditions);
nonadhd_condition_means = zeros(numel(nonadhd_indices), num_conditions);

for idx = 1:numel(adhd_indices)
    current_data = data_struct(adhd_indices(idx)); % Extract the participant index & data
    adhd_condition_means(idx, :) = [current_data.first_sacc.a{2}, current_data.first_sacc.b{2}, ...
                                    current_data.first_sacc.a_simple{2}, current_data.first_sacc.b_simple{2}];                              
end
for idx = 1:numel(nonadhd_indices)
    current_data = data_struct(nonadhd_indices(idx)); % Extract the participant index & data
    nonadhd_condition_means(idx, :) = [current_data.first_sacc.a{2}, current_data.first_sacc.b{2}, ...
                                       current_data.first_sacc.a_simple{2}, current_data.first_sacc.b_simple{2}];                              
end


adhd_median = median(adhd_condition_means, 1, 'omitnan');
adhd_sem = std(adhd_condition_means, 0, 1, 'omitnan')./ sqrt(size(adhd_condition_means, 1));
nonAdhd_median = median(nonadhd_condition_means, 1, 'omitnan');
nonAdhd_sem = std(nonadhd_condition_means, 0, 1, 'omitnan')./ sqrt(size(nonadhd_condition_means, 1));

plotADHDnonADHDandDiff('Time until Start of 1st Saccade of Trial',...
    adhd_condition_means, adhd_median, adhd_sem, 'ADHD', 'northeast', ...
    nonadhd_condition_means, nonAdhd_median, nonAdhd_sem, 'nonADHD', 'northeast', ...
    ids, condition_labels, 'Median RT (s)', group_labels, unique_conditions, color_map, color_map_individual, ...
    fullfile(comparison_results_folder, '02_start_first_saccade_allinone_median.png'));




% %% Bar Plot
% num_adhd_subjects = numel(adhd_indices);
% num_nonadhd_subjects = numel(nonadhd_indices);
% 
% adhd_colors = [linspace(0.3, 0, num_adhd_subjects)', linspace(0.8, 0.4, num_adhd_subjects)', linspace(0.3, 0.2, num_adhd_subjects)'];
% nonadhd_colors = [linspace(0.8, 0.4, num_nonadhd_subjects)', linspace(0.3, 0, num_nonadhd_subjects)', linspace(0.3, 0.2, num_nonadhd_subjects)'];
% 
% % Compute group means and standard deviations
% adhd_means = mean(adhd_condition_means, 1, 'omitnan');
% adhd_stds = std(adhd_condition_means, 0, 1, 'omitnan');
% nonadhd_means = mean(nonadhd_condition_means, 1, 'omitnan');
% nonadhd_stds = std(nonadhd_condition_means, 0, 1, 'omitnan');
% 
% figure; hold on;
% adhd_handles = gobjects(num_adhd_subjects, 1);% Scatter individual data points for ADHD
% for s = 1:num_adhd_subjects
%     adhd_handles(s) = scatter(1:num_conditions, adhd_condition_means(s, :), 30, ...
%         'MarkerEdgeColor', adhd_colors(s, :), 'MarkerFaceColor', adhd_colors(s, :), ...
%         'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.3, 'DisplayName', strcat('ID: ', num2str(s))); % ADHD subject scatter
% end
% 
% nonadhd_handles = gobjects(num_nonadhd_subjects, 1);% Scatter individual data points for Non-ADHD
% for s = 1:num_nonadhd_subjects
%     nonadhd_handles(s) = scatter(1:num_conditions, nonadhd_condition_means(s, :), 30, ...
%         'MarkerEdgeColor', nonadhd_colors(s, :), 'MarkerFaceColor', nonadhd_colors(s, :), ...
%         'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.3, 'DisplayName', strcat('ID: ', num2str(s + num_adhd_subjects))); % Non-ADHD subject scatter
% end
% 
% % Error bars for ADHD
% adhd_errorbar_handle = errorbar(1:num_conditions, adhd_means, adhd_stds, '-o', ...
%     'Color', color_map('ADHD'), 'MarkerFaceColor', color_map('ADHD'), ...
%     'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'ADHD Mean ± STD');
% 
% % Error bars for Non-ADHD
% nonadhd_errorbar_handle = errorbar(1:num_conditions, nonadhd_means, nonadhd_stds, '-o', ...
%     'Color', color_map('nonADHD'), 'MarkerFaceColor', color_map('nonADHD'), ...
%     'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'nonADHD Mean ± STD');
% 
% xticks(1:num_conditions); xticklabels(condition_labels);
% ylabel('First Saccade Start Time (s)');
% legend([adhd_errorbar_handle, nonadhd_errorbar_handle], 'Location', 'northeast'); % Legend for means
% title('Time from Stimulation Onset to Start of first Saccade');
% hold off;
% 
% 
% if safe == 1
%     set(gcf, 'Position', [50, 50, 1000, 600]);
%     saveas(gcf, fullfile(comparison_results_folder, '01_saccade_first_start_group_means_individual.png'));
% end

