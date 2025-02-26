function plotBarRTmeanParticipantsConditionPerGroup(data, color_map, comparison_results_folder)
% Identify ADHD and nonADHD groups based on 'group' field
ADHD_indices = strcmp({data.group}, 'ADHD');
nonADHD_indices = strcmp({data.group}, 'nonADHD');

% Extract data for each condition and group using logical indexing
nRTa_ADHD = arrayfun(@(x) x.nRTa, data(ADHD_indices));
nRTa_nonADHD = arrayfun(@(x) x.nRTa, data(nonADHD_indices));
nRTa_simple_ADHD = arrayfun(@(x) x.nRTasimple, data(ADHD_indices));
nRTa_simple_nonADHD = arrayfun(@(x) x.nRTasimple, data(nonADHD_indices));
nRTb_ADHD = arrayfun(@(x) x.nRTb, data(ADHD_indices));
nRTb_nonADHD = arrayfun(@(x) x.nRTb, data(nonADHD_indices));
nRTb_simple_ADHD = arrayfun(@(x) x.nRTbsimple, data(ADHD_indices));
nRTb_simple_nonADHD = arrayfun(@(x) x.nRTbsimple, data(nonADHD_indices));

% Calculate the means and SEMs for each condition and group
mean_values_ADHD = [mean(nRTa_ADHD, 'omitnan'), mean(nRTa_simple_ADHD, 'omitnan'), mean(nRTb_ADHD, 'omitnan'), mean(nRTb_simple_ADHD, 'omitnan')];
sem_values_ADHD = [std(nRTa_ADHD, 'omitnan') / sqrt(length(nRTa_ADHD)), std(nRTa_simple_ADHD, 'omitnan') / sqrt(length(nRTa_simple_ADHD)), std(nRTb_ADHD, 'omitnan') / sqrt(length(nRTb_ADHD)), std(nRTb_simple_ADHD, 'omitnan') / sqrt(length(nRTb_simple_ADHD))];

mean_values_nonADHD = [mean(nRTa_nonADHD, 'omitnan'), mean(nRTa_simple_nonADHD, 'omitnan'), mean(nRTb_nonADHD, 'omitnan'), mean(nRTb_simple_nonADHD, 'omitnan')];
sem_values_nonADHD = [std(nRTa_nonADHD, 'omitnan') / sqrt(length(nRTa_nonADHD)), std(nRTa_simple_nonADHD, 'omitnan') / sqrt(length(nRTa_simple_nonADHD)), std(nRTb_nonADHD, 'omitnan') / sqrt(length(nRTb_nonADHD)), std(nRTb_simple_nonADHD) / sqrt(length(nRTb_simple_nonADHD))];

% Define condition labels and groups
conditions = {'a', 'a_simple', 'b', 'b_simple'};
condition_labels = {'a', 'a simple', 'b', 'b simple'};

% Set up bar positions for each group
bar_width = 0.35;
group_positions = [(1:length(conditions)) - bar_width / 2; (1:length(conditions)) + bar_width / 2];

% Create figure
figure; hold on;

% Plot bars for ADHD and nonADHD groups with SEMs
for i = 1:length(conditions)
    % ADHD bars
    bar(group_positions(1, i), mean_values_ADHD(i), 'FaceColor', color_map('ADHD'), 'BarWidth', bar_width, 'EdgeColor', 'none');
    errorbar(group_positions(1, i), mean_values_ADHD(i), sem_values_ADHD(i), 'k', 'LineStyle', 'none', 'LineWidth', 1, 'HandleVisibility', 'off');
    
    % nonADHD bars
    bar(group_positions(2, i), mean_values_nonADHD(i), 'FaceColor', color_map('nonADHD'), 'BarWidth', bar_width, 'EdgeColor', 'none');
    errorbar(group_positions(2, i), mean_values_nonADHD(i), sem_values_nonADHD(i), 'k', 'LineStyle', 'none', 'LineWidth', 1, 'HandleVisibility', 'off');
end

% Plot individual data points and connect with lines by group color
num_participants_ADHD = length(nRTa_ADHD);
num_participants_nonADHD = length(nRTa_nonADHD);

for i = 1:num_participants_ADHD
    plot(group_positions(1, :), [nRTa_ADHD(i), nRTa_simple_ADHD(i), nRTb_ADHD(i), nRTb_simple_ADHD(i)], ...
        '-o', 'MarkerSize', 5, 'Color', color_map('ADHD'), 'LineWidth', 1, 'HandleVisibility', 'off');
end

for i = 1:num_participants_nonADHD
    plot(group_positions(2, :), [nRTa_nonADHD(i), nRTa_simple_nonADHD(i), nRTb_nonADHD(i), nRTb_simple_nonADHD(i)], ...
        '-o', 'MarkerSize', 5, 'Color', color_map('nonADHD'), 'LineWidth', 1, 'HandleVisibility', 'off');
end

% Organize data for significance testing
all_rt_data_ADHD = {nRTa_ADHD, nRTa_simple_ADHD, nRTb_ADHD, nRTb_simple_ADHD};
all_rt_data_nonADHD = {nRTa_nonADHD, nRTa_simple_nonADHD, nRTb_nonADHD, nRTb_simple_nonADHD};

    
    % Define p-values and add significance lines and stars
y_offset = max([mean_values_ADHD + sem_values_ADHD, mean_values_nonADHD + sem_values_nonADHD]) * 1.1;
for i = 1:length(conditions)
    % Perform t-tests for each condition between ADHD and nonADHD groups
    [~, p] = ttest2(all_rt_data_ADHD{i}, all_rt_data_nonADHD{i}, 'Alpha', 0.05);
    %if p < 0.05
        % Position lines slightly above the bars and adjust for each significance level
        line_y_pos = y_offset * (1 + (i-1) * 0.1);
        plot(group_positions(:, i), [line_y_pos, line_y_pos], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
        
        % Add significance stars based on p-value
        if p < 0.001
            text(mean(group_positions(:, i)), line_y_pos + 0.02 * y_offset, '***', 'FontSize', 12, 'HorizontalAlignment', 'center');
        elseif p < 0.01
            text(mean(group_positions(:, i)), line_y_pos + 02 * y_offset, '**', 'FontSize', 12, 'HorizontalAlignment', 'center');
        elseif p < 0.05
            text(mean(group_positions(:, i)), line_y_pos + 0.02 * y_offset, '*', 'FontSize', 12, 'HorizontalAlignment', 'center');
        else
            text(mean(group_positions(:, i)), line_y_pos + 0.02 * y_offset, num2str(p), 'FontSize', 6, 'HorizontalAlignment', 'center');            
        end
    %end
end

% Set axis labels, ticks, and legend
set(gca, 'XTick', 1:length(conditions), 'XTickLabel', condition_labels);
xlabel('Condition');
ylabel('Reaction Time (RT)');
title('Individual and Average Reaction Time by Condition and Group');
legend({'ADHD', 'nonADHD'}, 'Location', 'northeastoutside');
grid off;
hold off;

set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf, fullfile(comparison_results_folder, 'RT_condition_coloured_groups.png'));