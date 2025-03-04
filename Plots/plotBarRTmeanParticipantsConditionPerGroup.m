function plotBarRTmeanParticipantsConditionPerGroup(data, color_map, comparison_results_folder, color_map_individual)
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
mean_values_ADHD = [median(nRTa_ADHD, 'omitnan'),  median(nRTb_ADHD, 'omitnan'),median(nRTa_simple_ADHD, 'omitnan'), median(nRTb_simple_ADHD, 'omitnan')];
sem_values_ADHD = [std(nRTa_ADHD, 'omitnan') / sqrt(length(nRTa_ADHD)),  std(nRTb_ADHD, 'omitnan') / sqrt(length(nRTb_ADHD)), std(nRTa_simple_ADHD, 'omitnan') / sqrt(length(nRTa_simple_ADHD)),std(nRTb_simple_ADHD, 'omitnan') / sqrt(length(nRTb_simple_ADHD))];

mean_values_nonADHD = [median(nRTa_nonADHD, 'omitnan'),  median(nRTb_nonADHD, 'omitnan'),median(nRTa_simple_nonADHD, 'omitnan'), median(nRTb_simple_nonADHD, 'omitnan')];
sem_values_nonADHD = [std(nRTa_nonADHD, 'omitnan') / sqrt(length(nRTa_nonADHD)), std(nRTb_nonADHD, 'omitnan') / sqrt(length(nRTb_nonADHD)), std(nRTa_simple_nonADHD, 'omitnan') / sqrt(length(nRTa_simple_nonADHD)), std(nRTb_simple_nonADHD) / sqrt(length(nRTb_simple_nonADHD))];

% Define condition labels and groups
conditions = {'a', 'b', 'a_simple', 'b_simple'};
condition_labels = {'a', 'b','a simple',  'b simple'};

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
for i = 1:size(data, 2)
    info = color_map_individual(data(i).id);
    if strcmp(data(i).group, 'ADHD')
         plot(group_positions(1, :), [data(i).nRTa, data(i).nRTb,data(i).nRTasimple, data(i).nRTbsimple], ...
        'MarkerFaceColor', info.color, 'Marker', info.marker,...
            'MarkerEdgeColor', info.color, 'Color', info.color, ...
            'DisplayName', sprintf('ID%d', data(i).id), ...
            'MarkerSize', 5, 'LineWidth', 1);
    else %strcmp(group_labels{i}, 'ADHD')
         plot(group_positions(2, :), [data(i).nRTa,  data(i).nRTb, data(i).nRTasimple,data(i).nRTbsimple], ...
        'MarkerFaceColor', info.color, 'Marker', info.marker,...
            'MarkerEdgeColor', info.color, 'Color', info.color, ...
            'DisplayName', sprintf('ID%d', data(i).id), ...
            'MarkerSize', 5, 'LineWidth', 1);
    end
   
end

% Organize data for significance testing
all_rt_data_ADHD = {nRTa_ADHD,  nRTb_ADHD, nRTa_simple_ADHD,nRTb_simple_ADHD};
all_rt_data_nonADHD = {nRTa_nonADHD,  nRTb_nonADHD, nRTa_simple_nonADHD,nRTb_simple_nonADHD};

    
    % Define p-values and add significance lines and stars
y_offset = max([mean_values_ADHD + sem_values_ADHD, mean_values_nonADHD + sem_values_nonADHD]) * 1.1;
for i = 1:length(conditions)
    % Perform t-tests for each condition between ADHD and nonADHD groups
    [~, p] = ttest2(all_rt_data_ADHD{i}, all_rt_data_nonADHD{i}, 'Alpha', 0.05);
    if p < 0.05
        % Position lines slightly above the bars and adjust for each significance level
        line_y_pos = y_offset * (1 + (i-1) * 0.1);
        plot(group_positions(:, i), [line_y_pos, line_y_pos], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
        
        % Add significance stars based on p-value
        if p < 0.001
            text(mean(group_positions(:, i)), line_y_pos + 0.04 * y_offset, [num2str(p), '***'], 'FontSize', 8, 'HorizontalAlignment', 'center');
        elseif p < 0.01
            text(mean(group_positions(:, i)), line_y_pos + 004 * y_offset, [num2str(p), '**'], 'FontSize', 8, 'HorizontalAlignment', 'center');
        elseif p < 0.05
            text(mean(group_positions(:, i)), line_y_pos + 0.04 * y_offset, [num2str(p), '*'], 'FontSize', 8, 'HorizontalAlignment', 'center');
        else
            %text(mean(group_positions(:, i)), line_y_pos + 0.04 * y_offset, num2str(p), 'FontSize', 8, 'HorizontalAlignment', 'center');            
        end
    end
end

% Set axis labels, ticks, and legend
set(gca, 'XTick', 1:length(conditions), 'XTickLabel', condition_labels);
%xlabel('Condition');
ylabel('Button Press Reaction Time (s)');
title('Individual and Median Reaction Time by Condition and Group');
legend({'ADHD', 'nonADHD'}, 'Location', 'northeast');
grid off;

legend('boxoff');
hold off;

%set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
set(gcf, 'Position', [50, 50, 700, 700]); % Resize the figure window (x, y, width, height)
saveas(gcf, fullfile(comparison_results_folder, 'RT_condition_coloured_groups.png'));