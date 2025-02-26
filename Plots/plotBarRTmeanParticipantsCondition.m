function plotBarRTmeanParticipantsCondition(data, color_map, comparison_results_folder)


nRTa_all = arrayfun(@(x) x.nRTa, data);
nRTa_simple_all = arrayfun(@(x) x.nRTasimple, data);
nRTb_all = arrayfun(@(x) x.nRTb, data);
nRTb_simple_all = arrayfun(@(x) x.nRTbsimple, data);

% Calculate the mean and SEM for each condition
mean_rt_a = mean(nRTa_all, 'omitnan');
mean_rt_a_simple = mean(nRTa_simple_all, 'omitnan');
mean_rt_b = mean(nRTb_all, 'omitnan');
mean_rt_b_simple = mean(nRTb_simple_all, 'omitnan');

sem_rt_a = std(nRTa_all, 'omitnan') / sqrt(length(nRTa_all));
sem_rt_a_simple = std(nRTa_simple_all, 'omitnan') / sqrt(length(nRTa_simple_all));
sem_rt_b = std(nRTb_all, 'omitnan') / sqrt(length(nRTb_all));
sem_rt_b_simple = std(nRTb_simple_all, 'omitnan') / sqrt(length(nRTb_simple_all));

% Set up condition names and labels
conditions = {'a', 'a_simple', 'b', 'b_simple'};
condition_labels = {'a', 'a simple', 'b', 'b simple'};
mean_values = [mean_rt_a, mean_rt_a_simple, mean_rt_b, mean_rt_b_simple];
sem_values = [sem_rt_a, sem_rt_a_simple, sem_rt_b, sem_rt_b_simple];
individual_data = [nRTa_all; nRTa_simple_all; nRTb_all; nRTb_simple_all];  % Rows represent conditions

% Pairwise t-tests for significance
num_conditions = 4;
p_values = NaN(num_conditions);  % Matrix to hold p-values
for i = 1:num_conditions
    for j = i+1:num_conditions
        [~, p] = ttest(individual_data(i, :), individual_data(j, :), 'Alpha', 0.05, 'Tail', 'both');
        p_values(i, j) = p;
    end
end

% Define significance levels based on p-values
significance_levels = cell(size(p_values));
for i = 1:num_conditions
    for j = i+1:num_conditions
        if p_values(i, j) < 0.001
            significance_levels{i, j} = '***';  % Highly significant
        elseif p_values(i, j) < 0.01
            significance_levels{i, j} = '**';   % Significant
        elseif p_values(i, j) < 0.05
            significance_levels{i, j} = '*';    % Mildly significant
        end
    end
end

% Create the figure
figure; hold on;

% Plot the bars with SEMs
for i = 1:length(conditions)
    bar(i, mean_values(i), 'FaceColor', color_map(conditions{i}), 'EdgeColor', 'none', 'HandleVisibility', 'off');
    errorbar(i, mean_values(i), sem_values(i), 'k', 'LineStyle', 'none', 'LineWidth', 1, 'HandleVisibility', 'off');
end

% Assign unique colors to each participant and plot individual lines
num_participants = size(individual_data, 2);
cmap = lines(num_participants);  % Generate unique colors for participants
for i = 1:num_participants
    plot(1:4, individual_data(:, i), '-o', 'MarkerSize', 5, ...
        'MarkerFaceColor', cmap(i, :), 'Color', cmap(i, :), 'LineWidth', 1);
end

% Add significance stars with horizontal lines close to the x-axis
y_offset = max(mean_values + sem_values) * 1.5;  % Offset above the bars
for i = 1:num_conditions
    for j = i+1:num_conditions
        if ~isempty(significance_levels{i, j})
            % Calculate x positions for line
            x1 = i;
            x2 = j;
            
            % Define y position for the significance line, adjusting for overlap
            y_line_pos = + y_offset - (j - i) * 0.05 * abs(y_offset);  % Stagger lines by a small amount
            
            % Draw line and place significance star above the line
            plot([x1, x2], [y_line_pos, y_line_pos], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');  % Horizontal line
            text(mean([x1, x2]), y_line_pos + 0.02 * abs(y_offset), significance_levels{i, j}, ...
                 'FontSize', 12, 'HorizontalAlignment', 'center', 'HandleVisibility', 'off');
        end
    end
end

% Configure the plot
set(gca, 'XTick', 1:4, 'XTickLabel', condition_labels);
xlabel('Condition');
ylabel('Reaction Time (RT)');
title('Reaction Time by Condition');
legend('Location', 'northeastoutside');
grid off;
hold off;

set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf, fullfile(comparison_results_folder, 'RT_condition_coloured_participants.png'));
