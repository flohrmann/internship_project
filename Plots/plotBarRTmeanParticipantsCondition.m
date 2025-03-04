function plotBarRTmeanParticipantsCondition(data, color_map, comparison_results_folder, color_map_individual, ids)

ADHD_indices = strcmp({data.group}, 'ADHD');
nonADHD_indices = strcmp({data.group}, 'nonADHD');

%eye
prettyPlots(1, data(nonADHD_indices), color_map, fullfile(comparison_results_folder, 'RT_condition_coloured_participants_nonADHD_eye.png'),'nonADHD: Gaze Reaction Time',  color_map_individual, ids(nonADHD_indices));
prettyPlots(1, data(ADHD_indices), color_map, fullfile(comparison_results_folder, 'RT_condition_coloured_participants_ADHD_eye.png'), 'ADHD: Gaze Reaction Time',  color_map_individual, ids(ADHD_indices));
%button
prettyPlots(0, data(nonADHD_indices), color_map, fullfile(comparison_results_folder, 'RT_condition_coloured_participants_nonADHD_button.png'),'nonADHD: Button Press Reaction Time',  color_map_individual, ids(nonADHD_indices));
prettyPlots(0, data(ADHD_indices), color_map, fullfile(comparison_results_folder, 'RT_condition_coloured_participants_ADHD_button.png'), 'ADHD: Button Press Reaction Time',  color_map_individual, ids(ADHD_indices));



end
% Extract data for each condition and group using logical indexing









function prettyPlots(eye, data, color_map, safe_here, this_title, color_map_individual,id)

if eye == 0
    nRTa_all = arrayfun(@(x) x.nRTa, data);
    nRTa_simple_all = arrayfun(@(x) x.nRTasimple, data);
    nRTb_all = arrayfun(@(x) x.nRTb, data);
    nRTb_simple_all = arrayfun(@(x) x.nRTbsimple, data);
else
    nRTa_all = arrayfun(@(x) x.nRTa_eye, data);
    nRTa_simple_all = arrayfun(@(x) x.nRTasimple_eye, data);
    nRTb_all = arrayfun(@(x) x.nRTb_eye, data);
    nRTb_simple_all = arrayfun(@(x) x.nRTbsimple_eye, data);
end

% Calculate the mean and SEM for each condition
mean_rt_a = median(nRTa_all, 'omitnan');
mean_rt_a_simple = median(nRTa_simple_all, 'omitnan');
mean_rt_b = mean(nRTb_all, 'omitnan');
mean_rt_b_simple = median(nRTb_simple_all, 'omitnan');

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
num_participants = size(data, 2);
for i = 1:num_participants
    
    info = color_map_individual(id(i));
    plot(1:4, individual_data(:, i), 'MarkerEdgeColor', info.color, ...
            'MarkerFaceColor', info.color, 'Marker', info.marker,...
            'Color', info.color, ...
            'DisplayName', sprintf('ID%d', id(i)), ...
            'MarkerSize', 5, 'LineWidth', 1);
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
%xlabel('Condition');
ylabel('Reaction Time (s)');
title(this_title);

legend('boxoff'); 
legend('Location', 'northeast');
grid off;
hold off;

set(gcf, 'Position', [50, 50, 700, 500]); % Resize the figure window (x, y, width, height)
saveas(gcf, safe_here);
end