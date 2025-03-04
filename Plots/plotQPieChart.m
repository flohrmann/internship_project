function plotQPieChart(adhd_numeric, non_adhd_numeric, safe, comparison_results_folder)
% Function to plot a pie chart of response proportions for a group
%
% Inputs:
%   data_numeric: Numeric responses for the group
%   group_name: Name of the group ('ADHD' or 'non-ADHD')
%   color_map: Color mapping for groups


% Plot the pie chart
figure; sgtitle('ASRS Scores');

subplot(2, 2, 1);
plotPie(adhd_numeric, 'ADHD: all')
%plotQPieChart(adhd_numeric, 'ADHD', color_map, 'all', safe, comparison_results_folder);
subplot(2, 2, 2);
plotPie(non_adhd_numeric, 'nonADHD: all')
%plotQPieChart(non_adhd_numeric, 'nonADHD', color_map, 'all', safe, comparison_results_folder);

% Pie Chart for first 6 responses: proportions for (non)ADHD group
subplot(2, 2, 3);
plotPie(adhd_numeric(:,1:6), 'ADHD: first 6')
%plotQPieChart(adhd_numeric(:,1:6), 'ADHD', color_map, '6', safe, comparison_results_folder);
subplot(2, 2, 4);
plotPie(non_adhd_numeric(:,1:6), 'nonADHD: first 6')
%plotQPieChart(non_adhd_numeric(:,1:6), 'nonADHD', color_map, '6', safe, comparison_results_folder);


legend({'Never', 'Rarely', 'Sometimes', 'Often', 'Very Often'}, 'Location', 'bestoutside');
title('ASRS Response Proportions');

if safe == 1
    %set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf, fullfile(comparison_results_folder, strcat('quest_pie_answer_distr_all.png')));
end
end




function plotStandardTwoGroupsEachData(data, medians, sems, group, y_label, position, ids, group_labels, conditions, labels, color_map, participant_map, sharedYLimits)
    % Plot group medians and scatter individual participant data
    numConditions = length(medians);
    hold on;
    
    % Bar plot for medians with error bars
    for c = 1:numConditions
        bar(c, medians(c), 'FaceColor', color_map(conditions{c}), 'EdgeColor', 'none', ...
            'FaceAlpha', 0.3, 'DisplayName', group, 'HandleVisibility', 'off');
    end
    
    current_indices = strcmp(group_labels, group);
    group_ids = ids(current_indices);

    % Scatter individual participant data
    jitterAmount = 0.4; 
    for participant = 1:length(group_ids)
        participant_id = group_ids(participant);

        %if strcmp(group_labels{participant}, group)
            info = participant_map(participant_id); % to retrieve color and marker
            current_data = data(participant,:);
            jitteredX = (1:numConditions) + jitterAmount * (rand(1, numConditions) - 0.5); % Jitter x-coordinates
            
            scatter(jitteredX, current_data, 40, ... % Use scatter for points only
                'MarkerEdgeColor', info.color, ...
                'MarkerFaceColor', info.color, ...
                'Marker', info.marker, ...
                'DisplayName', sprintf('ID%d', participant_id));
        %end
    end
    
    for c = 1:numConditions % errorbar on top for better visibility
       errorbar(c, medians(c), sems(c), '.', 'LineWidth', 1.5, 'CapSize', 10, 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    end
    xticks(1:numConditions); xticklabels(labels);
    ylabel(y_label);  ylim(sharedYLimits);
    legend('boxoff');
    legend('Location', position, 'Orientation', 'horizontal', 'NumColumns', 2); % Two-column legend
    title(group);    hold off;
end





function plotPie(data_numeric, this_title)
    % Count occurrences of each response category
    title(this_title);
    counts_total = histcounts(data_numeric(:), 0.5:1:5.5)
    pie(counts_total);

end