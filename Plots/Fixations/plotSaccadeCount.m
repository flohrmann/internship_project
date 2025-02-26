function plotSaccadeCount(saccadeStats, group_labels, conditions, color_map, ids, color_map_individual, comparison_results_folder, safe)
    numConditions = length(conditions);
    
    % Initialize matrices to store group data
    adhdSaccadeCounts = zeros(numConditions, 0);
    nonAdhdSaccadeCounts = zeros(numConditions, 0);
    
    % Loop over participants to group data by ADHD status
    for participant = 1:length(saccadeStats)
        conditionSaccadeCounts = zeros(numConditions, 1);
        
        for c = 1:numConditions
            conditionSaccadeCounts(c) = saccadeStats(participant).conditions(c).avgSaccadeCount;
        end
        
        if strcmp(group_labels{participant}, 'ADHD')
            adhdSaccadeCounts = [adhdSaccadeCounts, conditionSaccadeCounts];
        else
            nonAdhdSaccadeCounts = [nonAdhdSaccadeCounts, conditionSaccadeCounts];
        end
    end
    
    % Calculate the mean and standard deviation for each condition and group
    adhdSaccadeCountMeans = mean(adhdSaccadeCounts, 2, 'omitnan');
    adhdSaccadeCountStds = std(adhdSaccadeCounts, 0, 2, 'omitnan');
    nonAdhdSaccadeCountMeans = mean(nonAdhdSaccadeCounts, 2, 'omitnan');
    nonAdhdSaccadeCountStds = std(nonAdhdSaccadeCounts, 0, 2, 'omitnan');
    
    % Plot Saccade Count Differences
    figure;
    hold on;
    
    
            % Add individual scatter points for non-ADHD participants
    nonadhd_indices = strcmp(group_labels, 'nonADHD');
    nonadhd_ids = ids(nonadhd_indices);
    
    
    %% todo add jitter
%     jitterAmount = 0.1; % Adjust for marker spread
%     for participant = 1:size(all_means_condition, 1)
%         for condition = 1:num_conditions
%             x = condition + (rand(1) - 0.5) * jitterAmount; % Jitter x-coordinates
%             
            
            
    for i = 1:length(nonadhd_ids)
        participant_id = nonadhd_ids(i);
        scatter(1:numConditions, nonAdhdSaccadeCounts(:, i), 50, ...
            'MarkerEdgeColor', color_map_individual(participant_id).color, ...
            'MarkerFaceColor', color_map_individual(participant_id).color, ...
            'Marker', color_map_individual(participant_id).marker, ...
            'DisplayName', strcat('ID ', num2str(participant_id)));
    end

    % Add individual scatter points for ADHD participants
    adhd_indices = strcmp(group_labels, 'ADHD');
    adhd_ids = ids(adhd_indices);
    for i = 1:length(adhd_ids)
        participant_id = adhd_ids(i);
        scatter(1:numConditions, adhdSaccadeCounts(:, i), 50, ...
            'MarkerEdgeColor', color_map_individual(participant_id).color, ...
            'MarkerFaceColor', color_map_individual(participant_id).color, ...
            'Marker', color_map_individual(participant_id).marker, ...
            'DisplayName', strcat('ID ', num2str(participant_id)));
    end
    
    
    % Plot non-ADHD data
    errorbar(1:numConditions, nonAdhdSaccadeCountMeans, nonAdhdSaccadeCountStds, '-o', ...
        'Color', color_map('nonADHD'), 'MarkerFaceColor', color_map('nonADHD'), ...
        'DisplayName', 'nonADHD');
    
    % Plot ADHD data
    errorbar(1:numConditions, adhdSaccadeCountMeans, adhdSaccadeCountStds, '-o', ...
        'Color', color_map('ADHD'), 'MarkerFaceColor', color_map('ADHD'), ...
        'DisplayName', 'ADHD');
    
    
    
    

    
    
    
    
    % Customize the plot
    xticks(1:numConditions);
    xticklabels(conditions);
    xlabel('Condition');
    ylabel('Average Saccade Count');
    legend('Location', 'Best');
    title('Saccade Count Differences between ADHD and nonADHD');
    hold off;
    
    if safe == 1
    set(gcf, 'Position', [100, 100, 1200, 600]); % Resize the figure window (x, y, width, height)
    saveas(gcf, strcat(comparison_results_folder, '\01_saccade_count_group.png')); % Save the figure
    end
    
    
end
