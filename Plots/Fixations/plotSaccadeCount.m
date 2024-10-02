function plotSaccadeCount(saccadeStats, group_labels, conditions, color_map)
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
end
