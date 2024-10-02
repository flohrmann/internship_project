function plotGroupFixationDurations(fixationStats, group_labels, conditions, color_map)
    % Initialize variables to store fixation data by group and condition
    numConditions = length(conditions);
    adhdData = zeros(numConditions, 0); % Store data for ADHD group
    nonAdhdData = zeros(numConditions, 0); % Store data for non-ADHD group
    
    % Loop over participants to group data by ADHD status
    for participant = 1:length(fixationStats)
        conditionAverages = zeros(numConditions, 1);
        for c = 1:numConditions
            conditionAverages(c) = fixationStats(participant).conditions(c).avgFixationDuration;
        end
        
        if strcmp(group_labels{participant}, 'ADHD')
            adhdData = [adhdData, conditionAverages];
        else
            nonAdhdData = [nonAdhdData, conditionAverages];
        end
    end
    
    % Calculate means and standard deviations
    adhdMeans = mean(adhdData, 2, 'omitnan');
    adhdStds = std(adhdData, 0, 2, 'omitnan');
    nonAdhdMeans = mean(nonAdhdData, 2, 'omitnan');
    nonAdhdStds = std(nonAdhdData, 0, 2, 'omitnan');
    
    % Plot the results
    figure;
    hold on;
    
    % Plot non-ADHD data
    errorbar(1:numConditions, nonAdhdMeans, nonAdhdStds, '-o', ...
        'Color', color_map('nonADHD'), 'MarkerFaceColor', color_map('nonADHD'), ...
        'DisplayName', 'nonADHD');
    
    % Plot ADHD data
    errorbar(1:numConditions, adhdMeans, adhdStds, '-o', ...
        'Color', color_map('ADHD'), 'MarkerFaceColor', color_map('ADHD'), ...
        'DisplayName', 'ADHD');
    
    % Customize the plot
    xticks(1:numConditions);
    xticklabels(conditions);
    xlabel('Condition');
    ylabel('Average Fixation Duration (ms)');
    legend('Location', 'Best');
    title('Average Fixation Duration per Group and Condition');
    %grid on;
    hold off;
end
