function plotFixationDifferences(fixationStats, group_labels, conditions, color_map)
    % Initialize variables to store the fixation data by group and condition
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
    
    % Calculate the average fixation durations per condition
    adhdMeans = mean(adhdData, 2, 'omitnan');
    nonAdhdMeans = mean(nonAdhdData, 2, 'omitnan');
    
end
