function plotSaccadeDifferences(saccadeStats, group_labels, conditions, color_map)
    numConditions = length(conditions);
    
    % Initialize matrices to store group data
    adhdAmplitudeData = zeros(numConditions, 0);
    nonAdhdAmplitudeData = zeros(numConditions, 0);
    adhdVelocityData = zeros(numConditions, 0);
    nonAdhdVelocityData = zeros(numConditions, 0);
    
    % Loop over participants to group data by ADHD status
    for participant = 1:length(saccadeStats)
        conditionAmplitudes = zeros(numConditions, 1);
        conditionVelocities = zeros(numConditions, 1);
        
        for c = 1:numConditions
            conditionAmplitudes(c) = saccadeStats(participant).conditions(c).avgSaccadeAmplitude;
            conditionVelocities(c) = saccadeStats(participant).conditions(c).avgSaccadeVelocity;
        end
        
        if strcmp(group_labels{participant}, 'ADHD')
            adhdAmplitudeData = [adhdAmplitudeData, conditionAmplitudes];
            adhdVelocityData = [adhdVelocityData, conditionVelocities];
        else
            nonAdhdAmplitudeData = [nonAdhdAmplitudeData, conditionAmplitudes];
            nonAdhdVelocityData = [nonAdhdVelocityData, conditionVelocities];
        end
    end
    
    % Calculate the mean and standard deviation for each condition and group
    adhdAmplitudeMeans = mean(adhdAmplitudeData, 2, 'omitnan');
    adhdAmplitudeStds = std(adhdAmplitudeData, 0, 2, 'omitnan');
    nonAdhdAmplitudeMeans = mean(nonAdhdAmplitudeData, 2, 'omitnan');
    nonAdhdAmplitudeStds = std(nonAdhdAmplitudeData, 0, 2, 'omitnan');
    
    adhdVelocityMeans = mean(adhdVelocityData, 2, 'omitnan');
    adhdVelocityStds = std(adhdVelocityData, 0, 2, 'omitnan');
    nonAdhdVelocityMeans = mean(nonAdhdVelocityData, 2, 'omitnan');
    nonAdhdVelocityStds = std(nonAdhdVelocityData, 0, 2, 'omitnan');
    
    % Plot Saccade Amplitude Differences
    figure;
    hold on;
    
    % Plot non-ADHD data
    errorbar(1:numConditions, nonAdhdAmplitudeMeans, nonAdhdAmplitudeStds, '-o', ...
        'Color', color_map('nonADHD'), 'MarkerFaceColor', color_map('nonADHD'), ...
        'DisplayName', 'nonADHD');
    
    % Plot ADHD data
    errorbar(1:numConditions, adhdAmplitudeMeans, adhdAmplitudeStds, '-o', ...
        'Color', color_map('ADHD'), 'MarkerFaceColor', color_map('ADHD'), ...
        'DisplayName', 'ADHD');
    
    % Customize the plot
    xticks(1:numConditions);
    xticklabels(conditions);
    xlabel('Condition');
    ylabel('Average Saccade Amplitude (pixels)');
    legend('Location', 'Best');
    title('Saccade Amplitude Differences between ADHD and nonADHD');
    grid on;
    hold off;

    % Plot Saccade Velocity Differences
    figure;
    hold on;
    
    % Plot non-ADHD data
    errorbar(1:numConditions, nonAdhdVelocityMeans, nonAdhdVelocityStds, '-o', ...
        'Color', color_map('nonADHD'), 'MarkerFaceColor', color_map('nonADHD'), ...
        'DisplayName', 'nonADHD');
    
    % Plot ADHD data
    errorbar(1:numConditions, adhdVelocityMeans, adhdVelocityStds, '-o', ...
        'Color', color_map('ADHD'), 'MarkerFaceColor', color_map('ADHD'), ...
        'DisplayName', 'ADHD');
    
    % Customize the plot
    xticks(1:numConditions);
    xticklabels(conditions);
    xlabel('Condition');
    ylabel('Average Saccade Velocity (pixels/ms)');
    legend('Location', 'Best');
    title('Saccade Velocity Differences between ADHD and nonADHD');
    grid on;
    hold off;
end
