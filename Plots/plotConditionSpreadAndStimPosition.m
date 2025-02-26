function plotConditionSpreadAndStimPosition(rand_trials, num_rows, num_columns, analysis_folder)
    targetPositions = table2array(rand_trials(:, {'TargetPosition'}));%, 'TargetPosition'}));

    conditions = categorical(rand_trials.Condition);
    uniqueConditions = unique(conditions);
    numConditions = length(uniqueConditions);
    maxPosition_y = num_rows;
    maxPosition_x = num_columns;
    
    countMatrices = cell(numConditions, 1);
    for i = 1:numConditions
        %countMatrices{i} = zeros(maxPosition_y, maxPosition_x);
        countMatrices{i} = zeros(maxPosition_y, maxPosition_x);
    end
    
    % Count the occurrences of each target position per condition
    for i = 1:height(rand_trials)
        xPos = targetPositions(i, 1);
        yPos = targetPositions(i, 2);
        conditionIdx = find(uniqueConditions == conditions(i));
        countMatrices{conditionIdx}(yPos, xPos) = countMatrices{conditionIdx}(yPos, xPos) + 1;
        %countMatrices{conditionIdx}(xPos, yPos) = countMatrices{conditionIdx}(xPos, yPos) + 1;
    end

    % Plot the data in a 2x2 subplot layout
    figure(1);
    for i = 1:numConditions
        subplot(2, 2, i);
        imagesc(countMatrices{i});
        colorbar;
        title(['Condition: ', char(uniqueConditions(i))]);
        xlabel('Y Position');
        ylabel('X Position');
        axis equal tight;
        saveas(gcf,strcat(analysis_folder, '\condition_spread.png'));

    end
        % Sum all conditions into one
    totalCountMatrix = sum(cat(3, countMatrices{:}), 3);

    % Plot the summed conditions
    figure(2);
    imagesc(totalCountMatrix);
    colorbar;
    title('Sum of All Conditions');
    xlabel('X Position');
    ylabel('Y Position');
    axis equal tight;
    set(gca, 'YDir', 'normal'); % Ensure the y-axis is in normal direction
    saveas(gcf,strcat(analysis_folder, '\condition_spread_sum.png'));

end