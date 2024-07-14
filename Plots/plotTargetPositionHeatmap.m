function plotTargetPositionHeatmap(rand_trials, num_rows, num_columns)
    % Extract target positions from the table
    num_trials = height(rand_trials);
    targetPositions = zeros(num_trials, 2);
    
    for i = 1:num_trials
        targetPositions(i, :) = cell2mat(rand_trials.TargetPosition(i));
    end
    
    % Initialize heatmap matrix
    heatmapMatrix = zeros(num_rows, num_columns);
    
    for i = 1:size(targetPositions, 1)
        xPos = targetPositions(i, 1);
        yPos = targetPositions(i, 2);
        heatmapMatrix(yPos, xPos) = heatmapMatrix(yPos, xPos) + 1;
    end
    
    % Plot the heatmap
    figure;
    imagesc(heatmapMatrix);
    colorbar;
    xlabel('X Position');
    ylabel('Y Position');
    title('Heatmap of Target Positions');
    axis equal tight;
    set(gca, 'YDir', 'normal'); % Ensure the y-axis is in normal direction
end
