function plotRTbyEccentricity(trial_results, eye_rt, screenXpixels, screenYpixels, n_rows, n_columns, analysis_folder)
    % Calculate center row and column
    centerRow = round(n_rows / 2);
    centerCol = round(n_columns / 2);
    
    % Define inner positions as a 3x3 block around the center
    
    innerPositions = [
        centerRow-1, centerCol-1;
        centerRow-1, centerCol;
        centerRow-1, centerCol+1;
        centerRow,   centerCol-1;
        centerRow,   centerCol;
        centerRow,   centerCol+1;
        centerRow+1, centerCol-1;
        centerRow+1, centerCol;
        centerRow+1, centerCol+1
    ];

    % Initialize arrays to store RTs for inner and outer circle positions
    RT_inner = [];
    RT_outer = [];

    % Loop through each trial in trial_results
    for trial = 1:size(trial_results, 1)
        % Extract the target position (row and column)
        targetRow = trial_results.TargetPosition(trial, 2);
        targetCol = trial_results.TargetPosition(trial, 1);
        
        % Extract the RT for the current trial
        rt_matlab = eye_rt.RTmatlab(trial);

        % Check if the target position is one of the inner positions
        if ismember([targetRow, targetCol], innerPositions, 'rows')
            RT_inner = [RT_inner; rt_matlab];
        else
            RT_outer = [RT_outer; rt_matlab];
        end
    end

    % Plot the RTs for inner vs. outer circle positions
    figure;
    hold on;
    scatter(ones(size(RT_inner)), RT_inner, 50, 'b', 'filled', 'DisplayName', 'Inner Positions');
    scatter(ones(size(RT_outer)) * 2, RT_outer, 50, 'r', 'filled', 'DisplayName', 'Outer Positions');
    xlim([0 3]);
    xticks([1 2]);
    xticklabels({'Inner Positions', 'Outer Positions'});
    title('Reaction Times by Eccentricity (Inner vs. Outer Positions)');
    xlabel('Eccentricity');
    ylabel('Reaction Time (s)');
    ylim([0 max([RT_inner; RT_outer]) + 0.5]); % Adjust y-limits based on data
    grid on;
    legend show;
    hold off;

    % Save the plot
    saveas(gcf, fullfile(analysis_folder, 'RT_by_eccentricity_positions.png'));
end
