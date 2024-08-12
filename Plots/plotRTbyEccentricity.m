function plotRTbyEccentricity(trial_results, eye_rt, screenXpixels, screenYpixels, analysis_folder)
    % Define a threshold for inner vs. outer circle (in pixels)
    centerX = screenXpixels / 2;
    centerY = screenYpixels / 2;
    distanceThreshold = min(centerX, centerY) * 0.3; % 30% of the distance to the edge for inner circle

    % Initialize arrays to store RTs for inner and outer circle positions
    RT_inner = [];
    RT_outer = [];

    % Loop through each trial in trial_results
    for trial = 1:size(trial_results, 1)
        % Extract the target position (center of the stimulus position)
        targetX = trial_results.x_centers{1}(trial_results.TargetPosition(trial, 2), trial_results.TargetPosition(trial, 1));
        targetY = trial_results.y_centers{1}(trial_results.TargetPosition(trial, 2), trial_results.TargetPosition(trial, 1));
        
        % Calculate the distance from the center of the screen
        distanceFromCenter = sqrt((targetX - centerX)^2 + (targetY - centerY)^2);

        % Extract the RT for the current trial
        rt_matlab = eye_rt.RTmatlab(trial);

        % Categorize RT based on the distance from the center
        if distanceFromCenter <= distanceThreshold
            RT_inner = [RT_inner; rt_matlab];
        else
            RT_outer = [RT_outer; rt_matlab];
        end
    end

    % Plot the RTs for inner vs. outer circle positions
    figure;
    hold on;
    scatter(ones(size(RT_inner)), RT_inner, 'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'Inner Circle');
    scatter(ones(size(RT_outer)) * 2, RT_outer, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'Outer Circle');
    xlim([0 3]);
    xticks([1 2]);
    xticklabels({'Inner Circle', 'Outer Circle'});
    title('Reaction Times by Eccentricity (Inner vs. Outer Circle)');
    xlabel('Eccentricity');
    ylabel('Reaction Time (s)');
    grid on;
    legend show;
    hold off;

    % Save the plot
    saveas(gcf, fullfile(analysis_folder, 'RT_by_eccentricity.png'));
end
