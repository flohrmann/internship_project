function plotDistanceVsRT(trial_results, eye_tracking_data, eye_rt, screenXpixels, screenYpixels, analysis_folder, color_map)
    % Initialize arrays to store distances, reaction times, and conditions
    distances = [];
    reactionTimes = [];
    conditions = [];

    % Loop through each trial in trial_results
    for trial = 1:size(trial_results, 1)
        % Extract the current trial data
        current_data = trial_results(trial, :);

        % Extract the relevant times
        stimulusOnsetTime = current_data.StimulusOnsetTime;
        trial_resp_time = current_data.rt;

        % Extract the target position (center of the stimulus position)
        targetX = current_data.x_centers{1}(current_data.TargetPosition(2), current_data.TargetPosition(1));
        targetY = current_data.y_centers{1}(current_data.TargetPosition(2), current_data.TargetPosition(1)); % Inverting y-axis

        % Extract system timestamps and gaze points for both eyes
        timestamps = double(eye_tracking_data.systemTimeStamp) / 1e6;

        % Find indices for the stimulus onset
        [~, idx_stimulus_onset] = min(abs(timestamps - stimulusOnsetTime));

        % Extract gaze points at stimulus onset for both eyes
        gazeX_right = eye_tracking_data.right.gazePoint.onDisplayArea(1, idx_stimulus_onset) * screenXpixels;
        gazeY_right = screenYpixels - (eye_tracking_data.right.gazePoint.onDisplayArea(2, idx_stimulus_onset) * screenYpixels); % Inverting y-axis

        gazeX_left = eye_tracking_data.left.gazePoint.onDisplayArea(1, idx_stimulus_onset) * screenXpixels;
        gazeY_left = screenYpixels - (eye_tracking_data.left.gazePoint.onDisplayArea(2, idx_stimulus_onset) * screenYpixels); % Inverting y-axis

        % Calculate distances to the target position
        distance_right = sqrt((gazeX_right - targetX).^2 + (gazeY_right - targetY).^2);
        distance_left = sqrt((gazeX_left - targetX).^2 + (gazeY_left - targetY).^2);

        % Calculate reaction time
        reactionTime = trial_resp_time;

        % Store the minimum distance, corresponding reaction time, and condition
        distances = [distances; min(distance_right, distance_left)];
        reactionTimes = [reactionTimes; reactionTime];
        conditions = [conditions; current_data.Condition];
    end

    % Plot the distances vs. reaction times with colors per condition
    figure;
    hold on;

    % Get unique conditions
    unique_conditions = unique(conditions);

    for i = 1:length(unique_conditions)
        cond = unique_conditions{i};
        idx = strcmp(conditions, cond);
        scatter(distances(idx), reactionTimes(idx), 50, 'filled', 'MarkerFaceColor', color_map(cond));
    end

    title('Distance to Target at Onset vs. Reaction Time');
    xlabel('Distance to Target at Onset (pixels)');
    ylabel('Reaction Time (s)');
    grid on;

    % Add a regression line to the plot
    p = polyfit(distances, reactionTimes, 1); % Linear fit
    yfit = polyval(p, distances);
    plot(distances, yfit, '-k', 'LineWidth', 1.5); % Regression line in black

    hold off;

    % Save the plot
    saveas(gcf, fullfile(analysis_folder, 'distance_start_stim_vs_rt.png'));
end
