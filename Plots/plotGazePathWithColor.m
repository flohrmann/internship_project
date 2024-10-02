function [right_eye_arrival_time, left_eye_arrival_time] = plotGazePathWithColor(trial_results, screenXpixels, screenYpixels, targetRow, targetCol)
    % Extract gaze points and system timestamps, convert to milliseconds
    
    a_r = trial_results.eyeTrial.right.gazePoint.onDisplayArea;
    x_r = a_r(1,:) * screenXpixels;
    y_r = screenYpixels - (a_r(2,:) * screenYpixels);  % Inverting y-axis
    t_r = double(trial_results.eyeTrial.systemTimeStamp) / 1e6;  % Converts timestamp to seconds

    a_l = trial_results.eyeTrial.left.gazePoint.onDisplayArea;
    x_l = a_l(1,:) * screenXpixels;
    y_l = screenYpixels - (a_l(2,:) * screenYpixels);  % Inverting y-axis
    t_l = double(trial_results.eyeTrial.systemTimeStamp) / 1e6;  % Converts timestamp to seconds

    % Calculate the mean gaze positions
    x_mean = mean([x_r; x_l], 1);
    y_mean = mean([y_r; y_l], 1);

    [saccade_onsets, saccade_offsets] = detectSaccades(x_mean, y_mean, t_r);  % Replace with your saccade detection function

    % Plot each saccade with a different color
    colormap_lines = lines(length(saccade_onsets));  % Generate a colormap with different colors

    for i = 1:length(saccade_onsets)
        plot(x_mean(saccade_onsets(i):saccade_offsets(i)), y_mean(saccade_onsets(i):saccade_offsets(i)), ...
             'Color', colormap_lines(i,:), 'LineWidth', 2);
    end

    % Draw a border around the screen dimensions
    rectangle('Position', [0, 0, screenXpixels, screenYpixels], 'EdgeColor', 'k', 'LineWidth', 2, 'LineStyle', '--');

    % Set the axis limits to slightly beyond the full screen dimensions
    xlim([-100 screenXpixels + 100]);
    ylim([-100 screenYpixels + 100]);

    % Annotate time for both eyes
    %[right_eye_arrival_time, left_eye_arrival_time] = annotateTime(trial_results, x_r, y_r, t_r, x_l, y_l, t_l, screenXpixels, screenYpixels, targetRow, targetCol);

    % Additional plot settings
    title('Gaze Path Visualization with Timing (ms) and Saccades Colored');
    xlabel('Horizontal Position (pixels)');
    ylabel('Vertical Position (pixels)');
    axis equal;
    grid on;
end

function [saccade_onsets, saccade_offsets] = detectSaccades(x_mean, y_mean, t)
    % Placeholder for saccade detection
    % You should replace this with your actual saccade detection logic
    velocity_threshold = 30;  % Example threshold in degrees/s
    velocity = sqrt(diff(x_mean).^2 + diff(y_mean).^2) ./ diff(t);  % Calculate velocity
    saccades = velocity > velocity_threshold;  % Detect saccades

    saccade_onsets = find(diff([0 saccades]) == 1);
    saccade_offsets = find(diff([saccades 0]) == -1);
end
