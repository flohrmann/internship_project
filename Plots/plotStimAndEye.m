function plotStimAndEye(analysis_folder, cutData, num_plots)
    close all;
    safe_name = strcat(analysis_folder, '\gaze_path_trial_start');

    %%% get stimulus onset time instead of trialstarttime
%     cutData = renamevars(cutData, 'eyeTrial', 'eyeTrial1');
%     cutData = renamevars(cutData, 'stimulusTrial', 'eyeTrial');
%     safe_name = strcat(analysis_folder, '\gaze_path_stim_onset_');
    %%%%
    
    screenXpixels = 3240;
    screenYpixels = 2160;

    for trial = 1:num_plots%size(cutData,1)
        current_data = cutData(trial,:);
        plot_lines_with_gaze(current_data, screenXpixels, screenYpixels, trial);
        saveas(gcf,strcat(safe_name, num2str(trial),'.png'));
    end
end

function plot_lines_with_gaze(trial_data, screenXpixels, screenYpixels, t)
    % Extracting data for the current trial
    AngleMatrix = trial_data.AngleMatrix{1};
    line_length = trial_data.line_length;
    line_width = 4; %trial_data.line_width;
    x_centers = trial_data.x_centers{1};
    y_centers = trial_data.y_centers{1};

    % Getting the size of the matrices
    [numRows, numCols] = size(AngleMatrix);

    % Extracting the target position
    targetPos = trial_data.TargetPosition;
    targetRow = targetPos(2); % Row index of the target
    targetCol = targetPos(1); % Column index of the target

    % Creating a new figure for each trial
    figure('Position', [0, 0, 3240, 2160]);
    hold on;
    axis equal;

    for i = 1:numRows
        for j = 1:numCols
            angles = AngleMatrix{i, j}; % Cell containing angles

            x_center = x_centers(i, j);
            y_center = screenYpixels - y_centers(i, j);

            % Determine the color based on whether this is the target position
            if i == targetRow && j == targetCol
                line_color = 'y'; % Yellow for the target position
            else
                line_color = 'k'; % Black for other positions
            end

            % Loop through each angle and plot the corresponding line
            for k = 1:length(angles)
                angle = angles(k);

                % Calculate the start and end points of the line
                x_start = x_center - (line_length / 2) * cosd(angle);
                y_start = y_center - (line_length / 2) * sind(angle);
                x_end = x_center + (line_length / 2) * cosd(angle);
                y_end = y_center + (line_length / 2) * sind(angle);

                % Plot the line with the determined color
                plot([x_start, x_end], [y_start, y_end], line_color, 'LineWidth', line_width);
            end
        end
    end

    % Overlay eyetracking data
    visualizeGazePathWithColor(trial_data, screenXpixels, screenYpixels, targetRow, targetCol);
    % Annotate RT from matlab
    text(screenXpixels - 2200, -40, sprintf('RT matlab: %.2f s', trial_data.rt), 'Color', 'black');

    hold off;
    title(['Trial ', num2str(t)]);
end

function visualizeGazePathWithColor(trial_results, screenXpixels, screenYpixels, targetRow, targetCol)
    % Extract gaze points and system timestamps, convert to milliseconds
    
    a_r = trial_results.eyeTrial.right.gazePoint.onDisplayArea;
    x_r = a_r(1,:) * screenXpixels;
    y_r = screenYpixels - (a_r(2,:) * screenYpixels);  % Inverting y-axis
    t_r = double(trial_results.eyeTrial.systemTimeStamp) / 1e6;  % Converts timestamp to seconds

    a_l = trial_results.eyeTrial.left.gazePoint.onDisplayArea;
    x_l = a_l(1,:) * screenXpixels;
    y_l = screenYpixels - (a_l(2,:) * screenYpixels);  % Inverting y-axis
    t_l = double(trial_results.eyeTrial.systemTimeStamp) / 1e6;  % Converts timestamp to seconds

    % Create color gradient maps for both eyes
    numPoints = max([length(x_r), length(x_l)]);
    cm = colormap(jet(numPoints));  % 'jet' ranges from blue to red

    % Plotting both eyes with color gradients
    plotGazePaths(x_r, y_r, cm);
    plotGazePaths(x_l, y_l, cm);

    % Draw a border around the screen dimensions
    rectangle('Position', [0, 0, screenXpixels, screenYpixels], 'EdgeColor', 'k', 'LineWidth', 2, 'LineStyle', '--');

    % Set the axis limits to slightly beyond the full screen dimensions
    xlim([-100 screenXpixels + 100]);
    ylim([-100 screenYpixels + 100]);

    % Annotate time for both eyes
    targetX = trial_results.x_centers{1}(targetRow, targetCol);
    targetY = screenYpixels - trial_results.y_centers{1}(targetRow, targetCol);  % Inverting y-axis for display
    annotateTime(trial_results, x_r, y_r, t_r, x_l, y_l, t_l, screenXpixels, screenYpixels, targetRow, targetCol);


    % Additional plot settings
    title('Gaze Path Visualization with Timing (ms)');
    xlabel('Horizontal Position (pixels)');
    ylabel('Vertical Position (pixels)');
    axis equal;
    grid on;
    
end

function plotGazePaths(x, y, colorMap)
    % Helper function to plot lines with interpolation
    segments = find(~isnan(x) & ~isnan(y));
    if isempty(segments)
        return;  % If all are NaNs, skip
    end

    for i = 1:length(segments)-1
        if segments(i+1) - segments(i) == 1  % Adjacent points are valid
            line(x(segments(i):segments(i+1)), y(segments(i):segments(i+1)), 'Color', colorMap(segments(i),:), 'LineWidth', 2);
        else  % There's a gap
            xi = [x(segments(i)), x(segments(i+1))];
            yi = [y(segments(i)), y(segments(i+1))];
            plot(xi, yi, '--', 'Color', [0.5, 0.5, 0.5]);  % Gray dashed line
        end
    end
end

function annotateTime(trial_results, x_r, y_r, t_r, x_l, y_l, t_l, screenXpixels, screenYpixels, targetRow, targetCol)
    proximityThreshold = 70; % Pixels within which the gaze is considered to reach the target
    targetX = trial_results.x_centers{1}(targetRow, targetCol);
    targetY = screenYpixels - trial_results.y_centers{1}(targetRow, targetCol);

    targetIndex_r = find(sqrt((x_r - targetX).^2 + (y_r - targetY).^2) < proximityThreshold, 1, 'first');
    targetIndex_l = find(sqrt((x_l - targetX).^2 + (y_l - targetY).^2) < proximityThreshold, 1, 'first');

    % Change stimulus color if either eye reaches the target
    if ~isempty(targetIndex_r) || ~isempty(targetIndex_l)
        changeStimulusColor(trial_results, targetRow, targetCol, 'g');
        
        if ~isempty(targetIndex_r)
            timeToTarget_r = t_r(targetIndex_r) - t_r(1);  % Time from trial start to target reach in s
            text(screenXpixels - 500, -40, sprintf('Right Eye Found: %.2f s', timeToTarget_r), 'Color', 'black');
            plot(x_r(targetIndex_r), y_r(targetIndex_r), 'o', 'MarkerFaceColor', 'g');
        end
        
        if ~isempty(targetIndex_l)
            timeToTarget_l = t_l(targetIndex_l) - t_l(1);  % Time from trial start to target reach in s
            text(screenXpixels - 500, -80, sprintf('Left Eye Found: %.2f s', timeToTarget_l), 'Color', 'black');
            plot(x_l(targetIndex_l), y_l(targetIndex_l), 'o', 'MarkerFaceColor', 'g');
        end
    end

    % Determine the last valid data point time for both eyes
    lastValidIndex_r = find(~isnan(x_r) & ~isnan(y_r), 1, 'last');
    lastValidIndex_l = find(~isnan(x_l) & ~isnan(y_l), 1, 'last');

    if lastValidIndex_r
        timeToEnd_r = t_r(lastValidIndex_r) - t_r(1);  % Total time from start to last data point in s
        plot(x_r(lastValidIndex_r), y_r(lastValidIndex_r), 'o', 'MarkerFaceColor', 'k');
        text(0, -40, sprintf('Right Eye Total time: %.2f s', timeToEnd_r), 'Color', 'black');
    end

    if lastValidIndex_l
        timeToEnd_l = t_l(lastValidIndex_l) - t_l(1);  % Total time from start to last data point in s
        plot(x_l(lastValidIndex_l), y_l(lastValidIndex_l), 'o', 'MarkerFaceColor', 'k');
        text(0, -80, sprintf('Left Eye Total time: %.2f s', timeToEnd_l), 'Color', 'black');
    end
end

function changeStimulusColor(trial_results, targetRow, targetCol, newColor)
    % Extracting data for the current trial
    AngleMatrix = trial_results.AngleMatrix{1};
    line_length = trial_results.line_length;
    line_width = 4; % trial_data.line_width;
    x_centers = trial_results.x_centers{1};
    y_centers = trial_results.y_centers{1};
    screenYpixels = 2160;

    % Getting the size of the matrices
    [numRows, numCols] = size(AngleMatrix);

    for i = 1:numRows
        for j = 1:numCols
            angles = AngleMatrix{i, j}; % Cell containing angles

            x_center = x_centers(i, j);
            y_center = screenYpixels - y_centers(i, j);

            % Change the color of the unique stimulus
            if i == targetRow && j == targetCol
                line_color = newColor; % Change color for the target position
            else
                line_color = 'k'; % Black for other positions
            end

            % Loop through each angle and plot the corresponding line
            for k = 1:length(angles)
                angle = angles(k);

                % Calculate the start and end points of the line
                x_start = x_center - (line_length / 2) * cosd(angle);
                y_start = y_center - (line_length / 2) * sind(angle);
                x_end = x_center + (line_length / 2) * cosd(angle);
                y_end = y_center + (line_length / 2) * sind(angle);

                % Plot the line with the determined color
                plot([x_start, x_end], [y_start, y_end], line_color, 'LineWidth', line_width);
            end
        end
    end
end
