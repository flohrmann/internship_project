function plot_this(trial_results)
close all;


%% use once for getting screen width/height
% PsychDebugWindowConfiguration(0,0.5); % TODO disable, makes window transparent
% Screen('Preference', 'SkipSyncTests', 1); % TODO disable, skips synchronization tests
% [window, windowRect] = PsychImaging('OpenWindow', 0, WhiteIndex(0));
% [screenXpixels, screenYpixels] = Screen('WindowSize', window);
% disp(screenXpixels) % 1920
% disp(screenYpixels) % 1080
%% plot
screenXpixels = 1920;
screenYpixels = 1080;
visualizeGazePathWithColorFix(trial_results(1,:).sampFixAll, screenXpixels, screenYpixels);
visualizeGazePathWithColorFix(trial_results(2,:).sampFixAll, screenXpixels, screenYpixels);
visualizeGazePathWithColorFix(trial_results(3,:).sampFixAll, screenXpixels, screenYpixels);
visualizeGazePathWithColorFix(trial_results(4,:).sampFixAll, screenXpixels, screenYpixels);

visualizeGazePathWithColor(trial_results(1,:), screenXpixels, screenYpixels);
visualizeGazePathWithColor(trial_results(2,:), screenXpixels, screenYpixels);
visualizeGazePathWithColor(trial_results(3,:), screenXpixels, screenYpixels);
visualizeGazePathWithColor(trial_results(4,:), screenXpixels, screenYpixels);

end

function visualizeGazePathWithColor(trial_results, screenXpixels, screenYpixels)
    % Extract gaze points and system timestamps, convert to milliseconds
    a_r = trial_results.samp.right.gazePoint.onDisplayArea;
    x_r = a_r(1,:) * screenXpixels;
    y_r = screenYpixels - (a_r(2,:) * screenYpixels);  % Inverting y-axis
    t_r = trial_results.samp.systemTimeStamp / 1e3;  % Converts timestamp to milliseconds

    a_l = trial_results.samp.left.gazePoint.onDisplayArea;
    x_l = a_l(1,:) * screenXpixels;
    y_l = screenYpixels - (a_l(2,:) * screenYpixels);  % Inverting y-axis
    t_l = trial_results.samp.systemTimeStamp / 1e3;  % Converts timestamp to milliseconds

    % Create a figure
    figure;
    hold on;

    % Create color gradient maps for both eyes
    numPoints = max([length(x_r), length(x_l)]);
    cm = colormap(jet(numPoints));  % 'jet' ranges from blue to red

    % Plotting both eyes with color gradients
    plotGazePaths(x_r, y_r, cm);
    plotGazePaths(x_l, y_l, cm);

    % Draw a star at the target position
    targetPos = trial_results.TargetPosition;  
    targetRow = targetPos(1); % Row index of the target
    targetCol = targetPos(2);     % Column index of the target
    targetX = trial_results.x_centers{1}(targetRow, targetCol);
    targetY = screenYpixels - trial_results.y_centers{1}(targetRow, targetCol);  % Inverting y-axis for display
    plot(targetX, targetY, 'p', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y'); % plot star
    
    % Annotate time and proximity for both eyes
    annotateTime(x_r, y_r, t_r, screenXpixels,screenYpixels, 'green', targetX, targetY); % only for right eye should be enough
    %annotateTime(x_l, y_l, t_l, 'red', targetX, targetY);


    % Draw a border around the screen dimensions
    rectangle('Position', [0, 0, screenXpixels, screenYpixels], 'EdgeColor', 'k', 'LineWidth', 2, 'LineStyle', '--');

    % Draw a border around the screen dimensions
    rectangle('Position', [0, 0, screenXpixels, screenYpixels], 'EdgeColor', 'k', 'LineWidth', 2, 'LineStyle', '--');

    % Set the axis limits to slightly beyond the full screen dimensions
    xlim([-100 screenXpixels + 100]);
    ylim([-100 screenYpixels + 100]);

    % Additional plot settings
    title('Gaze Path Visualization with Timing (ms)');
    xlabel('Horizontal Position (pixels)');
    ylabel('Vertical Position (pixels)');
    axis equal;
    grid on;
    hold off;
end



function visualizeGazePathWithColorFix(samp, screenXpixels, screenYpixels)
    % Extract gaze points and system timestamps, convert to milliseconds
    a_r = samp.right.gazePoint.onDisplayArea;
    x_r = a_r(1,:) * screenXpixels;
    y_r = screenYpixels - (a_r(2,:) * screenYpixels);  % Inverting y-axis
    t_r = samp.systemTimeStamp / 1e3;  % Converts timestamp to milliseconds

    a_l = samp.left.gazePoint.onDisplayArea;
    x_l = a_l(1,:) * screenXpixels;
    y_l = screenYpixels - (a_l(2,:) * screenYpixels);  % Inverting y-axis
    t_l = samp.systemTimeStamp / 1e3;  % Converts timestamp to milliseconds

    % Create a figure
    figure;
    hold on;

    % Create color gradient maps for both eyes
    numPoints = max([length(x_r), length(x_l)]);
    cm = colormap(jet(numPoints));  % 'jet' ranges from blue to red

    % Plotting both eyes with color gradients
    plotGazePaths(x_r, y_r, cm);
    plotGazePaths(x_l, y_l, cm);

    % Draw a border around the screen dimensions
    rectangle('Position', [0, 0, screenXpixels, screenYpixels], 'EdgeColor', 'k', 'LineWidth', 2, 'LineStyle', '--');

    % Draw a border around the screen dimensions
    rectangle('Position', [0, 0, screenXpixels, screenYpixels], 'EdgeColor', 'k', 'LineWidth', 2, 'LineStyle', '--');

        lastValidIndex = find(~isnan(x_r) & ~isnan(y_r), 1, 'last');
    if lastValidIndex
        timeToEnd = t_r(lastValidIndex) - t_r(1);  % Total time from start to last data point in ms
        plot(x_r(lastValidIndex), y_r(lastValidIndex), 'o', 'MarkerFaceColor', 'green');
        %text(x(lastValidIndex), y(lastValidIndex) + 20, sprintf('Last Time: %.2f ms', t(lastValidIndex)), 'Color', eyeColor);
        text(0, -40, sprintf('Total time: %.2f ms', timeToEnd), 'Color', 'black');
        %return;
    end

    % Set the axis limits to slightly beyond the full screen dimensions
    xlim([-100 screenXpixels + 100]);
    ylim([-100 screenYpixels + 100]);

    % Additional plot settings
    title('Gaze Path Visualization with Timing (ms)');
    xlabel('Horizontal Position (pixels)');
    ylabel('Vertical Position (pixels)');
    axis equal;
    grid on;
    hold off;
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


function annotateTime(x, y, t, screenXpixels, screenYpixels, eyeColor, targetX, targetY)
    proximityThreshold = 20; % Pixels within which the gaze is considered to reach the target
    targetIndex = find(sqrt((x - targetX).^2 + (y - targetY).^2) < proximityThreshold, 1, 'first');

    if ~isempty(targetIndex)
        timeToTarget = t(targetIndex) - t(1);  % Time from trial start to target reach in ms
        plot(targetX, targetY, 'p', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', eyeColor);
        %text(targetX, targetY - 20, sprintf('Reached at: %.2f ms', t(targetIndex)), 'Color', eyeColor);
        text(screenXpixels -500, -40, sprintf('Found: %.2f ms', timeToTarget), 'Color', 'black');
        %return;  % Exit function after plotting first time target is found
    end

    % Determine the last valid data point time
    lastValidIndex = find(~isnan(x) & ~isnan(y), 1, 'last');
    if lastValidIndex
        timeToEnd = t(lastValidIndex) - t(1);  % Total time from start to last data point in ms
        plot(x(lastValidIndex), y(lastValidIndex), 'o', 'MarkerFaceColor', eyeColor);
        %text(x(lastValidIndex), y(lastValidIndex) + 20, sprintf('Last Time: %.2f ms', t(lastValidIndex)), 'Color', eyeColor);
        text(0, -40, sprintf('Total time: %.2f ms', timeToEnd), 'Color', 'black');
        %return;
    end
end