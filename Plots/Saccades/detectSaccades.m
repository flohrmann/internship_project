function [saccade_onsets, saccade_offsets] = detectSaccades(dx, dy, timestamps, method, threshold)
%     % Calculate the difference in x and y coordinates between consecutive points
%     dx = diff(x);
%     dy = diff(y);

    % Calculate the Euclidean distance between consecutive points
    distance = sqrt(dx.^2 + dy.^2);

    if strcmp(method, 'distance')
        % Saccade detection based on distance
        saccades = distance > threshold;
    elseif strcmp(method, 'velocity')
        % Calculate the time differences (dt in milliseconds, convert to seconds)
        dt = diff(timestamps) / 1000;  % Convert milliseconds to seconds

        % Calculate the velocity (distance per second)
        velocity = distance ./ dt;

        % Saccade detection based on velocity
        saccades = velocity > threshold;
    else
        error('Unknown method. Use ''distance'' or ''velocity''.');
    end

    % Ensure saccades is a column vector
    saccades = saccades(:);

    % Identify the onset and offset of saccades
    saccade_onsets = find(diff(saccades) == 1);
    saccade_offsets = find(diff(saccades) == -1);

    % Adjust for the last data point if it is part of a saccade
    if saccades(end) == 1
        saccade_offsets = [saccade_offsets; length(saccades)];
    end
end
