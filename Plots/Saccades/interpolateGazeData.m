function [x_interp, y_interp, new_timestamps] = interpolateGazeData(x, y, timestamps, step_size)
    % Interpolate gaze data to a new temporal resolution
    % Inputs:
    %   x - Gaze data for the x-coordinate
    %   y - Gaze data for the y-coordinate
    %   timestamps - Original timestamps (in milliseconds)
    %   step_size - Step size for new temporal resolution (in milliseconds)
    % Outputs:
    %   x_interp - Interpolated x-coordinate data
    %   y_interp - Interpolated y-coordinate data
    %   new_timestamps - New timestamps for interpolated data

    % Check the range of your timestamps
    start_time = timestamps(1);
    end_time = timestamps(end);
    time_range = end_time - start_time;

    % Ensure the time range is sufficient
    if time_range < step_size
        error('Not enough data or interval too short for interpolation.');
    end

    % Generate new timestamps for interpolation
    new_timestamps = start_time:step_size:end_time;

    % Interpolate the gaze coordinates to match the new temporal resolution
    x_interp = interp1(timestamps, x, new_timestamps, 'linear');
    y_interp = interp1(timestamps, y, new_timestamps, 'linear');
end
