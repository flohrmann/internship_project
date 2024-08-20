function [x_smoothed, y_smoothed] = smoothGazeData(x, y, window_size)
    % Smooth the gaze data using a Gaussian kernel
    % Inputs:
    %   x - Gaze data for the x-coordinate
    %   y - Gaze data for the y-coordinate
    %   window_size - Size of the Gaussian window for smoothing
    % Outputs:
    %   x_smoothed - Smoothed x-coordinate data
    %   y_smoothed - Smoothed y-coordinate data

    % Define the Gaussian kernel
    gaussian_kernel = gausswin(window_size) / sum(gausswin(window_size));

    % Apply Gaussian smoothing to the x and y coordinates
    x_smoothed = conv(x, gaussian_kernel, 'same');
    y_smoothed = conv(y, gaussian_kernel, 'same');
end
