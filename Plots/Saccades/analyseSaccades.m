function analyseSaccades(data, screenXpixels, screenYpixels, safe, comparison_results_folder)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Saccades


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Saccades (try 2)
% size of stimulation screen for plotting
%% TODO make loop instead of placeholder
trial = 1;

% 1. Preprocessing the Eye-Tracking Data
% Interpolate NaN values for left eye
x_left = data(1).eyeTrial(trial).left.gazePoint.onDisplayArea(1, :) * screenXpixels;
y_left = screenYpixels -  data(1).eyeTrial(trial).left.gazePoint.onDisplayArea(2, :) * screenYpixels; % Inverting y-axis
x_left_interpolated = fillmissing(x_left, 'linear');
y_left_interpolated = fillmissing(y_left, 'linear');

% Interpolate NaN values for right eye
x_right = data(1).eyeTrial(trial).right.gazePoint.onDisplayArea(1, :) * screenXpixels;
y_right = screenYpixels -  data(1).eyeTrial(trial).right.gazePoint.onDisplayArea(2, :) * screenYpixels; % Inverting y-axis
x_right_interpolated = fillmissing(x_right, 'linear');
y_right_interpolated = fillmissing(y_right, 'linear');

% mean gaze coordinates between the left and right eyes
x_mean = mean([x_left_interpolated; x_right_interpolated], 1);
y_mean = mean([y_left_interpolated; y_right_interpolated], 1);

% interpolate to double the data
%timestamps = double(data(1).eyeTrial(1).systemTimeStamp) / 1000000; % Convert to milliseconds
timestamps = double(data(1).eyeTrial(trial).systemTimeStamp) / 1000;
% Generate new timestamps that are halfway between each of the original timestamps
new_timestamps = linspace(timestamps(1), timestamps(end), 2 * length(timestamps) - 1);
% Interpolate the gaze coordinates to match the new temporal resolution
x_mean_interp = interp1(timestamps, x_mean, new_timestamps, 'linear');
y_mean_interp = interp1(timestamps, y_mean, new_timestamps, 'linear');







%% 2. Saccade detection
% Define a velocity threshold for saccade detection
% For 30 Hz data, a good detection threshold might be around 30°/s to 50°/s.
% Saccades at this sampling rate will have larger apparent velocities due to fewer data points representing each movement.
threshold_vel = 100;
threshold_dist = 30;


% --- smooth vs non-smoothed data (unused) --- 
% timestamps = double(data(1).eyeTrial(1).systemTimeStamp) / 1000;
% step_size = 8; threshold_vel = 100; threshold_dist = 30;
% % --- unsmoothed data ---
% [x_mean_interp, y_mean_interp, new_timestamps] = interpolateGazeData(x_mean, y_mean, timestamps, step_size);
% % Original data saccades detection
% [saccade_onsets_dist, saccade_offsets_dist] = detectSaccades(x_mean_interp, y_mean_interp, new_timestamps, 'distance', threshold_dist);
% [saccade_onsets_vel, saccade_offsets_vel] = detectSaccades(x_mean_interp, y_mean_interp, new_timestamps, 'velocity', threshold_vel);
% % --- smoothed data ---
% window_size = 5;  % Define the Gaussian window size
% [x_mean_smoothed, y_mean_smoothed] = smoothGazeData(x_mean, y_mean, window_size);
% % smoothed and interpolated
% [x_mean_smooth_interp, y_smooth_mean_interp, new_smooth_timestamps] = interpolateGazeData(x_mean_smoothed, y_mean_smoothed, timestamps, step_size);
% % Smoothed data saccades detection
% [saccade_onsets_dist_smoothed, saccade_offsets_dist_smoothed] = detectSaccades(x_mean_smooth_interp, y_smooth_mean_interp, new_smooth_timestamps, 'distance', threshold_dist);
% [saccade_onsets_vel_smoothed, saccade_offsets_vel_smoothed] = detectSaccades(x_mean_smooth_interp, y_smooth_mean_interp, new_smooth_timestamps, 'velocity', threshold_vel);
% figure;
% subplot(2, 2, 1); % Distance-Based Saccades (Original Data) 
% plotSaccades(x_mean_interp, y_mean_interp, saccade_onsets_dist, saccade_offsets_dist, sprintf('Original Data: Distance-Based Saccades (Threshold: %d pixels)', threshold_dist));
% subplot(2, 2, 2); % Velocity-Based Saccades (Original Data)
% plotSaccades(x_mean_interp, y_mean_interp, saccade_onsets_vel, saccade_offsets_vel, sprintf('Original Data: Velocity-Based Saccades (Threshold: %d deg/s)', threshold_vel));
% subplot(2, 2, 3); % Distance-Based Saccades (Smoothed Data)
% plotSaccades(x_mean_smooth_interp, y_smooth_mean_interp, saccade_onsets_dist_smoothed, saccade_offsets_dist_smoothed, sprintf('Smoothed Data: Distance-Based Saccades (Threshold: %d pixels)', threshold_dist));
% subplot(2, 2, 4); % Velocity-Based Saccades (Smoothed Data)
% plotSaccades(x_mean_smooth_interp, y_smooth_mean_interp, saccade_onsets_vel_smoothed, saccade_offsets_vel_smoothed, sprintf('Smoothed Data: Velocity-Based Saccades (Threshold: %d deg/s)', threshold_vel));
% sgtitle('Comparison of Saccade Detection Methods for Original and Smoothed Data'); saveas(gcf, 'saccade_comparison_smoothed_plot.png');  % Save the figure as a PNG file


%% --- distance and velocity based saccades --- 

% --- Distances and velocities for non-interpolated data ---
dx = diff(x_mean);  % Difference in x-coordinates (non-interpolated)
dy = diff(y_mean);  % Difference in y-coordinates (non-interpolated)
distance_non_interp = sqrt(dx.^2 + dy.^2);  % Euclidean distance

dt_non_interp = diff(timestamps) / 1000;  % Convert time differences to seconds
velocity_non_interp = distance_non_interp ./ dt_non_interp;  % Velocity

% --- Distances and velocities for interpolated data ---
dx_interp = diff(x_mean_interp);  % Difference in x-coordinates (interpolated)
dy_interp = diff(y_mean_interp);  % Difference in y-coordinates (interpolated)
distance_interp = sqrt(dx_interp.^2 + dy_interp.^2);  % Euclidean distance

dt_interp = diff(new_timestamps) / 1000;  % Convert time differences to seconds
velocity_interp = distance_interp ./ dt_interp;  % Velocity



% ---  Create histograms and log y-axis ---
figure;
% Subplot 1: Histogram of Distances (Non-Interpolated Data)
subplot(2, 2, 1);
histogram(distance_non_interp, 'Normalization', 'probability');
xlabel('Distance (pixels)');
ylabel('Probability');
title('Non-Interpolated Data: Distance');
set(gca, 'YScale', 'log');  % Set y-axis to logarithmic scale
% Subplot 2: Histogram of Velocities (Non-Interpolated Data)
subplot(2, 2, 2);
histogram(velocity_non_interp, 'Normalization', 'probability');
xlabel('Velocity (pixels/second)');
ylabel('Probability');
title('Non-Interpolated Data: Velocity');
set(gca, 'YScale', 'log');  % Set y-axis to logarithmic scale
% Subplot 3: Histogram of Distances (Interpolated Data)
subplot(2, 2, 3);
histogram(distance_interp, 'Normalization', 'probability');
xlabel('Distance (pixels)');
ylabel('Probability');
title('Interpolated Data: Distance');
set(gca, 'YScale', 'log');  % Set y-axis to logarithmic scale
% Subplot 4: Histogram of Velocities (Interpolated Data)
subplot(2, 2, 4);
histogram(velocity_interp, 'Normalization', 'probability');
xlabel('Velocity (pixels/second)');
ylabel('Probability');
title('Interpolated Data: Velocity');
set(gca, 'YScale', 'log');  % Set y-axis to logarithmic scale
sgtitle('Histograms of Distances and Velocities');
saveas(gcf, 'histogram_comparison_plot_log_y.png');  % Save the figure as a PNG file



%% get threshholds for distances/velocities
% --- Calculate the KDE/global minimum thresholds ---
% distributions are bimodal in both cases (slow non saccades vs fast saccades)
% find global minimum and use it as threshhold for detecting saccades
smoothness = 0.4;
smoothness_dist = 0.1;
% Non-interpolated data - Distance
% Fit a Kernel Density Estimate (KDE)
[dni_f, dni_xi] = ksdensity(distance_non_interp, 'Bandwidth', smoothness_dist); % Adjust bandwidth for smoothness
% Find the global minimum in the KDE
[dni_global_min_value, dni_idx] = min(dni_f); % Find the minimum value of f and its index
threshold_dist_non_interp = dni_xi(dni_idx);  % The corresponding x value (threshold)

% Non-interpolated data - Velocity
[vni_f, vni_xi] = ksdensity(velocity_non_interp, 'Bandwidth', smoothness); 
[vni_global_min_value, vni_idx] = min(vni_f);  % Find the global minimum in the KDE
threshold_vel_non_interp = vni_xi(vni_idx);  % The corresponding x value (threshold)

% Interpolated data - Distance
[di_f, di_xi] = ksdensity(distance_interp, 'Bandwidth', smoothness_dist); 
[di_global_min_value, di_idx] = min(di_f);  % Find the global minimum in the KDE
threshold_dist_interp = di_xi(di_idx);  % The corresponding x value (threshold)

% Interpolated data - Velocity
[vi_f, vi_xi] = ksdensity(velocity_interp, 'Bandwidth', smoothness); 
[vi_global_min_value, vi_idx] = min(vi_f);  % Find the global minimum in the KDE
threshold_vel_interp = vi_xi(vi_idx);  % The corresponding x value (threshold)

fprintf('Non-Interpolated Data: Global Minimum Thresholds:\nDistance: %.2f pixels\nVelocity: %.2f pixels/second\n', threshold_dist_non_interp, threshold_vel_non_interp);
fprintf('Interpolated Data: Global Minimum Thresholds:\nDistance: %.2f pixels\nVelocity: %.2f pixels/second\n', threshold_dist_interp, threshold_vel_interp);

% 2x2 plot: KDEs and the global minimum
plotKDEWithGlobalMin(dni_f, dni_xi, threshold_dist_non_interp, dni_global_min_value, ...
                     vni_f, vni_xi, threshold_vel_non_interp, vni_global_min_value, ...
                     di_f, di_xi, threshold_dist_interp, di_global_min_value, ...
                     vi_f, vi_xi, threshold_vel_interp, vi_global_min_value, ...
                     safe, comparison_results_folder, trial);
% TODO change folder above to participant folder or sth 


% --- Detect saccades using calculated thresholds ---
% Non-interpolated data
[saccade_onsets_dist, saccade_offsets_dist] = detectSaccades(dx, dy, timestamps, 'distance', threshold_dist_non_interp);
[saccade_onsets_vel, saccade_offsets_vel] = detectSaccades(dx, dy, timestamps, 'velocity', threshold_vel_non_interp );

% Interpolated data
[saccade_onsets_dist_interp, saccade_offsets_dist_interp] = detectSaccades(dx_interp, dy_interp, new_timestamps, 'distance', threshold_dist_interp);
[saccade_onsets_vel_interp, saccade_offsets_vel_interp] = detectSaccades(dx_interp, dy_interp, new_timestamps, 'velocity', threshold_vel_interp);

% --- Step 6: Plot saccades detected for both datasets ---
figure;
% Subplot 1: Distance-Based Saccades (Non-Interpolated Data)
subplot(2, 2, 1);
plotSaccades(x_mean, y_mean, saccade_onsets_dist, saccade_offsets_dist, sprintf('Non-Interpolated Data: Distance-Based Saccades (Threshold: %.2f pixels)', threshold_dist));
% Subplot 2: Velocity-Based Saccades (Non-Interpolated Data)
subplot(2, 2, 2);
plotSaccades(x_mean, y_mean, saccade_onsets_vel, saccade_offsets_vel, sprintf('Non-Interpolated Data: Velocity-Based Saccades (Threshold: %.2f deg/s)', threshold_vel));
% Subplot 3: Distance-Based Saccades (Interpolated Data)
subplot(2, 2, 3);
plotSaccades(x_mean_interp, y_mean_interp, saccade_onsets_dist_interp, saccade_offsets_dist_interp, sprintf('Interpolated Data: Distance-Based Saccades (Threshold: %.2f pixels)', threshold_dist_interp));
% Subplot 4: Velocity-Based Saccades (Interpolated Data)
subplot(2, 2, 4);
plotSaccades(x_mean_interp, y_mean_interp, saccade_onsets_vel_interp, saccade_offsets_vel_interp, sprintf('Interpolated Data: Velocity-Based Saccades (Threshold: %.2f deg/s)', threshold_vel_interp));
sgtitle('Comparison of Saccade Detection Methods (Interpolated vs Non-Interpolated Data)');
saveas(gcf, 'saccade_comparison_plot.png');  % Save the figure as a PNG file



%% TODO OR just take the min from the histogram instead of the smoothed histogram/KDE



%% --- OR: Calculate the 95th percentile thresholds ---
% Interpolated data
%perc = 85;
%threshold_dist_interp = prctile(distance_interp, perc);
%threshold_vel_interp = prctile(velocity_interp, perc);
threshold_dist_interp = 200; threshold_vel_interp  = 3000;
fprintf('Interpolated Data: %.2ft h Percentile Thresholds:\nDistance: %.2f pixels\nVelocity: %.2f pixels/second\n', perc, threshold_dist_interp, threshold_vel_interp);

% Non-interpolated data
%threshold_dist = prctile(distance_non_interp, perc);
%threshold_vel = prctile(velocity_non_interp, perc);
threshold_dist = 400; threshold_vel  = 3000;
fprintf('Non-Interpolated Data: %.2f th Percentile Thresholds:\nDistance: %.2f pixels\nVelocity: %.2f pixels/second\n', perc, threshold_dist, threshold_vel);

% --- Step 5: Detect saccades using calculated thresholds ---
% Non-interpolated data
[saccade_onsets_dist, saccade_offsets_dist] = detectSaccades(x_mean, y_mean, timestamps, 'distance', threshold_dist);
[saccade_onsets_vel, saccade_offsets_vel] = detectSaccades(x_mean, y_mean, timestamps, 'velocity', threshold_vel);

% Interpolated data
[saccade_onsets_dist_interp, saccade_offsets_dist_interp] = detectSaccades(x_mean_interp, y_mean_interp, new_timestamps, 'distance', threshold_dist_interp);
[saccade_onsets_vel_interp, saccade_offsets_vel_interp] = detectSaccades(x_mean_interp, y_mean_interp, new_timestamps, 'velocity', threshold_vel_interp);

% --- Step 6: Plot saccades detected for both datasets ---
figure;
% Subplot 1: Distance-Based Saccades (Non-Interpolated Data)
subplot(2, 2, 1);
plotSaccades(x_mean, y_mean, saccade_onsets_dist, saccade_offsets_dist, sprintf('Non-Interpolated Data: Distance-Based Saccades (Threshold: %.2f pixels)', threshold_dist));
% Subplot 2: Velocity-Based Saccades (Non-Interpolated Data)
subplot(2, 2, 2);
plotSaccades(x_mean, y_mean, saccade_onsets_vel, saccade_offsets_vel, sprintf('Non-Interpolated Data: Velocity-Based Saccades (Threshold: %.2f deg/s)', threshold_vel));
% Subplot 3: Distance-Based Saccades (Interpolated Data)
subplot(2, 2, 3);
plotSaccades(x_mean_interp, y_mean_interp, saccade_onsets_dist_interp, saccade_offsets_dist_interp, sprintf('Interpolated Data: Distance-Based Saccades (Threshold: %.2f pixels)', threshold_dist_interp));
% Subplot 4: Velocity-Based Saccades (Interpolated Data)
subplot(2, 2, 4);
plotSaccades(x_mean_interp, y_mean_interp, saccade_onsets_vel_interp, saccade_offsets_vel_interp, sprintf('Interpolated Data: Velocity-Based Saccades (Threshold: %.2f deg/s)', threshold_vel_interp));
sgtitle('Comparison of Saccade Detection Methods (Interpolated vs Non-Interpolated Data)');
saveas(gcf, 'saccade_comparison_plot.png');  % Save the figure as a PNG file





%% interpolated saccade data
% Example timestamps vector (your actual data will be different)
timestamps = double(data(1).eyeTrial(1).systemTimeStamp) / 1000;  % Convert to milliseconds
% Check the range of your timestamps
start_time = timestamps(1);
end_time = timestamps(end);
time_range = end_time - start_time;
% Decide on step size based on range
if time_range < 8
   % step_size = 1;  % Use 1 ms steps if the range is too small
   disp('not enough data/too short interval/trial');
else
    step_size = 8;  % Otherwise, use 8 ms steps
end
% Generate new timestamps for interpolation
new_timestamps = start_time:step_size:end_time;
% Display the size and contents of new_timestamps to confirm
disp(['new_timestamps has ', num2str(length(new_timestamps)), ' elements.']);
% Interpolate the mean gaze coordinates to match the new temporal resolution
x_mean_interp = interp1(timestamps, x_mean, new_timestamps, 'linear');
y_mean_interp = interp1(timestamps, y_mean, new_timestamps, 'linear');
% plot interpolated gaze data
figure; hold on;
plot(x_mean, y_mean, 'bo-', 'MarkerSize', 3, 'DisplayName', 'Original Gaze Path');
plot(x_mean_interp, y_mean_interp, 'r.-', 'MarkerSize', 10, 'DisplayName', 'Interpolated Gaze Path');
xlabel('Horizontal Gaze Position (pixels)');
ylabel('Vertical Gaze Position (pixels)');
title('Interpolated Gaze Path Visualization');
legend('Original', 'Interpolated');
grid on; hold off;
