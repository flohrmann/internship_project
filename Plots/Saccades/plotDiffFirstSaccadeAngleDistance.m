function plotDiffFirstSaccadeAngleDistance(trial_metrics, compare_folder, safe)


first_saccade_angles = trial_metrics.first_saccade_angle;       % Actual first saccade angles
optimal_angles = trial_metrics.optimal_angle;                   % Optimal angles
first_saccade_distances = trial_metrics.first_saccade_distance; % Actual first saccade distances
optimal_distances = trial_metrics.optimal_distance;             % Optimal distances
num_trials = length(first_saccade_angles);

figure;
% Angle comparison
subplot(2, 1, 1);
plot(1:num_trials, rad2deg(optimal_angles), '-b', 'DisplayName', 'Optimal Angle'); hold on;
plot(1:num_trials, rad2deg(first_saccade_angles), '-r', 'DisplayName', 'First Saccade Angle');
plot(1:num_trials, abs(rad2deg(optimal_angles - first_saccade_angles)), '--k', 'DisplayName', 'Angle Difference');
xlabel('Trial Number'); ylabel('Angle (degrees)'); title('Angle Comparison');
legend('show', 'Location', 'best'); grid on; hold off;

% Distance comparison
subplot(2, 1, 2);
plot(1:num_trials, optimal_distances, '-b', 'DisplayName', 'Optimal Distance'); hold on;
plot(1:num_trials, first_saccade_distances, '-r', 'DisplayName', 'First Saccade Distance');
plot(1:num_trials, abs(optimal_distances - first_saccade_distances), '--k', 'DisplayName', 'Distance Difference');
xlabel('Trial Number'); ylabel('Distance (pixels)'); title('Distance Comparison');
legend('show', 'Location', 'best');grid on; hold off;

if safe == 1
    print(gcf, fullfile(compare_folder,  strcat('saccade_dist_angle_per_trial_id_',num2str(trial_metrics.id),'.svg')), '-dsvg');
end