function plotPupilSizeAroundEvents(cutData)
    % Define the time window around the events (in seconds)
    pre_event_time = 0.2; % Start from StimulusOnsetTime
    post_event_time = 0; % 2 seconds after the event (adjust as needed)
    sampling_rate = 60; % Assuming 60 Hz sampling rate

    % Define colors for plotting
    color_blue = [0, 0.4470, 0.7410]; % Blue
    color_green = [0, 0.5, 0]; % Green

    % Initialize arrays to hold resampled data
    time_base = pre_event_time:(1/sampling_rate):post_event_time; % Common time base
    resampled_pupil_right = nan(size(cutData, 1), length(time_base));
    resampled_pupil_left = nan(size(cutData, 1), length(time_base));

    % Loop through each trial in cutData
    for trial = 1:size(cutData, 1)
        % Extract the current trial data
        current_data = cutData(trial, :);

        % Extract relevant times (in seconds)
        StimulusOnsetTime = current_data.StimulusOnsetTime / 1e6;
        trialEndTime = current_data.trialEndTime / 1e6;

        % Extract system timestamps and pupil diameters for both eyes
        timestamps = double(current_data.eyeTrial.systemTimeStamp) / 1e6;
        pupil_right = current_data.eyeTrial.right.pupil.diameter;
        pupil_left = current_data.eyeTrial.left.pupil.diameter;

        % Resample the data from StimulusOnsetTime to trialEndTime
        idx_segment = timestamps >= StimulusOnsetTime & timestamps <= trialEndTime;
        time_segment = timestamps(idx_segment) - StimulusOnsetTime;
        resampled_pupil_right(trial, :) = interp1(time_segment, pupil_right(idx_segment), time_base, 'linear', NaN);
        resampled_pupil_left(trial, :) = interp1(time_segment, pupil_left(idx_segment), time_base, 'linear', NaN);
    end

    % Calculate the mean pupil diameter across trials
    mean_pupil_right = nanmean(resampled_pupil_right, 1);
    mean_pupil_left = nanmean(resampled_pupil_left, 1);

    % Create a new figure
    figure;
    hold on;
    title('Average Pupil Diameter from Stimulus Onset to Trial End');
    xlabel('Time (s)');
    ylabel('Pupil Diameter (mm)');

    % Plot the mean pupil diameter
    plot(time_base, mean_pupil_right, 'Color', color_blue, 'LineWidth', 2, 'DisplayName', 'Right Eye');
    plot(time_base, mean_pupil_left, 'Color', color_green, 'LineWidth', 2, 'DisplayName', 'Left Eye');

    % Add vertical lines for key events
    xline(0, '--k', 'Stimulus Onset');
    xline(trialEndTime - StimulusOnsetTime, '--r', 'Trial End');

    legend show;
    hold off;
end

