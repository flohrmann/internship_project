function [lookedAtFixation] = checkFixation(trial_results, eye_tracking_data, screenXpixels, screenYpixels, fixationThreshold, analysis_folder)
    % Initialize the output array
    lookedAtFixation = false(size(trial_results, 1), 2); % Columns for right and left eyes

    % Define the center of the screen
    centerX = screenXpixels / 2;
    centerY = screenYpixels / 2;

    % Loop through each trial in trial_results
    for trial = 1:size(trial_results, 1)
        % Extract the current trial data
        current_data = trial_results(trial, :);

        % Extract the relevant times
        fixStartTime = current_data.fixStartTime;
        blankStartTime = current_data.blankStartTime;

        % Extract system timestamps and gaze points for both eyes
        timestamps = double(eye_tracking_data.systemTimeStamp) / 1000000;

        % Find indices for the fixation period
        idx_fixation = timestamps >= fixStartTime & timestamps < blankStartTime;

        % Extract gaze points during fixation for both eyes
        gazeX_right = eye_tracking_data.right.gazePoint.onDisplayArea(1, idx_fixation) * screenXpixels;
        gazeY_right = screenYpixels - (eye_tracking_data.right.gazePoint.onDisplayArea(2, idx_fixation) * screenYpixels); % Inverting y-axis

        gazeX_left = eye_tracking_data.left.gazePoint.onDisplayArea(1, idx_fixation) * screenXpixels;
        gazeY_left = screenYpixels - (eye_tracking_data.left.gazePoint.onDisplayArea(2, idx_fixation) * screenYpixels); % Inverting y-axis

        % Check if gaze points are within the fixation threshold for the right eye
        distance_right = sqrt((gazeX_right - centerX).^2 + (gazeY_right - centerY).^2);
        if any(distance_right < fixationThreshold)
            lookedAtFixation(trial, 1) = true;
        end

        % Check if gaze points are within the fixation threshold for the left eye
        distance_left = sqrt((gazeX_left - centerX).^2 + (gazeY_left - centerY).^2);
        if any(distance_left < fixationThreshold)
            lookedAtFixation(trial, 2) = true;
        end
    end


    fixation_summary = any(lookedAtFixation, 2);
    fixationSummary = table();
    fixationSummary.fix_count = sum(fixation_summary(:,1));
    fixationSummary.fix_perc  = fixationSummary.fix_count / size(lookedAtFixation,1) * 100;
    fixationSummary.fix_l     = {double(lookedAtFixation(:, 2))};
    fixationSummary.fix_r     = {double(lookedAtFixation(:, 1))};
    fixationSummary.fix_any   = {double(fixation_summary(:))};
    save(fullfile(analysis_folder, 'fixationSummary.mat'), 'fixationSummary');
    % Display the results
%     for trial = 1:size(trial_results, 1)
%         fprintf('Trial %d: Right Eye - %s, Left Eye - %s\n', trial, ...
%             lookedAtFixation(trial, 1) * 'Looked' + ~lookedAtFixation(trial, 1) * 'Did not look', ...
%             lookedAtFixation(trial, 2) * 'Looked' + ~lookedAtFixation(trial, 2) * 'Did not look');
%     end
end
