function [result_table, tnf] = findTargetAndExtractData(cutData, screenXpixels, screenYpixels, analysis_folder, tolerance, num_before, num_after)
    result_table = {}; 
    tnf = 0;  % Counter for trials where the target not found

    % Loop through each trial in the cutData
    for trial = 1:size(cutData, 1)
        trial_data = cutData(trial, :);
        condition = trial_data.Condition{1};  % Extract the condition for the trial
        % get the pupil diameter during blank screen for baseline
        blank_start = cutData.blankStartTime(trial);
        timestamps = double(trial_data.eyeTrial.systemTimeStamp) / 1e6;
        [~, idx_blank] = min(abs(timestamps - blank_start)); % get start of blank screen
        blank_d_l = trial_data.eyeTrial.left.pupil.diameter;
        blank_d_r = trial_data.eyeTrial.right.pupil.diameter; 
        blank_d_avg = nanmean([blank_d_l; blank_d_r], 1); % avg over both eyes
        blank_diam_mean = nanmean(blank_d_avg(idx_blank:idx_blank+9)); % avg over time of blank screen (10 datapoints)
        blank_diam_median = nanmedian(blank_d_avg(idx_blank:idx_blank+9)); % avg over time of blank screen (10 datapoints)

        % Extract target position    
        targetPos = trial_data.TargetPosition;
        targetRow = targetPos(2); % Row index of the target
        targetCol = targetPos(1); % Column index of the target

        % Get the target X and Y coordinates
        targetX = trial_data.x_centers{1}(targetRow, targetCol);
        targetY = screenYpixels - trial_data.y_centers{1}(targetRow, targetCol);

        % Right eye data
        a_r = trial_data.stimulusTrial.right.gazePoint.onDisplayArea;
        x_r = a_r(1,:) * screenXpixels;
        y_r = screenYpixels - (a_r(2,:) * screenYpixels);  % Inverting y-axis
        t_r = double(trial_data.stimulusTrial.systemTimeStamp) / 1e6;  % Converts timestamp to seconds
        d_r = trial_data.stimulusTrial.right.pupil.diameter; % pupil diameter
        % Left eye data
        a_l = trial_data.stimulusTrial.left.gazePoint.onDisplayArea;
        x_l = a_l(1,:) * screenXpixels;
        y_l = screenYpixels - (a_l(2,:) * screenYpixels);  % Inverting y-axis
        t_l = double(trial_data.stimulusTrial.systemTimeStamp) / 1e6;  % Converts timestamp to seconds
        d_l = trial_data.stimulusTrial.left.pupil.diameter; % pupil diameter
        
        % Calculate the average gaze position between both eyes
        x_avg = nanmean([x_r; x_l], 1);
        y_avg = nanmean([y_r; y_l], 1);
        d_avg = nanmean([d_r; d_l], 1);

        % Determine the index where the eyes first "find" the target
        found_target_idx = find(sqrt((x_avg - targetX).^2 + (y_avg - targetY).^2) <= tolerance, 1, 'first');

        % Check if the target was found
        if isempty(found_target_idx)
            fprintf('Target not found in trial %d. Skipping this trial.\n', trial);
            tnf = tnf + 1;
            continue;  % Skip to the next trial
        end

        % Extract the last 30 points before finding the target
        if found_target_idx >= num_before
            points_before_stim = d_avg(found_target_idx - num_before + 1 : found_target_idx);
        else
            % If there are fewer than `num_before` points, pad with NaNs
            points_before_stim = [NaN(1, num_before - found_target_idx), d_avg(1:found_target_idx)];
        end

        % Extract the 10 points after finding the target
        if found_target_idx + num_after <= length(x_avg)
            points_after_stim = d_avg(found_target_idx + 1 : found_target_idx + num_after);
        else
            % If there are fewer than `num_after` points, pad with NaNs
            points_after_stim = [d_avg(found_target_idx + 1 : end), NaN(1, num_after - (length(x_avg) - found_target_idx))];
        end

        % Ensure that both before and after arrays are the correct length
        if length(points_before_stim) > num_before
            points_before_stim = points_before_stim(1:num_before);
        end
        if length(points_after_stim) > num_after
            points_after_stim = points_after_stim(1:num_after);
        end

        % Save the trial number, condition, and datapoints in a row
        result_table = [result_table; {trial, condition, points_before_stim, points_after_stim, blank_diam_mean, blank_diam_median}];
    end

    % Convert the cell array to a table for easier processing
    result_table = cell2table(result_table, 'VariableNames', {'TrialNumber', 'Condition', 'DataPointsBefore', 'DataPointsAfter', 'MeanDiamBlank', 'MedianDiamBlank'});

    % Save the table to a file
    save(fullfile(analysis_folder, '\pupil_before_after_finding_stim.mat'), 'result_table');
    disp(tnf);
    disp('Trials where the target was not found with the eyes.');
end
