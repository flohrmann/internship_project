function trial_metrics = calcGazeNonTargetSaccadesDistanceDirection(cutData, fixations, saccades, tolerance, num_before, num_after, baseline_length, screenXpixels, screenYpixels)

    % Initialize trial_metrics structure
    trial_metrics = struct();
    trial_metrics.id = fixations.id;
    
    % For each trial, calculate metrics
    for trial = 1:size(cutData.stimulusTrial, 1)
        trial_data = cutData(trial, :);
        current_saccades = saccades(trial).saccades;
        
        trial_metrics.trial(trial, :)       = trial;

        % Baseline pupil diameters
        blank_start = cutData.blankStartTime(trial);
        timestamps = double(trial_data.eyeTrial.systemTimeStamp) / 1e6;
        [~, idx_blank] = min(abs(timestamps - blank_start));
        
        blank_diam_l = trial_data.eyeTrial.left.pupil.diameter(idx_blank:idx_blank+baseline_length-1);
        blank_diam_r = trial_data.eyeTrial.right.pupil.diameter(idx_blank:idx_blank+baseline_length-1);
        blank_diam_avg = nanmean([blank_diam_l; blank_diam_r],1);
        
        trial_start_diam_l = trial_data.stimulusTrial.left.pupil.diameter(1:baseline_length);
        trial_start_diam_r = trial_data.stimulusTrial.right.pupil.diameter(1:baseline_length);
        trial_start_diam_avg = nanmean([trial_start_diam_l; trial_start_diam_r],1);
        
            % Coordinates of the first gaze of the trial
        first_gaze_x = trial_data.stimulusTrial.left.gazePoint.onDisplayArea(1, 1) * screenXpixels;
        first_gaze_y = screenYpixels - (trial_data.stimulusTrial.left.gazePoint.onDisplayArea(2, 1) * screenYpixels);
    
    
        % Target coordinates
        target_x = fixations.targetCenters(1, trial);
        target_y = fixations.targetCenters(2, trial);

        % Saccade properties for non-target saccades
        nts_distances = [];
        nts_angles = [];
        nts_direction_diffs = [];
        nts_start_times = [];
        nts_diam_before = [];
        nts_diam_after = [];
        nts_start = [];
        nts_end = [];
        
        % Iterate through each saccade
        for s = 1:length(current_saccades)
            % Coordinates of the current saccade start and end
            saccade_start_x = current_saccades(s).startCenter(1);
            saccade_start_y = current_saccades(s).startCenter(2);
            saccade_end_x   = current_saccades(s).endCenter(1);
            saccade_end_y   = current_saccades(s).endCenter(2);
            % Distance to target and direction difference
            distance_to_target = sqrt((saccade_end_x - target_x)^2 + (saccade_end_y - target_y)^2);
            saccade_angle = atan2(saccade_end_x - saccade_start_x, saccade_end_y - saccade_start_y);
            optimal_angle = atan2(target_x - saccade_start_x, target_y - saccade_start_y);
            direction_difference = rad2deg(abs(optimal_angle - saccade_angle));

            % Check if saccade does not go toward the target
            if distance_to_target > tolerance
                % Store metrics for non-target saccade
                nts_distances = [nts_distances; current_saccades(s).distance];
                nts_angles = [nts_angles; saccade_angle];
                nts_direction_diffs = [nts_direction_diffs; direction_difference];
                nts_start_times = [nts_start_times; current_saccades(s).saccStartAfterStimOnset];
                
                % Pupil data before and after non-target saccade
                nts_start = [nts_start; current_saccades(s).startIdx];
                nts_end = [nts_end; current_saccades(s).endIdx];
                nts_diam_before = [nts_diam_before; get_pupil_data(trial_data, current_saccades(s).startIdx, num_before, 'before')];
                nts_diam_after = [nts_diam_after; get_pupil_data(trial_data, current_saccades(s).endIdx, num_after, 'after')];
            end
        end

        % Store all non-target saccade metrics in trial_metrics struct
        trial_metrics.optimal_distance(trial,:)      = sqrt((target_x - first_gaze_x)^2 + (target_y - first_gaze_y)^2);
        trial_metrics.blank_diam_avg(trial, :)       = blank_diam_avg;
        trial_metrics.trial_start_diam_avg(trial, :) = trial_start_diam_avg;

        trial_metrics.nts_distances{trial,:}         = nts_distances;
        trial_metrics.nts_angles{trial,:}            = nts_angles;
        trial_metrics.nts_direction_diffs{trial,:}   = nts_direction_diffs;
        trial_metrics.nts_start_times{trial,:}       = nts_start_times;
        trial_metrics.nts_start{trial,:}             = nts_start;       
        trial_metrics.nts_diam_before{trial,:}       = nts_diam_before;
        trial_metrics.nts_diam_after{trial,:}        = nts_diam_after;
    end
end

% Helper function to extract pupil data around a saccade
function data_segment = get_pupil_data(trial_data, saccade_idx, num_points, direction)
    d_avg = nanmean([trial_data.stimulusTrial.right.pupil.diameter; trial_data.stimulusTrial.left.pupil.diameter], 1);
    if strcmp(direction, 'before')
        if saccade_idx > num_points
            data_segment = d_avg(saccade_idx - num_points:saccade_idx - 1);
        else
            data_segment = [NaN(1, num_points - saccade_idx+1), d_avg(1:saccade_idx - 1)];
        end
    else
        if saccade_idx + num_points <= length(d_avg)
            data_segment = d_avg(saccade_idx:saccade_idx + num_points - 1);
        else
            data_segment = [d_avg(saccade_idx:end), NaN(1, num_points - (length(d_avg) - saccade_idx + 1))];
        end
    end
end
