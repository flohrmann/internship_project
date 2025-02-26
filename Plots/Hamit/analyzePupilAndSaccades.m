function [epoch_data, confirmatory_epoch_data] = analyzePupilAndSaccades(data, participant_saccades, screenXpixels, screenYpixels, baseline_window, post_saccade_window, gaze_shift_threshold)
    % Function to analyze pupil size around saccades and confirmatory analysis
    % Inputs:
    %   - data: data structure containing stimulusTrial information
    %   - participant_saccades: structure containing saccade information
    %   - screenXpixels, screenYpixels: screen dimensions
    %   - baseline_window: time (in ms) before saccade for baseline calculation
    %   - post_saccade_window: time (in ms) after saccade for pupil analysis
    %   - gaze_shift_threshold: maximum allowed gaze shift distance (in pixels) for confirmatory analysis
    %
    % Outputs:
    %   - epoch_data: structure containing pupil data epochs around saccades
    %   - confirmatory_epoch_data: structure containing confirmatory epochs for trials with minimal gaze shifts

    % Initialize structures to store epoch data and confirmatory epoch data
    epoch_data = struct();
    confirmatory_epoch_data = struct();

    sampling_rate = 60;  % Hz, adjust this if your data has a different sampling rate

    % Loop through participants
    for participant = 1:size(data, 2)
        epoch_data(participant).id = data(participant).id;
        confirmatory_epoch_data(participant).id = data(participant).id;

        % Loop through each trial for the current participant
        for trialIdx = 1:size(data(participant).stimulusTrial, 1)
            % Extract pupil and saccade data
            pupil_data = data(participant).stimulusTrial(trialIdx).left.pupilDiameter;  % Pupil diameter data
            pupil_timestamps = data(participant).stimulusTrial(trialIdx).timestamps;  % Timestamps
            saccades = participant_saccades(participant).saccades(trialIdx).saccades;  % Saccade data
            
            % Skip if no saccades in this trial
            if isempty(saccades)
                continue;
            end

            % Coordinates of gaze (to compute distances)
            xGaze = (data(participant).stimulusTrial(trialIdx).left.gazePoint.onDisplayArea(1, :)) * screenXpixels;
            yGaze = screenYpixels - (data(participant).stimulusTrial(trialIdx).left.gazePoint.onDisplayArea(2, :)) * screenYpixels;

            % Target coordinates
            target_x = participant_saccades(participant).targetCenters(1, trialIdx);
            target_y = participant_saccades(participant).targetCenters(2, trialIdx);

            % Compute optimal path (Euclidean distance) between first gaze and target
            first_gaze_x = xGaze(1);  % X coordinate of first gaze point
            first_gaze_y = yGaze(1);  % Y coordinate of first gaze point
            optimal_path = sqrt((target_x - first_gaze_x)^2 + (target_y - first_gaze_y)^2);

            %% Epoch Data Around Saccades
            for saccade_idx = 1:length(saccades)
                % Start time of the current saccade
                saccade_start_time = saccades(saccade_idx).startTime;
                saccade_start_idx = find(pupil_timestamps >= saccade_start_time, 1);  % Find index of saccade start
                
                % Define the baseline and post-saccade windows
                baseline_start_idx = max(1, saccade_start_idx - baseline_window);
                post_saccade_end_idx = min(length(pupil_data), saccade_start_idx + post_saccade_window);
                
                % Extract baseline and post-saccade pupil data
                baseline_pupil = pupil_data(baseline_start_idx:saccade_start_idx-1);
                post_saccade_pupil = pupil_data(saccade_start_idx:post_saccade_end_idx);

                % Compute baseline correction
                baseline_mean = mean(baseline_pupil);
                pupil_change = post_saccade_pupil - baseline_mean;

                % Store the epoch data in the structure
                epoch_data(participant).trial(trialIdx).saccade(saccade_idx).baseline = baseline_pupil;
                epoch_data(participant).trial(trialIdx).saccade(saccade_idx).post_saccade = pupil_change;

                % First saccade path (Euclidean distance)
                first_saccade_distance = saccades(1).distance;
                
                % Direction analysis (difference between optimal and first saccade direction)
                first_saccade_end_x = saccades(1).endCenter(1,1) * screenXpixels;
                first_saccade_end_y = screenYpixels - (saccades(1).endCenter(1,2) * screenYpixels);
                first_saccade_angle = atan2(first_saccade_end_y - first_gaze_y, first_saccade_end_x - first_gaze_x);
                optimal_angle = atan2(target_y - first_gaze_y, target_x - first_gaze_x);
                direction_difference = rad2deg(abs(optimal_angle - first_saccade_angle));
                
                % Save optimal path and first saccade metrics
                epoch_data(participant).trial(trialIdx).optimal_path = optimal_path;
                epoch_data(participant).trial(trialIdx).first_saccade_distance = first_saccade_distance;
                epoch_data(participant).trial(trialIdx).direction_difference = direction_difference;
            end

            %% Confirmatory Analysis: Target Found in First Saccade and Minimal Gaze Shift
            if trial_metrics(participant).saccades_until_target(trialIdx) == 1
                % Calculate gaze shifts between saccades
                gaze_shifts = [saccades.distance];
                
                % Confirm if gaze shifts are below the threshold
                if all(gaze_shifts <= gaze_shift_threshold)
                    % Epoch data time-locked to first saccade
                    first_saccade_time = saccades(1).startTime;
                    first_saccade_idx = find(pupil_timestamps >= first_saccade_time, 1);
                    
                    % Define baseline and post-saccade windows for confirmatory analysis
                    baseline_start_idx = max(1, first_saccade_idx - baseline_window);
                    post_saccade_end_idx = min(length(pupil_data), first_saccade_idx + post_saccade_window);

                    % Extract baseline and post-saccade pupil data
                    baseline_pupil = pupil_data(baseline_start_idx:first_saccade_idx-1);
                    post_saccade_pupil = pupil_data(first_saccade_idx:post_saccade_end_idx);

                    % Perform baseline correction
                    baseline_mean = mean(baseline_pupil);
                    pupil_change = post_saccade_pupil - baseline_mean;

                    % Store confirmatory epoch data
                    confirmatory_epoch_data(participant).trial(trialIdx).first_saccade.baseline = baseline_pupil;
                    confirmatory_epoch_data(participant).trial(trialIdx).first_saccade.post_saccade = pupil_change;
                end
            end
        end
    end

    % Save the data for analysis
    save('epoch_data.mat', 'epoch_data');
    save('confirmatory_epoch_data.mat', 'confirmatory_epoch_data');
end
