function trial_metrics = calcGazeDistanceDirection(data, all_fixations, all_saccades, tolerance, screenXpixels, screenYpixels)

trial_metrics = struct();
for participant = 1:size(data, 2)
    trial_metrics(participant).id = data(participant).id;
    start_at_target = 0;
    for trial = 1:size(data(participant).stimulusTrial, 1)
        participant_data = data(participant).stimulusTrial(trial);
        participant_fixations = all_fixations(participant).fixations.stimulusFixations(trial).fixations;
        participant_saccades = all_saccades(participant).saccades(trial).saccades;
        
        % Coordinates of the first gaze of the trial
        first_gaze_x = participant_data.left.gazePoint.onDisplayArea(1, 1) * screenXpixels;
        first_gaze_y = screenYpixels - (participant_data.left.gazePoint.onDisplayArea(2, 1) * screenYpixels);
        
        % Coordinates of the target
        target_x = all_fixations(participant).fixations.targetCenters(1, trial);
        target_y = all_fixations(participant).fixations.targetCenters(2, trial);
        
        % Optimal path between first gaze and target (Euclidean distance)
        optimal_distance = sqrt((target_x - first_gaze_x)^2 + (target_y - first_gaze_y)^2);
        % Direction of the optimal path (angle in radians)
        optimal_angle = atan2(target_x - first_gaze_x, target_y - first_gaze_y);   
        
        % Get first saccade data (if available)
        if ~isempty(participant_saccades)
            % Coordinates of the first saccade start and end
            first_saccade_start_x = participant_saccades(1).startCenter(1,1) * screenXpixels;
            first_saccade_start_y = screenYpixels - (participant_saccades(1).startCenter(1,2) * screenYpixels);
            first_saccade_end_x = participant_saccades(1).endCenter(1,1) * screenXpixels;
            first_saccade_end_y = screenYpixels - (participant_saccades(1).endCenter(1,2) * screenYpixels);
            
            % First saccade path (Euclidean distance)
            first_saccade_distance = participant_saccades(1).distance;
                  
            % Direction of the first saccade (angle in radians)
            first_saccade_angle = atan2(first_saccade_end_x - first_saccade_start_x, ...
                                        first_saccade_end_y - first_saccade_start_y);
                                        
            % Difference between optimal path and first saccade direction (in degrees)
            direction_difference = rad2deg(abs(optimal_angle - first_saccade_angle));
            
            % time until first saccade starts
            first_saccade_start = participant_saccades(1).saccStartAfterStimOnset;
        
            % Check if the saccade direction was optimal (i.e., minimal difference)
            if direction_difference <= 10  % Example threshold of 10 degrees
                direction_good = 1;
            else
                direction_good = 0;
            end
            
            
        else
            %disp('No saccades in this trial.');
        end
        
        % Number of saccades in trial
        n_sac = size(participant_saccades, 2);
        % Total distance of saccades (sum of all saccade distances)
        total_dist_sac = sum([participant_saccades.distance]);
        % Mean saccade distance
        mean_dist_sac = mean([participant_saccades.distance]);
        
        % Find the number of saccades until the target is found
        gaze_found_target_idx = find(sqrt((participant_data.left.gazePoint.onDisplayArea(1,:) * screenXpixels - target_x).^2 + ...
                                          (screenYpixels - (participant_data.left.gazePoint.onDisplayArea(2,:) * screenYpixels) - target_y).^2) <= tolerance, 1);
                   
        if ~isempty(gaze_found_target_idx)
            % Number of saccades until the target is found
            saccades_until_target = find([participant_saccades.startIdx] <= gaze_found_target_idx, 1, 'last');
            
            % Number of saccades after target is found
            if isempty(saccades_until_target) % gaze started at target
                saccades_until_target = 0;
                start_at_target = start_at_target +1;
            else
                saccades_after_target = n_sac - saccades_until_target;
                % time of saccade where subject finds target 
                relevant_saccade_start = participant_saccades(saccades_until_target).saccStartAfterStimOnset;
                relevant_saccade_end   = participant_saccades(saccades_until_target).saccStartAfterStimOnset;
            end
        else
            saccades_until_target  = NaN;  % Target not found
            saccades_after_target  = NaN;  % No data for after target found
            relevant_saccade_start = NaN; 
            relevant_saccade_end   = NaN; 
        end
        

        % Store the calculated metrics for each trial
        trial_metrics(participant).optimal_distance(trial)       = optimal_distance;
        trial_metrics(participant).optimal_angle(trial)          = optimal_angle;
        
        trial_metrics(participant).first_saccade_distance(trial) = first_saccade_distance;
        trial_metrics(participant).first_saccade_angle(trial)    = first_saccade_angle;
        
        % first saccade timepoint start
        trial_metrics(participant).first_saccade_start(trial)    = first_saccade_start;
        trial_metrics(participant).relevant_saccade_start(trial) = relevant_saccade_start;
        trial_metrics(participant).relevant_saccade_end(trial)   = relevant_saccade_end;
        
        trial_metrics(participant).direction_difference(trial)   = direction_difference;
        trial_metrics(participant).direction_good(trial)         = direction_good;
        trial_metrics(participant).total_dist_sac(trial)         = total_dist_sac;
        trial_metrics(participant).n_sac(trial)                  = n_sac;
        trial_metrics(participant).mean_dist_sac(trial)          = mean_dist_sac;
        trial_metrics(participant).saccades_until_target(trial)  = saccades_until_target;
        trial_metrics(participant).saccades_after_target(trial)  = saccades_after_target;
    end
    trial_metrics(participant).start_at_target = start_at_target;
end
