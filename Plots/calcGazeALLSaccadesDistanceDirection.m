function trial_metrics = calcGazeALLSaccadesDistanceDirection(cutData, fixations, saccades, tolerance, num_before, num_after, baseline_length, screenXpixels, screenYpixels, min_length)

trial_metrics = struct();
trial_metrics.id = fixations.id;
start_at_target = 0;

for trial = 1:size(cutData.stimulusTrial, 1)
    trial_data = cutData(trial,:);
    current_fixations = fixations.stimulusFixations(trial).fixations;
    current_saccades = saccades(trial).saccades;
    trial_metrics.trial(trial, :)       = trial;
    %  baselines
    %blank_start = cutData.blankStartTime(trial);
    %timestamps = double(trial_data.eyeTrial.systemTimeStamp) / 1e6;
    %[~, idx_blank] = min(abs(timestamps - blank_start)); % get start of blank screen
    idx_blank = cutData.blankIdx(trial) -cutData.startIdx(trial);
    try
        blank_diam_l = trial_data.eyeTrial.left.pupil.diameter(idx_blank:idx_blank+baseline_length-1);
    catch
        disp('trial with too little datapoints/misclick')
        
        trial_metrics.nts_distances{trial,:}         = NaN;
        trial_metrics.nts_angles{trial,:}            = NaN;
        trial_metrics.nts_direction_diffs{trial,:}   = NaN;
        trial_metrics.nts_start_times{trial,:}       = NaN;
        trial_metrics.nts_start{trial,:}             = NaN;
        trial_metrics.ntss_diam_before{trial,:}      = NaN;
        trial_metrics.ntss_diam_after{trial,:}       = NaN;
        trial_metrics.nts_end{trial,:}               = NaN;
        trial_metrics.ntse_diam_before{trial,:}      = NaN;
        trial_metrics.ntse_diam_after{trial,:}       = NaN;
        %start_at_target = start_at_target +1;
        trial_metrics.saccades_idx{trial,:}                = NaN;
        trial_metrics.saccades_after_target{trial,:}       = NaN;
        trial_metrics.ts_start{trial,:}                    = NaN;
        trial_metrics.ts_end{trial,:}                      = NaN;
        trial_metrics.tss_diam_before{trial, :}            = NaN;
        trial_metrics.tss_diam_after{trial, :}             = NaN;
        trial_metrics.tse_diam_before{trial, :}            = NaN;
        trial_metrics.tse_diam_after{trial, :}             = NaN;
        trial_metrics.ts_distance_size{trial,:}            = NaN;
        trial_metrics.ts_distance_to_target_start{trial,:} = NaN;
        trial_metrics.ts_distance_to_target_end{trial,:}   = NaN;
        trial_metrics.ts_angle{trial,:}                    = NaN;
        trial_metrics.ts_optimal_angle{trial,:}            = NaN;
        trial_metrics.nts_distances{trial,:}         = NaN;
        trial_metrics.nts_angles{trial,:}            = NaN;
        trial_metrics.nts_direction_diffs{trial,:}   = NaN;
        trial_metrics.nts_start_times{trial,:}       = NaN;
        trial_metrics.nts_start{trial,:}             = NaN;
        trial_metrics.ntss_diam_before{trial,:}      = NaN;
        trial_metrics.ntss_diam_after{trial,:}       = NaN;
        trial_metrics.nts_end{trial,:}               = NaN;
        trial_metrics.ntse_diam_before{trial,:}      = NaN;
        trial_metrics.ntse_diam_after{trial,:}       = NaN;
        trial_metrics.first_saccade_distance(trial, :)   = NaN;
        trial_metrics.first_saccade_angle(trial, :)      = NaN;
        trial_metrics.direction_difference(trial, :)     = NaN;
        trial_metrics.first_saccade_start_time(trial, :) = NaN;
        trial_metrics.optimal_distance(trial, :)     = NaN;
        trial_metrics.optimal_angle(trial, :)        = NaN;
        trial_metrics.blank_diam_avg{trial, :}       = NaN;
        trial_metrics.trial_start_diam_avg{trial, :} = NaN;
        trial_metrics.direction_good(trial, :)       = NaN;
        trial_metrics.n_sac(trial, :)                = NaN;
        trial_metrics.total_dist_sac(trial, :)       = NaN;
        trial_metrics.mean_dist_sac(trial, :)        = NaN;
        continue
    end
    
    blank_diam_r = trial_data.eyeTrial.right.pupil.diameter(idx_blank:idx_blank+baseline_length-1);
    blank_diam_avg = nanmean([blank_diam_l; blank_diam_r],1);
    num_data = size(trial_data.stimulusTrial.left.pupil.diameter,2);
    if baseline_length <= num_data
        trial_start_diam_l = trial_data.stimulusTrial.left.pupil.diameter(1:baseline_length);
        trial_start_diam_r = trial_data.stimulusTrial.right.pupil.diameter(1:baseline_length);
        trial_start_diam_avg = nanmean([trial_start_diam_l; trial_start_diam_r],1);
    else % in case trial is too short fill with nans
        trial_start_diam_l = NaN(1,15);
        trial_start_diam_l(1:num_data) = trial_data.stimulusTrial.left.pupil.diameter(1:end);
        trial_start_diam_r = NaN(1,15);
        trial_start_diam_r(1:num_data) = trial_data.stimulusTrial.right.pupil.diameter(1:end);
        trial_start_diam_avg = nanmean([trial_start_diam_l; trial_start_diam_r],1);
    end
    
    
    % Coordinates of the first gaze of the trial GAZE(x,y)
    %first_gaze_x = trial_data.stimulusTrial.left.gazePoint.onDisplayArea(1, 1) * screenXpixels;
    %first_gaze_y = screenYpixels - (trial_data.stimulusTrial.left.gazePoint.onDisplayArea(2, 1) * screenYpixels);
    % Coordinates of the target TARGET(x,y)
    target_x = fixations.targetCenters(1, trial);
    target_y = fixations.targetCenters(2, trial);
%     % Optimal path between first gaze and target (Euclidean distance)
%     optimal_distance = sqrt((target_x - first_gaze_x)^2 + (target_y - first_gaze_y)^2);
%     % Direction of the optimal path (angle in radians)
%     optimal_angle = atan2(target_y - first_gaze_y, target_x - first_gaze_x);
    
    if ~isempty(current_saccades)
        
        % Coordinates of the first saccade start and end CENTERS(x,y)
        first_saccade_start_x = current_saccades(1).startCenter(1);% * screenXpixels;
        first_saccade_start_y = current_saccades(1).startCenter(2); %screenYpixels - (current_saccades(1).startCenter(1,2) * screenYpixels);
        first_saccade_end_x   = current_saccades(1).endCenter(1); %current_saccades(1).endCenter(1,1) * screenXpixels;
        first_saccade_end_y   = current_saccades(1).endCenter(2); %screenYpixels - (current_saccades(1).endCenter(1,2) * screenYpixels);
        
        % Optimal path between first saccade start and target (Euclidean distance)
        optimal_distance = sqrt((target_x - first_saccade_start_x)^2 + (target_y - first_saccade_start_y)^2);
        % Direction of the optimal path (angle in radians)
        optimal_angle = atan2(target_x - first_saccade_start_x,...
                               target_y - first_saccade_start_y);
           
        % First saccade path (Euclidean distance)
        first_saccade_distance = current_saccades(1).distance;
        % Direction of the first saccade (angle in radians)
         first_saccade_angle = atan2(first_saccade_end_x - first_saccade_start_x, ...
                                     first_saccade_end_y - first_saccade_start_y);
        % Difference between optimal path and first saccade direction (in degrees)
        direction_difference = rad2deg(abs(optimal_angle - first_saccade_angle));
        
        
        
        
        
        
        
        % time until first saccade starts
        first_saccade_start_time = current_saccades(1).saccStartAfterStimOnset;
        % Check if the saccade direction was optimal (i.e., minimal difference)
        if direction_difference <= 10  % Example threshold of 10 degrees
            direction_good = 1;
        else
            direction_good = 0;
        end
        
        % Number of saccades in trial
        n_sac = size(current_saccades, 2);
        % Total distance of saccades (sum of all saccade distances)
        total_dist_sac = sum([current_saccades.distance]);
        % Mean saccade distance
        mean_dist_sac = mean([current_saccades.distance]);
        
        found_target = false;
        nts = false;
        
        
        %%%
        saccades_idx = [];
        saccades_after = [];
        ts_distances_size = [];
        ts_angles = [];
        ts_optimal_angles = [];
        ts_direction_diff_end = [];
        ts_direction_diff_start = [];
        ts_start_times = [];
        ts_end_times = [];
        tss_diam_before = [];
        tss_diam_after = [];
        ts_end_times = [];
        tse_diam_before = [];
        tse_diam_after = [];
        
        %%%
        nts_distances = [];
        nts_angles = [];
        nts_direction_diffs = [];
        nts_start_times = [];
        nts_start = [];
        ntss_diam_before = [];
        ntss_diam_after = [];
        nts_end = [];
        ntse_diam_before = [];
        ntse_diam_after = [];
        %%%
        
        for s = 1:length(current_saccades)
            trial_metrics.trial(trial, :)       = trial;
            
            % Coordinates of the current saccade start and end
            % trial_data.stimulusTrial.left.gazePoint.onDisplayArea(1,:) * screenXpixels
            saccade_start_x = current_saccades(s).startCenter(1);
            saccade_start_y = current_saccades(s).startCenter(2);
            saccade_end_x   = current_saccades(s).endCenter(1);
            saccade_end_y   = current_saccades(s).endCenter(2);
            % Distance to target and direction difference
            distance_to_target_start = sqrt((saccade_end_x - target_x)^2 + (saccade_end_y - target_y)^2);
            distance_to_target_end   = sqrt((saccade_end_x - target_x)^2 + (saccade_end_y - target_y)^2);
            saccade_angle            = atan2(saccade_end_y - saccade_start_y, saccade_end_x - saccade_start_x); % atan2(Y,X)
            optimal_angle            = atan2(target_y - saccade_start_y, target_x - saccade_start_x);
            %direction_difference     = rad2deg(abs(optimal_angle - saccade_angle));
            
            % Check if saccade goes toward the target
            if distance_to_target_end <= tolerance % only safes first saccade towards target so far
                %gaze_found_target_idx = s;
                % if ~found_target % first saccade reaching target
                saccades_until_target = s; %find([current_saccades.startIdx] <= gaze_found_target_idx, 1, 'last');
                
                % Number of saccades after target is found
                if isempty(saccades_until_target) % gaze started at target
                    saccades_until_target = 0;
                    start_at_target = start_at_target +1;
                else
                    saccades_after_target = n_sac - saccades_until_target;
                end
                
                found_target = 1; % Mark that the target was reached
                %                 trial_metrics.saccades_until_target(trial, :) = saccades_until_target; % first target saccade
                %                 trial_metrics.saccades_after_target(trial, :) = saccades_after_target; % first target saccade
                
                saccades_idx = [saccades_idx; saccades_until_target]; % first target saccade
                saccades_after = [saccades_after; saccades_after_target]; % first target saccade
                
                ts_distances_size = [ts_distances_size; current_saccades(s).distance];
                ts_angles = [ts_angles; saccade_angle];
                ts_direction_diff_start = [ts_direction_diff_start; distance_to_target_start];
                ts_direction_diff_start = [ts_direction_diff_start; distance_to_target_end];
                ts_optimal_angles = [ts_optimal_angles, optimal_angle];
                
                if  size(trial_data.stimulusTrial.right.pupil.diameter,2) >= (current_saccades(s).startIdx)+min_length % if enough datapoints
                    % Pupil data before and after target saccade STARTS  TSS
                    ts_start_times = [ts_start_times; current_saccades(s).startIdx];
                    tss_diam_before = [tss_diam_before; get_pupil_data(trial_data, current_saccades(s).startIdx, num_before, 'before')];
                    tss_diam_after = [tss_diam_after; get_pupil_data(trial_data, current_saccades(s).startIdx, num_after, 'after')];
                    % Pupil data before and after target saccade ENDS TSE
                    ts_end_times = [ts_end_times; current_saccades(s).endIdx];
                    tse_diam_before = [tse_diam_before; get_pupil_data(trial_data, current_saccades(s).endIdx, num_before, 'before')];
                    tse_diam_after = [tse_diam_after; get_pupil_data(trial_data, current_saccades(s).endIdx, num_after, 'after')];
                else
                    ts_start_times = [ts_start_times; current_saccades(s).startIdx];
                    trial_metrics.tss_diam_before{trial, :}            = NaN;
                    trial_metrics.tss_diam_after{trial, :}             = NaN;
                    ts_end_times = [ts_end_times; current_saccades(s).endIdx];
                    trial_metrics.tse_diam_before{trial, :}            = NaN;
                    trial_metrics.tse_diam_after{trial, :}             = NaN;
                    
                end
            else
                nts = true;
                nts_distances = [nts_distances; current_saccades(s).distance];
                nts_angles = [nts_angles; saccade_angle];
                nts_direction_diffs = [nts_direction_diffs; direction_difference];
                nts_start_times = [nts_start_times; current_saccades(s).saccStartAfterStimOnset];
                
                % Pupil data before and after non-target saccade start
                nts_start = [nts_start; current_saccades(s).startIdx];
                ntss_diam_before = [ntss_diam_before; get_pupil_data(trial_data, current_saccades(s).startIdx, num_before, 'before')];
                ntss_diam_after = [ntss_diam_after; get_pupil_data(trial_data, current_saccades(s).endIdx, num_after, 'after')];
                % Pupil data before and after non-target saccade end
                nts_end = [nts_end; current_saccades(s).endIdx];
                ntse_diam_before = [ntse_diam_before; get_pupil_data(trial_data, current_saccades(s).endIdx, num_before, 'before')];
                ntse_diam_after = [ntse_diam_after; get_pupil_data(trial_data, current_saccades(s).endIdx, num_after, 'after')];
            end
        end
        % first saccade of trial
        trial_metrics.first_saccade_distance(trial, :)   = first_saccade_distance;
        trial_metrics.first_saccade_angle(trial, :)      = first_saccade_angle;
        trial_metrics.direction_difference(trial, :)     = direction_difference;
        trial_metrics.first_saccade_start_time(trial, :) = first_saccade_start_time;
        
        % general trial metrics
        trial_metrics.optimal_distance(trial, :)     = optimal_distance;
        trial_metrics.optimal_angle(trial, :)        = optimal_angle;
        trial_metrics.blank_diam_avg{trial, :}       = blank_diam_avg;
        trial_metrics.trial_start_diam_avg{trial, :} = trial_start_diam_avg;
        trial_metrics.direction_good(trial, :)       = direction_good;
        trial_metrics.n_sac(trial, :)                = n_sac;
        trial_metrics.total_dist_sac(trial, :)       = total_dist_sac;
        trial_metrics.mean_dist_sac(trial, :)        = mean_dist_sac;
        
        if nts % if any non target saccade didnt go towards target
            trial_metrics.nts_distances{trial,:}         = nts_distances;
            trial_metrics.nts_angles{trial,:}            = nts_angles;
            trial_metrics.nts_direction_diffs{trial,:}   = nts_direction_diffs;
            trial_metrics.nts_start_times{trial,:}       = nts_start_times;
            trial_metrics.nts_start{trial,:}             = nts_start;
            trial_metrics.ntss_diam_before{trial,:}      = ntss_diam_before;
            trial_metrics.ntss_diam_after{trial,:}       = ntss_diam_after;
            trial_metrics.nts_end{trial,:}               = nts_end;
            trial_metrics.ntse_diam_before{trial,:}      = ntse_diam_before;
            trial_metrics.ntse_diam_after{trial,:}       = ntse_diam_after;
        else
            trial_metrics.nts_distances{trial,:}         = NaN;
            trial_metrics.nts_angles{trial,:}            = NaN;
            trial_metrics.nts_direction_diffs{trial,:}   = NaN;
            trial_metrics.nts_start_times{trial,:}       = NaN;
            trial_metrics.nts_start{trial,:}             = NaN;
            trial_metrics.ntss_diam_before{trial,:}      = NaN;
            trial_metrics.ntss_diam_after{trial,:}       = NaN;
            trial_metrics.nts_end{trial,:}               = NaN;
            trial_metrics.ntse_diam_before{trial,:}      = NaN;
            trial_metrics.ntse_diam_after{trial,:}       = NaN;
        end
        
        % If no saccade reached the target, mark with NaNs
        if ~found_target
            trial_metrics.saccades_idx{trial,:}                = NaN;
            trial_metrics.saccades_after_target{trial,:}       = NaN;
            trial_metrics.ts_start{trial,:}                    = NaN;
            trial_metrics.ts_end{trial,:}                      = NaN;
            trial_metrics.tss_diam_before{trial, :}            = NaN;
            trial_metrics.tss_diam_after{trial, :}             = NaN;
            trial_metrics.tse_diam_before{trial, :}            = NaN;
            trial_metrics.tse_diam_after{trial, :}             = NaN;
            trial_metrics.ts_distance_size{trial,:}            = NaN;
            trial_metrics.ts_distance_to_target_start{trial,:} = NaN;
            trial_metrics.ts_distance_to_target_end{trial,:}   = NaN;
            trial_metrics.ts_angle{trial,:}                    = NaN;
            trial_metrics.ts_optimal_angle{trial,:}            = NaN;
        else
            trial_metrics.saccades_idx{trial,:}                = saccades_idx;
            trial_metrics.saccades_after_target{trial,:}       = saccades_after;
            trial_metrics.ts_start{trial,:}                    = ts_start_times;
            trial_metrics.ts_end{trial,:}                      = ts_end_times;
            trial_metrics.tss_diam_before{trial,:}             = tss_diam_before;
            trial_metrics.tss_diam_after{trial, :}             = tss_diam_after;
            trial_metrics.tse_diam_before{trial, :}            = tse_diam_before;
            trial_metrics.tse_diam_after{trial, :}             = tse_diam_after;
            trial_metrics.ts_distance_size{trial,:}            = ts_distances_size;
            trial_metrics.ts_distance_to_target_start{trial,:} = ts_direction_diff_start;
            trial_metrics.ts_distance_to_target_end{trial,:}   = ts_direction_diff_end;
            trial_metrics.ts_angle{trial,:}                    = ts_angles;
            trial_metrics.ts_optimal_angle{trial,:}            = ts_optimal_angles;
            
        end
        
    else % no saccades in trial at all
        start_at_target = start_at_target +1;
        trial_metrics.saccades_idx{trial,:}                = NaN;
        trial_metrics.saccades_after_target{trial,:}       = NaN;
        trial_metrics.ts_start{trial,:}                    = NaN;
        trial_metrics.ts_end{trial,:}                      = NaN;
        trial_metrics.tss_diam_before{trial, :}            = NaN;
        trial_metrics.tss_diam_after{trial, :}             = NaN;
        trial_metrics.tse_diam_before{trial, :}            = NaN;
        trial_metrics.tse_diam_after{trial, :}             = NaN;
        trial_metrics.ts_distance_size{trial,:}            = NaN;
        trial_metrics.ts_distance_to_target_start{trial,:} = NaN;
        trial_metrics.ts_distance_to_target_end{trial,:}   = NaN;
        trial_metrics.ts_angle{trial,:}                    = NaN;
        trial_metrics.ts_optimal_angle{trial,:}            = NaN;
        trial_metrics.nts_distances{trial,:}         = NaN;
        trial_metrics.nts_angles{trial,:}            = NaN;
        trial_metrics.nts_direction_diffs{trial,:}   = NaN;
        trial_metrics.nts_start_times{trial,:}       = NaN;
        trial_metrics.nts_start{trial,:}             = NaN;
        trial_metrics.ntss_diam_before{trial,:}      = NaN;
        trial_metrics.ntss_diam_after{trial,:}       = NaN;
        trial_metrics.nts_end{trial,:}               = NaN;
        trial_metrics.ntse_diam_before{trial,:}      = NaN;
        trial_metrics.ntse_diam_after{trial,:}       = NaN;
        % maybe move this????
        trial_metrics.first_saccade_distance(trial, :)   = NaN;
        trial_metrics.first_saccade_angle(trial, :)      = NaN;
        trial_metrics.direction_difference(trial, :)     = NaN;
        trial_metrics.first_saccade_start_time(trial, :) = NaN;
        
        trial_metrics.optimal_distance(trial, :)     = NaN;
        trial_metrics.optimal_angle(trial, :)        = NaN;
        trial_metrics.blank_diam_avg{trial, :}       = NaN;
        trial_metrics.trial_start_diam_avg{trial, :} = NaN;
        trial_metrics.direction_good(trial, :)       = NaN;
        trial_metrics.n_sac(trial, :)                = 0;
        trial_metrics.total_dist_sac(trial, :)       = NaN;
        trial_metrics.mean_dist_sac(trial, :)        = NaN;
%         trial_metrics.first_saccade_distance(trial, :)   = first_saccade_distance;
%         trial_metrics.first_saccade_angle(trial, :)      = first_saccade_angle;
%         trial_metrics.direction_difference(trial, :)     = direction_difference;
%         trial_metrics.first_saccade_start_time(trial, :) = first_saccade_start_time;
%         
%         trial_metrics.optimal_distance(trial, :)     = optimal_distance;
%         trial_metrics.optimal_angle(trial, :)        = optimal_angle;
%         trial_metrics.blank_diam_avg(trial, :)       = blank_diam_avg;
%         trial_metrics.trial_start_diam_avg(trial, :) = trial_start_diam_avg;
%         trial_metrics.direction_good(trial, :)       = direction_good;
%         trial_metrics.n_sac(trial, :)                = n_sac;
%         trial_metrics.total_dist_sac(trial, :)       = total_dist_sac;
%         trial_metrics.mean_dist_sac(trial, :)        = mean_dist_sac;
    end
end
trial_metrics.start_at_target = start_at_target;




end




%% Helper function to extract pupil data around a saccade
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