function [diam_t0, tnf, tts] = findTargetAndExtractData(cutData, screenXpixels, screenYpixels, analysis_folder, tolerance, num_before, num_after, min_length, baseline_length)
% get baseline starting at blank screen or stimulation screen for each trial

result_table = {};
tnf = 0;  % Counter for trials where the target not found
tts = 0; % Counter where trials too short
condition_data = struct();

for trial = 1:size(cutData, 1)
    trial_data = cutData(trial, :);
    condition = trial_data.Condition{1};  % Extract the condition for the trial
    
    if size(trial_data.stimulusTrial.left.pupil.available, 2)>=min_length
        
        % --- baselines ---
        % get the pupil diameter during blank screen for baseline
        blank_start = cutData.blankStartTime(trial);
        timestamps = double(trial_data.eyeTrial.systemTimeStamp) / 1e6;
        [~, idx_blank] = min(abs(timestamps - blank_start)); % get start of blank screen
        blank_d_l = trial_data.eyeTrial.left.pupil.diameter;
        blank_d_r = trial_data.eyeTrial.right.pupil.diameter;
        %blank_d_avg = nanmean([blank_d_l; blank_d_r], 1); % avg over both eyes
        
        % cut only data used for blank baseline
        blank_diam_l = trial_data.eyeTrial.left.pupil.diameter(idx_blank:idx_blank+baseline_length-1);
        blank_diam_r = trial_data.eyeTrial.right.pupil.diameter(idx_blank:idx_blank+baseline_length-1);
        blank_diam_avg = nanmean([blank_diam_l; blank_diam_r],1);
        blank_diam_mean =  nanmean(blank_diam_avg); 
        blank_diam_median = nanmedian(blank_diam_avg); 
        
        % get the pupil diameter during first datapoints of trial for baseline
        start_d_l = trial_data.stimulusTrial.left.pupil.diameter;
        start_d_r = trial_data.stimulusTrial.right.pupil.diameter;
        %start_d_avg = nanmean([start_d_l; start_d_r], 1); % avg over both eyes
        
        trial_start_diam_l = trial_data.stimulusTrial.left.pupil.diameter(1:baseline_length);
        trial_start_diam_r = trial_data.stimulusTrial.right.pupil.diameter(1:baseline_length);
        trial_start_diam_avg = nanmean([trial_start_diam_l; trial_start_diam_r],1); 
        %trial_start_diam_mean = nanmean(trial_start_diam_avg); % avg over time of trial start screen 
        %trial_start_diam_median = nanmedian(trial_start_diam_avg); % avg over time of trial start screen
        
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
        % CHANGE made it the last time!
        %found_target_idx = find(sqrt((x_avg - targetX).^2 + (y_avg - targetY).^2) <= tolerance, 1, 'last');
        
        % Check if the target was found
        if isempty(found_target_idx)
            %fprintf('Target not found in trial %d. Skipping this trial.\n', trial);
            tnf = tnf + 1;
            continue;  % Skip to the next trial
        end
        
        %% Extract the last X points before finding the target
        %found_target_idx = found_target_idx -1; % found target counts as 'after finding target already'
        if found_target_idx > num_before
            points_before_stim = d_avg(found_target_idx - num_before: found_target_idx-1);
        else
            % If there are fewer than `num_before` points, pad with NaNs
            points_before_stim = [NaN(1, num_before - found_target_idx), d_avg(1:found_target_idx)];
        end
        
        %% Extract the X points after finding the target
        if found_target_idx + num_after < length(x_avg)
            points_after_stim = d_avg(found_target_idx : found_target_idx + num_after-1);
        else
            % If there are fewer than `num_after` points, pad with NaNs
            points_after_stim = [d_avg(found_target_idx +1: end), NaN(1, num_after - (length(x_avg) - found_target_idx))];
        end
        
        % Ensure that both before and after arrays are the correct length
        if length(points_before_stim) > num_before
            points_before_stim = points_before_stim(1:num_before);
        end
        if length(points_after_stim) > num_after
            points_after_stim = points_after_stim(1:num_after);
        end
        % Get X points around finding the target
        %points_around_stim = [points_before_stim, points_after_stim];
%         
%         %% plot the raw cut data
%             %y_limit_avg = []; y_limits = [];
%             time_vector = (-num_before:num_after - 1) * sr; % Time in seconds
%             num_rows = 1;
%             num_col = 2;
%             figure; sgtitle('Pupil Response Before and After Gaze Reached Target (Raw)')
%             subplot(num_rows, num_col, 1); hold on; %y_limits = [];
%             for i = 1:4
%                 condition = conditions{i}; 
%                 condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
%                 for j = 1:height(condition_trials)
%                     current_diam = [condition_trials.DataPointsBefore(j, :), condition_trials.DataPointsAfter(j, :)];
%                     plot(time_vector, current_diam, 'Color', color_map(condition)); % Plot each trial in light grey        y_limits = [y_limits; ylim];  % Store y-axis limits
%                 end % y_limits = [y_limits; ylim];  % Store y-axis limits
%             end
%             title('Individual Trials'); xlabel('Time (s)'); ylabel('Pupil Diameter (mm)');
%             hold off;
% 
%             subplot(num_rows, num_col, 2); hold on; 
%             for i = 1:4
%                 condition = conditions{i};
%                 condition_trials = diam_t0(strcmp(diam_t0.Condition, condition), :);
%                 current_diams = [condition_trials.DataPointsBefore, condition_trials.DataPointsAfter];
%                 avg_response = nanmean(current_diams);
%                 plot(time_vector, avg_response, 'Color', color_map(condition), 'DisplayName', condition, 'LineWidth', 1.5);
%             end
%             hold off;
%             title('Average per Condition'); xlabel('Time (s)'); ylabel('Pupil Diameter (mm)');
%             legend('show'); legend('boxoff'); % y_limit_avg = [y_limit_avg; ylim];

        % Save the trial number, condition, and datapoints in a row
        result_table = [result_table; {trial, condition, points_before_stim, points_after_stim, ...
                        blank_diam_l, blank_diam_r, blank_diam_avg, ... % blank_diam_mean, blank_diam_median, ...
                        trial_start_diam_l, trial_start_diam_r, trial_start_diam_avg}];%, trial_start_diam_mean, trial_start_diam_median}];
    else
        tts = tts +1;
        continue
    end
end

% Convert the cell array to a table for easier processing
diam_t0 = cell2table(result_table, 'VariableNames', {'TrialNumber', 'Condition', 'DataPointsBefore', 'DataPointsAfter', ...
    'DiamLBlank', 'DiamRBlank', 'DiamBothBlank', ... % 'MeanDiamBlank', 'MedianDiamBlank', ...
    'DiamLStart', 'DiamRStart', 'DiamBothStart'});%, 'MeanDiamStart', 'MedianDiamStart'});



save(fullfile(analysis_folder, strcat('\pupil_', num2str(num_before), 'before_', num2str(num_after), 'after_finding_stim.mat')), 'diam_t0');
save(fullfile(analysis_folder, strcat('\tts_pupil_', num2str(num_before), 'before_', num2str(num_after), 'after_finding_stim.mat')), 'tts');
save(fullfile(analysis_folder, strcat('\tnf_pupil_', num2str(num_before), 'before_', num2str(num_after), 'after_finding_stim.mat')), 'tnf');

fprintf('Trials where the target was not found with the eyes: %d \n', tnf);
fprintf('Trials shorter than %d datapoints: %d \n', min_length, tts);
fprintf('Trials left %d \n',height(diam_t0));

end
