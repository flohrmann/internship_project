function cutData = cutEyeTrackingData(folder_name, trial_results, eye_tracking_data)
    % Get the system timestamps
    timestamps = double(eye_tracking_data.systemTimeStamp)/1000000; % Convert to milliseconds

    % Number of trials
    numTrials = size(trial_results, 1);

    % Initialize the output struct array
    cutData = struct([]);

    for trial = 1:numTrials
        % Get the start and end times for the current trial
        current_data = trial_results(trial,:);
        fixStartTime = current_data.fixStartTime;
        trialEndTime = current_data.trialEndTime;
        StimulusOnsetTime = current_data.StimulusOnsetTime;
        
        % Find the closest timestamps to fixStartTime, StimulusOnsetTime, and trialEndTime
        [~, startIdx] = min(abs(timestamps - fixStartTime));
        [~, stimulusIdx] = min(abs(timestamps - StimulusOnsetTime));
        [~, endIdx] = min(abs(timestamps - trialEndTime));
        
        %%%
%         figure(1); clf; plot(eye_tracking_data.systemTimeStamp/1e3, '.');
%         hold on;
%         plot(trial_results.fixStartTime*1000, '*');
%         hold on;
%         plot(trial_results.StimulusOnsetTime*1000, '.');
% 
%         figure(2); clf;
%         plot(trial_results.fixStartTime*1000, 'r*'); hold on;
%         plot(trial_results.StimulusOnsetTime*1000, 'b.');

        
        %% plot test
%         FixationStart_times = trial_results.fixStartTime; % seconds
%         ButtonPress_times = trial_results.trialEndTime;
%         StimulusStart_time  = trial_results.StimulusOnsetTime; % seconds
%         
%         eye_tracker_times = eye_tracking_data.systemTimeStamp; % microseconds
%         eye_tracker_times_inSeconds  = double(eye_tracker_times)/1000000;
%         NTrials = length(FixationStart_times);
%         Fixation_Index_InEyeTrackerTimeStamps = zeros(1, NTrials);
%         StimulusOnset_Index_InEyeTrackerTimeStamps = zeros(1, NTrials);
%         ButtonPressend_Index_MatlabTimeStamp = zeros(1, NTrials);
% 
%         for t = 1:NTrials
%                 [~, Fixation_Index_InEyeTrackerTimeStamps(t)] = min(abs(eye_tracker_times_inSeconds -FixationStart_times(t)));  
%                 [~, StimulusOnset_Index_InEyeTrackerTimeStamps(t)] = min(abs(eye_tracker_times_inSeconds -StimulusStart_time(t)));
%                 [~, ButtonPressend_Index_MatlabTimeStamp(t)] = min(abs(eye_tracker_times_inSeconds - ButtonPress_times(t)));
%         end
%         figure(10); clf; plot(eye_tracker_times_inSeconds, 'r.'); hold;
%         plot(Fixation_Index_InEyeTrackerTimeStamps, eye_tracker_times_inSeconds(Fixation_Index_InEyeTrackerTimeStamps), 'b*');
%         plot(StimulusOnset_Index_InEyeTrackerTimeStamps, eye_tracker_times_inSeconds(StimulusOnset_Index_InEyeTrackerTimeStamps), 'g>');
%         plot(ButtonPressend_Index_MatlabTimeStamp, eye_tracker_times_inSeconds(ButtonPressend_Index_MatlabTimeStamp), 'ks');
%         
%         
%         
%         
        cutData(trial).eyeTrial = extractEyeTrackingData(eye_tracking_data, startIdx, endIdx);
        cutData(trial).stimulusTrial = extractEyeTrackingData(eye_tracking_data, stimulusIdx, endIdx);
        
        cutData(trial).startIdx = startIdx;
        cutData(trial).stimulusIdx = stimulusIdx;
        cutData(trial).endIdx = endIdx;
%         
        % Only save for trials with eye tracking data
        if endIdx == startIdx
            cutData(trial).withTracking = 0;
        else
            cutData(trial).withTracking = 1;
        end
        % Add trial-specific data
        cutData(trial).AngleMatrix = current_data.AngleMatrix;
        cutData(trial).TargetPosition = current_data.TargetPosition;
        cutData(trial).targetSide = current_data.TargetSide;
        cutData(trial).Condition = current_data.Condition;
        cutData(trial).conditionType = current_data.ConditionType;
        cutData(trial).cell_width = current_data.cell_width;
        cutData(trial).cell_height = current_data.cell_height;
        cutData(trial).line_length = current_data.line_length;
        cutData(trial).line_width = current_data.line_width;
        cutData(trial).x_jitters = current_data.x_jitters{1};
        cutData(trial).y_jitters = current_data.y_jitters{1};
        cutData(trial).x_centers = current_data.x_centers{1};
        cutData(trial).y_centers = current_data.y_centers{1};
        cutData(trial).response = current_data.response;
        cutData(trial).trialStartTime = current_data.trialStartTime;
        cutData(trial).blankStartTime = current_data.blankStartTime;
        cutData(trial).fixStartTime = fixStartTime;
        cutData(trial).StimulusOnsetTime = StimulusOnsetTime;
        cutData(trial).trialEndTime = trialEndTime;
        cutData(trial).rt = current_data.rt;
        cutData(trial).fixTime = current_data.fixTime;     
    end

    % Convert struct array to table
    cutData = struct2table(cutData);
    
    % Save the table to a .mat file
    results_file = [folder_name, '\cut_trials_' , datestr(now, 'yyyymmdd_HHMM'), '.mat'];
    save(results_file, 'cutData');
    disp('Eye tracking data cut into smaller structs per trial and saved as a table.');
end

function eyeTrackingData = extractEyeTrackingData(eye_tracking_data, startIdx, endIdx)
    eyeTrackingData = struct( ...
        'systemTimeStamp', eye_tracking_data.systemTimeStamp(startIdx:endIdx), ...
        'deviceTimeStamp', eye_tracking_data.deviceTimeStamp(startIdx:endIdx), ...
        'left', struct( ...
            'gazePoint', struct( ...
                'onDisplayArea', eye_tracking_data.left.gazePoint.onDisplayArea(:, startIdx:endIdx), ...
                'inUserCoords', eye_tracking_data.left.gazePoint.inUserCoords(:, startIdx:endIdx), ...
                'valid', eye_tracking_data.left.gazePoint.valid(:, startIdx:endIdx), ...
                'available', eye_tracking_data.left.gazePoint.available(:, startIdx:endIdx)), ...
            'pupil', struct( ...
                'diameter', eye_tracking_data.left.pupil.diameter(:, startIdx:endIdx), ...
                'valid', eye_tracking_data.left.pupil.valid(:, startIdx:endIdx), ...
                'available', eye_tracking_data.left.pupil.available(:, startIdx:endIdx)), ...
            'gazeOrigin', struct( ...
                'inUserCoords', eye_tracking_data.left.gazeOrigin.inUserCoords(:, startIdx:endIdx), ...
                'inTrackBoxCoords', eye_tracking_data.left.gazeOrigin.inTrackBoxCoords(:, startIdx:endIdx), ...
                'valid', eye_tracking_data.left.gazeOrigin.valid(:, startIdx:endIdx), ...
                'available', eye_tracking_data.left.gazeOrigin.available(:, startIdx:endIdx)), ...
            'eyeOpenness', struct( ...
                'diameter', eye_tracking_data.left.eyeOpenness.diameter(:, startIdx:endIdx), ...
                'valid', eye_tracking_data.left.eyeOpenness.valid(:, startIdx:endIdx), ...
                'available', eye_tracking_data.left.eyeOpenness.available(:, startIdx:endIdx))), ...
        'right', struct( ...
            'gazePoint', struct( ...
                'onDisplayArea', eye_tracking_data.right.gazePoint.onDisplayArea(:, startIdx:endIdx), ...
                'inUserCoords', eye_tracking_data.right.gazePoint.inUserCoords(:, startIdx:endIdx), ...
                'valid', eye_tracking_data.right.gazePoint.valid(:, startIdx:endIdx), ...
                'available', eye_tracking_data.right.gazePoint.available(:, startIdx:endIdx)), ...
            'pupil', struct( ...
                'diameter', eye_tracking_data.right.pupil.diameter(:, startIdx:endIdx), ...
                'valid', eye_tracking_data.right.pupil.valid(:, startIdx:endIdx), ...
                'available', eye_tracking_data.right.pupil.available(:, startIdx:endIdx)), ...
            'gazeOrigin', struct( ...
                'inUserCoords', eye_tracking_data.right.gazeOrigin.inUserCoords(:, startIdx:endIdx), ...
                'inTrackBoxCoords', eye_tracking_data.right.gazeOrigin.inTrackBoxCoords(:, startIdx:endIdx), ...
                'valid', eye_tracking_data.right.gazeOrigin.valid(:, startIdx:endIdx), ...
                'available', eye_tracking_data.right.gazeOrigin.available(:, startIdx:endIdx)), ...
            'eyeOpenness', struct( ...
                'diameter', eye_tracking_data.right.eyeOpenness.diameter(:, startIdx:endIdx), ...
                'valid', eye_tracking_data.right.eyeOpenness.valid(:, startIdx:endIdx), ...
                'available', eye_tracking_data.right.eyeOpenness.available(:, startIdx:endIdx))));
end
