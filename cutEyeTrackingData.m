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
        % trialStartTime: 'press any button to continue' screen
        trialStartTime = current_data.trialStartTime;
        blankStartTime = current_data.blankStartTime;
        % fixSTartTime: 'fixation cross' screen
        fixStartTime = current_data.fixStartTime;
        % StimulusOnsetTime: 'stimulation' screen (0.7s after fixSTartTime)
        StimulusOnsetTime = current_data.StimulusOnsetTime;
        % trialEndTime: time L/R button was pressed to end trial
        trialEndTime = current_data.trialEndTime;
        
        % Find the closest timestamps to fixStartTime, StimulusOnsetTime, and trialEndTime
        [~, trialStartIdx] = min(abs(timestamps - trialStartTime));
        [~, startIdx] = min(abs(timestamps - fixStartTime));
        [~, blankIdx] = min(abs(timestamps - blankStartTime));
        [~, stimulusIdx] = min(abs(timestamps - StimulusOnsetTime));
        [~, endIdx] = min(abs(timestamps - trialEndTime));
        
        %% plot test
        FixationStart_times = trial_results.fixStartTime; % seconds
        ButtonPress_times = trial_results.trialEndTime;
        StimulusStart_time  = trial_results.StimulusOnsetTime; % seconds
        
        eye_tracker_times = eye_tracking_data.systemTimeStamp; % microseconds
        eye_tracker_times_inSeconds  = double(eye_tracker_times)/1000000;
        NTrials = length(FixationStart_times);
        Fixation_Index_InEyeTrackerTimeStamps = zeros(1, NTrials);
        StimulusOnset_Index_InEyeTrackerTimeStamps = zeros(1, NTrials);
        ButtonPressend_Index_MatlabTimeStamp = zeros(1, NTrials);

        for t = 1:NTrials
                [~, Fixation_Index_InEyeTrackerTimeStamps(t)] = min(abs(eye_tracker_times_inSeconds -FixationStart_times(t)));  
                [~, StimulusOnset_Index_InEyeTrackerTimeStamps(t)] = min(abs(eye_tracker_times_inSeconds -StimulusStart_time(t)));
                [~, ButtonPressend_Index_MatlabTimeStamp(t)] = min(abs(eye_tracker_times_inSeconds - ButtonPress_times(t)));
        end
       
        %% cut the data         
        % complete trial
        cutData(trial).startIdxWithInstructionScreen = trialStartIdx;
        cutData(trial).completeTrial = extractEyeTrackingData(eye_tracking_data, trialStartIdx, endIdx);
        % trial after user presses start trial button: fixation, blank screen and stimulation screen
        cutData(trial).startIdx = startIdx;
        cutData(trial).eyeTrial = extractEyeTrackingData(eye_tracking_data, startIdx, endIdx);
        
        if size(cutData(trial).eyeTrial.systemTimeStamp, 2) < 25
            disp('less than 25 datapoints in this trial');
            %disp(size(cutData(trial).eyeTrial.systemTimeStamp, 2));
        else
        end
        cutData(trial).blankIdx = blankIdx; % blank screen index
        % only data during stimulation presentatiopn
        cutData(trial).stimulusIdx = stimulusIdx;
        cutData(trial).stimulusTrial = extractEyeTrackingData(eye_tracking_data, stimulusIdx, endIdx);
        cutData(trial).endIdx = endIdx;
  
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
    %results_file = [folder_name, '\cut_trials_' , datestr(now, 'yyyymmdd_HHMM'), '.mat'];
    results_file = [folder_name, '\cut_trials.mat'];
    save(results_file, 'cutData');
    disp('Eye tracking data cut into smaller structs per trial and saved as a table.');
end



%% restructure the struct into original eyetracker format
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
