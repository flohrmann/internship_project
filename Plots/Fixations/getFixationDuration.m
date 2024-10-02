function fixationDurations  = getFixationDuration(all_fixations, data, comparison_results_folder)
    fixationDurations = struct();  % Structure to store fixation durations and counts for each participant
    
    % Loop over each participant
    for participant = 1:length(all_fixations)
        fixationDurations(participant).id = all_fixations(participant).id;  % Store participant ID
        fixationDurations(participant).trials = struct();  % Initialize trials
        
        % Loop over each trial
        for trial = 1:length(all_fixations(participant).fixations.stimulusFixations)
            fixations = all_fixations(participant).fixations.stimulusFixations(trial).fixations;  % Get fixations for this trial
            
            trialDurations = [];  % Initialize storage for fixation durations in this trial
            numFixations = length(fixations);  % Count the number of fixations in this trial
            
            firstFixationEndTime = NaN;  % Initialize search start time (time when first fixation is left)

            % Get stimulus onset time for the current trial
            stimulusOnsetTime = data(participant).StimulusOnsetTime(trial); 
            condition = data(participant).Condition(trial);
            
            % Loop over each fixation
            for fIdx = 1:numFixations
                % Get the start and end indices of the fixation
                startIdx = fixations(fIdx).startIdx;
                endIdx = fixations(fIdx).endIdx;
                
                % Retrieve the corresponding timestamps from the eye-tracking data
                startTime = double(data(participant).stimulusTrial(trial).systemTimeStamp(startIdx)) / 1e6;
                endTime = double(data(participant).stimulusTrial(trial).systemTimeStamp(endIdx)) / 1e6;
                
                % Calculate the fixation duration
                fixationDuration = endTime - startTime;
                
                % Store the duration
                trialDurations = [trialDurations, fixationDuration];
                
                % Store the time when the first fixation is left (i.e., end time of first fixation)
                if fIdx == 1
                    firstFixationEndTime = endTime;
                end
            end
            
            % Calculate the time it took the participant to start searching
            if ~isnan(firstFixationEndTime)
                searchStartTime = firstFixationEndTime - stimulusOnsetTime;  % Time from stimulus onset to end of first fixation
            else
                searchStartTime = NaN;  % If there's no fixation, return NaN
            end
            
            % Store the fixation durations, number of fixations, and search start time
            fixationDurations(participant).trials(trial).conditions = condition;
            fixationDurations(participant).trials(trial).durations = trialDurations;
            fixationDurations(participant).trials(trial).numFixations = numFixations;  % Store the number of fixations
            fixationDurations(participant).trials(trial).searchStartTime = searchStartTime;  % Store the time to start searching
        end
    end
    
    % Save the fixation durations and counts data
    save(fullfile(comparison_results_folder, 'fixationDurations.mat'), 'fixationDurations');
end
