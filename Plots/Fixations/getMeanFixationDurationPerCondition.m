function meanFixationDurations = getMeanFixationDurationPerCondition(fixationStats, conditions)
    % Initialize structure to store the mean fixation durations per participant and condition
    meanFixationDurations = struct();
    
    % Get the number of conditions
    numConditions = length(conditions);
    
    % Loop through each participant
    for participant = 1:length(fixationStats)
        % Get the participant ID and initialize storage for condition means
        meanFixationDurations(participant).id = fixationStats(participant).id;
        meanFixationDurations(participant).conditionMeans = NaN(1, numConditions);  % Preallocate means array

        % Loop through each condition
        for c = 1:numConditions
            % Find the fixations for the current condition
            conditionName = conditions{c};
            conditionIdx = find(strcmp({fixationStats(participant).trials.conditions}, conditionName));
            
            if ~isempty(conditionIdx)
                % Extract fixation durations for this condition
                fixationDurations = fixationStats(participant).trials(conditionIdx).fixationDurations;
                
                if ~isempty(fixationDurations)
                    % Calculate the mean fixation duration for this condition
                    meanFixationDurations(participant).conditionMeans(c) = mean(fixationDurations, 'omitnan');
                else
                    meanFixationDurations(participant).conditionMeans(c) = NaN;
                end
            end
        end
    end
    
    % Display the results (optional)
    disp('Mean Fixation Durations per Participant and Condition:');
    disp(meanFixationDurations);
end
