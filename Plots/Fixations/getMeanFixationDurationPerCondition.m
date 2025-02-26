function meanFixationDurations = getMeanFixationDurationPerCondition(fixationStats, conditions)
    % Initialize structure to store the mean fixation durations per participant and condition
    meanFixationDurations = struct();
    
    % Get the number of conditions
    numConditions = length(conditions);
    
    % Loop through each participant
    for participant = 1:size(fixationStats,2)
        % Get the participant ID and initialize storage for condition means
        meanFixationDurations(participant).id = fixationStats(participant).id;
        meanFixationDurations(participant).conditionMeans = NaN(1, numConditions);  % Preallocate means array

        conditions_list = []; % to get a mask per condition
        for trial = 1:length(fixationStats(participant).trials)
            conditions_list = [conditions_list; string(fixationStats(participant).trials(trial).conditions)];
        end
            
        % Loop through each condition
        for c = 1:numConditions
            % Find the fixations for the current condition
            condition = conditions{c};
            conditionIdx = strcmp(condition, conditions_list);
            
            if ~isempty(conditionIdx) % Extract fixation durations for this condition
                for trial = 1:length(fixationStats(participant).trials)
                    fixationDurations = fixationStats(participant).trials(conditionIdx(trial, 1)).durations;
                end
                
                if ~isempty(fixationDurations)
                    % Calculate the mean fixation duration for this condition
                    meanFixationDurations(participant).conditionMeans(c) = nanmean(fixationDurations, 'omitnan');
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
