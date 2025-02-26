function searchLatencyStats = analyzeSearchLatency(all_fixations, data, conditions, comparison_results_folder)
    searchLatencyStats = struct();
    
    for participant = 1:length(all_fixations)
        participant_fixations = all_fixations(1, participant).fixations;
        searchLatencyStats(participant).id = all_fixations(participant).id;
        searchLatencyStats(participant).conditions = struct('name', [], 'searchLatencies', [], ...
                                                            'searchDurations', []);

        % Loop over each trial
        for trial = 1:length(participant_fixations.stimulusFixations)
            fixations = participant_fixations.stimulusFixations(trial).fixations;
            condition = data(participant).Condition{trial};
            conditionIdx = find(strcmp(conditions, condition));
            
            % Initialize latency and duration metrics
            searchLatency = NaN;
            searchDuration = NaN;
            
            % Assuming that `rt_eye` in `data` indicates the time of eye movements
            stimulusOnset = data(participant).StimulusOnsetTime(trial);
            eyeMovementTimes = data(participant).rt_eye(trial);
            
            if ~isempty(eyeMovementTimes) && length(fixations) > 1
                % Search latency: Time from stimulus onset to the first saccade/fixation shift
                searchLatency = eyeMovementTimes(1) - stimulusOnset;
                
                % Search duration: Time from first to last significant saccade/fixation
                searchDuration = eyeMovementTimes(end) - eyeMovementTimes(1);
            end
            
            % Store the latency and duration per trial
            searchLatencyStats(participant).conditions(conditionIdx).name = conditions{conditionIdx};
            searchLatencyStats(participant).conditions(conditionIdx).searchLatencies = ...
                [searchLatencyStats(participant).conditions(conditionIdx).searchLatencies, searchLatency];
            searchLatencyStats(participant).conditions(conditionIdx).searchDurations = ...
                [searchLatencyStats(participant).conditions(conditionIdx).searchDurations, searchDuration];
        end
    end
    
    % Save the searchLatencyStats struct to a .mat file
    save(fullfile(comparison_results_folder, 'searchLatencyStats.mat'), 'searchLatencyStats');
end
