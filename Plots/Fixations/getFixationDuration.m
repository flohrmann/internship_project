function fixationDurations  = getFixationDuration(all_fixations, data, conditions, comparison_results_folder)
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
        
        % avg fixation duration per trial
        durations = arrayfun(@(x) x.fixEnd-x.fixStart, fixations);
        avg_fixation = mean(durations);
        median_fixation = median(durations);
        
        % Store the fixation durations, number of fixations, and search start time
        fixationDurations(participant).trials(trial).conditions = condition;
        fixationDurations(participant).trials(trial).durations = trialDurations;
        fixationDurations(participant).trials(trial).meanFixDuration = avg_fixation;
        fixationDurations(participant).trials(trial).medianFixDuration = median_fixation;       
        fixationDurations(participant).trials(trial).numFixations = numFixations;  % Store the number of fixations
        fixationDurations(participant).trials(trial).searchStartTime = searchStartTime;  % Store the time to start searching
        
    end
    
    
    conditions_list = []; % to get a mask per condition
    for trial = 1:length(fixationDurations(participant).trials)
        conditions_list = [conditions_list; string(data(participant).Condition(trial))];
    end
    
    condition_mask = strcmp(conditions{1}, conditions_list);
    mean_fix_durations = arrayfun(@(x) x.meanFixDuration, fixationDurations(participant).trials(condition_mask));
    fixationDurations(participant).a_mean = nanmean(mean_fix_durations);
    median_fix_durations = arrayfun(@(x) x.medianFixDuration, fixationDurations(participant).trials(condition_mask));
    fixationDurations(participant).a_median = nanmedian(median_fix_durations);
    fixationDurations(participant).a_sem = std(mean_fix_durations, 'omitnan');
    
    condition_mask = strcmp(conditions{2}, conditions_list);
    mean_fix_durations = arrayfun(@(x) x.meanFixDuration, fixationDurations(participant).trials(condition_mask));
    fixationDurations(participant).as_mean = nanmean(mean_fix_durations);
    median_fix_durations = arrayfun(@(x) x.medianFixDuration, fixationDurations(participant).trials(condition_mask));
    fixationDurations(participant).as_median = nanmedian(median_fix_durations);
    fixationDurations(participant).as_sem = std(mean_fix_durations, 'omitnan');
    
    condition_mask = strcmp(conditions{3}, conditions_list);
    mean_fix_durations = arrayfun(@(x) x.meanFixDuration, fixationDurations(participant).trials(condition_mask));
    fixationDurations(participant).b_mean = nanmean(mean_fix_durations);
    median_fix_durations = arrayfun(@(x) x.medianFixDuration, fixationDurations(participant).trials(condition_mask));
    fixationDurations(participant).b_median = nanmedian(median_fix_durations);
    fixationDurations(participant).b_sem = std(mean_fix_durations, 'omitnan');
    
    condition_mask = strcmp(conditions{4}, conditions_list);
    mean_fix_durations = arrayfun(@(x) x.meanFixDuration, fixationDurations(participant).trials(condition_mask));
    fixationDurations(participant).bs_mean = nanmean(mean_fix_durations);
    median_fix_durations = arrayfun(@(x) x.medianFixDuration, fixationDurations(participant).trials(condition_mask));
    fixationDurations(participant).bs_median = nanmedian(median_fix_durations);
    fixationDurations(participant).bs_sem = std(mean_fix_durations, 'omitnan');
    
end

% Save the fixation durations and counts data
save(fullfile(comparison_results_folder, 'fixationDurations.mat'), 'fixationDurations');
end
