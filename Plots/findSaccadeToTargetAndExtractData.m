function [diam_t0_s, tnf, tts] = findSaccadeToTargetAndExtractData(cutData, trial_metrics, saccades, analysis_folder, num_before, num_after, min_length)
% get baseline starting at blank screen or stimulation screen for each trial

result_table = {};
tnf = 0;  % Counter for trials where the target not found
tts = 0; % Counter where trials too short
for trial = 1:size(cutData, 1)
    trial_data = cutData(trial, :);
    condition = trial_data.Condition{1};
    if size(trial_data.stimulusTrial.left.pupil.available, 2)>=min_length

            points_before_stim = trial_metrics.ts_diam_before(trial,:);
            points_after_stim = trial_metrics.ts_diam_after(trial,:);

            DiamBothStart = trial_metrics.trial_start_diam_avg(trial,:);
            DiamBothBlank = trial_metrics.blank_diam_avg(trial,:);

            saccade_idx = trial_metrics.ts_start(trial,:);

            % index of saccade towards target
            found_target_idx_saccade = trial_metrics.saccades_until_target(trial);
            % Check if the target was found
            if isempty(found_target_idx_saccade)
                %fprintf('Target not found in trial %d. Skipping this trial.\n', trial);
                tnf = tnf + 1;
            else
                % Save the trial number, condition, and datapoints in a row
                result_table = [result_table; {trial, condition, saccade_idx, points_before_stim, points_after_stim, DiamBothStart, DiamBothBlank}];
                continue;  % Skip to the next trial
            end
    else
        tts = tts +1;
        continue
    end
end

% Convert the cell array to a table for easier processing
diam_t0_s = cell2table(result_table, 'VariableNames', {'TrialNumber', 'Condition', 'TS_index','DataPointsBefore', 'DataPointsAfter', 'DiamBothStart', 'DiamBothBlank'});

    save(fullfile(analysis_folder, strcat('\pupil_', num2str(num_before), 'before_', num2str(num_after), 'after_saccade_to_target.mat')), 'diam_t0_s');
    save(fullfile(analysis_folder, strcat('\tts_pupil_', num2str(num_before), 'before_', num2str(num_after), 'after_saccade_to_target.mat')), 'tts');
    save(fullfile(analysis_folder, strcat('\tnf_pupil_', num2str(num_before), 'before_', num2str(num_after), 'after_saccade_to_target.mat')), 'tnf');

    fprintf('Trials where the target was not found with the eyes: %d \n', tnf);
    fprintf('Trials shorter than %d datapoints: %d \n', min_length, tts);
    fprintf('Trials left %d \n',height(diam_t0_s));

end
