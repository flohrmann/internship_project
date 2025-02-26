function [result_table, tnf, tts] = findTargetAndExtractData(cutData, screenXpixels, screenYpixels, analysis_folder, tolerance, num_before, num_after, min_length, baseline_length)
% get baseline starting at blank screen or stimulation screen for each trial

result_table = {};
tnf = 0;  % Counter for trials where the target not found
tts = 0;  % Counter for trials too short
condition_data = containers.Map;  % To store pupil data per condition

% Loop through each trial in the cutData
for trial = 1:size(cutData, 1)
    trial_data = cutData(trial, :);
    condition = trial_data.Condition{1};  % Extract the condition for the trial

    if size(trial_data.stimulusTrial.left.pupil.available, 2) >= min_length
        % --- baselines ---
        % get the pupil diameter during blank screen for baseline
        blank_start = cutData.blankStartTime(trial);
        timestamps = double(trial_data.eyeTrial.systemTimeStamp) / 1e6;
        [~, idx_blank] = min(abs(timestamps - blank_start)); % get start of blank screen
        blank_d_l = trial_data.eyeTrial.left.pupil.diameter;
        blank_d_r = trial_data.eyeTrial.right.pupil.diameter;
        blank_d_avg = nanmean([blank_d_l; blank_d_r], 1); % avg over both eyes
        blank_diam_mean = nanmean(blank_d_avg(idx_blank:idx_blank+baseline_length)); 
        blank_diam_median = nanmedian(blank_d_avg(idx_blank:idx_blank+baseline_length)); 
        
        % get the pupil diameter during first datapoints of trial for baseline
        start_d_l = trial_data.stimulusTrial.left.pupil.diameter;
        start_d_r = trial_data.stimulusTrial.right.pupil.diameter;
        start_d_avg = nanmean([start_d_l; start_d_r], 1); % avg over both eyes
        trial_start_diam_mean = nanmean(start_d_avg(1:baseline_length));
        trial_start_diam_median = nanmedian(start_d_avg(1:baseline_length));

        % Right eye data
        d_r = trial_data.stimulusTrial.right.pupil.diameter;
        % Left eye data
        d_l = trial_data.stimulusTrial.left.pupil.diameter;

        % Calculate the average pupil diameter across both eyes
        d_avg = nanmean([d_r; d_l], 1);

        % Save the data by condition for plotting
        if isKey(condition_data, condition)
            condition_data(condition) = [condition_data(condition); d_avg];
        else
            condition_data(condition) = d_avg;
        end

        % Save the trial data in result_table with all previous values
        result_table = [result_table; {trial, condition, blank_diam_mean, blank_diam_median, ...
            trial_start_diam_mean, trial_start_diam_median, d_avg}];
    else
        tts = tts + 1;
        continue
    end
end

% Convert the cell array to a table for easier processing
result_table = cell2table(result_table, 'VariableNames', {'TrialNumber', 'Condition', 'MeanDiamBlank', ...
    'MedianDiamBlank', 'MeanDiamStart', 'MedianDiamStart', 'PupilDiameterData'});

% Save the data
save(fullfile(analysis_folder, strcat('\pupil_', num2str(num_before), 'before_', num2str(num_after), 'after_finding_stim.mat')), 'result_table');
save(fullfile(analysis_folder, strcat('\tts_pupil_', num2str(num_before), 'before_', num2str(num_after), 'after_finding_stim.mat')), 'tts');
save(fullfile(analysis_folder, strcat('\tnf_pupil_', num2str(num_before), 'before_', num2str(num_after), 'after_finding_stim.mat')), 'tnf');

fprintf('Trials where the target was not found with the eyes: %d \n', tnf);
fprintf('Trials shorter than %d datapoints: %d \n', min_length, tts);
fprintf('Trials left %d \n', height(result_table));

% Plot all trials per condition and average response
plotAllTrialsAndAverage(condition_data);

end

function plotAllTrialsAndAverage(condition_data)
    % Extract unique conditions
    conditions = keys(condition_data);
    num_conditions = length(conditions);

    figure;
    % Subplots for each condition showing individual trials
    for i = 1:num_conditions
        condition = conditions{i};
        data = condition_data(condition);

        subplot(num_conditions, 1, i);
        hold on;
        for j = 1:size(data, 1)
            plot(data(j, :), 'Color', [0.8, 0.8, 0.8]); % Plot individual trials in light grey
        end
        avg_response = nanmean(data, 1); % Calculate average response for condition
        plot(avg_response, 'k', 'LineWidth', 1.5); % Plot average in black
        hold off;
        title(sprintf('Condition: %s', condition));
        ylabel('Pupil Diameter');
    end
    xlabel('Time (s)');

    % Separate plot for average response per condition
    figure;
    hold on;
    for i = 1:num_conditions
        condition = conditions{i};
        data = condition_data(condition);
        avg_response = nanmean(data, 1);
        plot(avg_response, 'DisplayName', condition, 'LineWidth', 1.5);
    end
    hold off;
    title('Average Pupil Response per Condition');
    xlabel('Time (s)');
    ylabel('Pupil Diameter');
    legend('show');
end
