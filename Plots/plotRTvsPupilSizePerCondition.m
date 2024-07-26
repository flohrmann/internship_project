function plotRTvsPupilSizePerCondition(trial_results, eye_tracking_data, analysis_folder, color_map)
    % Extract unique conditions from trial_results
    uniqueConditions = categorical(trial_results.Condition);
    uniqueConditions = categories(uniqueConditions);
    
    % Initialize figure for subplots
    figure;
    
    % Loop through each condition and create subplots
    for c = 1:length(uniqueConditions)
        condition = uniqueConditions{c};
        
        % Initialize arrays to store pupil sizes and reaction times
        pupilSizes = [];
        reactionTimes = [];

        % Loop through each trial in trial_results
        for trial = 1:size(trial_results, 1)
            % Extract the current trial data
            current_data = trial_results(trial, :);

            % Check if the current trial matches the condition
            if strcmp(current_data.Condition, condition)
                % Extract the relevant times
                stimulusOnsetTime = current_data.StimulusOnsetTime;
                trial_resp_time = current_data.rt;

                % Extract system timestamps and pupil sizes for both eyes
                timestamps = double(eye_tracking_data.systemTimeStamp) / 1e6;

                % Find indices for the stimulus onset
                [~, idx_stimulus_onset] = min(abs(timestamps - stimulusOnsetTime));

                % Extract pupil sizes at stimulus onset for both eyes
                pupilSize_right = eye_tracking_data.right.pupil.diameter(idx_stimulus_onset);
                pupilSize_left = eye_tracking_data.left.pupil.diameter(idx_stimulus_onset);

                % Calculate reaction time
                reactionTime = trial_resp_time;

                % Store the average pupil size and corresponding reaction time
                pupilSizes = [pupilSizes; mean([pupilSize_right, pupilSize_left])];
                reactionTimes = [reactionTimes; reactionTime];
            end
        end

        % Remove rows with NaN values
        valid_idx = ~isnan(pupilSizes) & ~isnan(reactionTimes);
        pupilSizes = pupilSizes(valid_idx);
        reactionTimes = reactionTimes(valid_idx);

        % Determine the subplot position
        subplot(2, 2, c);
        
        % Plot the pupil sizes vs. reaction times for the current condition
        scatter(pupilSizes, reactionTimes, 'filled', 'MarkerFaceColor', color_map(condition));
        title(['RT vs. Pupil Size - Condition: ', condition]);
        xlabel('Pupil Size (mm)');
        ylabel('Reaction Time (s)');
        grid on;

        % Add a regression line to the plot if there are enough data points
        if length(pupilSizes) > 1 && range(pupilSizes) > 0
            hold on;
            p = polyfit(pupilSizes, reactionTimes, 1); % Linear fit
            yfit = polyval(p, pupilSizes);
            plot(pupilSizes, yfit, 'Color', color_map(condition), 'LineWidth', 1.5);
            hold off;
        end
    end

    % Save the combined figure
    saveas(gcf, fullfile(analysis_folder, 'rt_vs_pupil_size_per_condition.png'));
end
