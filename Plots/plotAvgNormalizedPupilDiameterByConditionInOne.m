function avg_results = plotAvgNormalizedPupilDiameterByConditionInOne(id, cutData, conditions, analysis_folder, color_map)
    % Close any existing figures
    close all;

    % Define consistent length for the stimulus segment
    num_stimulus_points = 100;  % Stimulus to 100 points

    % Initialize structure to hold averaged results
    avg_results = struct();

    % Create a new figure for the plot
    figure;
    hold on;

    % Define sliding window size for smoothing
    window_size = 5; % Adjust this value as needed

    % Loop through each condition to plot all in one figure
    for condIdx = 1:length(conditions)
        condition = conditions{condIdx};
        
        % Initialize stimulus data per condition
        stim_pupil_avg = [];
        
        % Loop through each trial in cutData and select based on condition
        for trial = 1:size(cutData, 1)
            if strcmp(cutData.Condition{trial}, condition) % Check condition match
                % Extract the current trial data
                current_data = cutData(trial, :);
                
                % Extract the relevant times
                StimulusOnsetTime = current_data.StimulusOnsetTime;
                trialEndTime = current_data.trialEndTime;
                
                % Extract system timestamps and pupil diameters for both eyes
                timestamps = double(current_data.eyeTrial.systemTimeStamp) / 1e6;
                pupil_right = current_data.eyeTrial.right.pupil.diameter;
                pupil_left = current_data.eyeTrial.left.pupil.diameter;
                
                % Get the indices corresponding to the key time points
                [~, idx_stimulus] = min(abs(timestamps - StimulusOnsetTime));
                [~, idx_end] = min(abs(timestamps - trialEndTime));
                
                % Apply sliding average for smoothing
                pupil_right_smooth = movmean(pupil_right, window_size);
                pupil_left_smooth = movmean(pupil_left, window_size);
                
                % Average the two eyes
                pupil_avg = nanmean([pupil_right_smooth, pupil_left_smooth], 1);
                % Handle NaN by taking the mean of the previous and next values
                for i = 2:length(pupil_avg)-1
                    if isnan(pupil_avg(i))
                        pupil_avg(i) = nanmean([pupil_avg(i-1), pupil_avg(i+1)]);
                    end
                end
                
                % Resample the stimulus segment to ensure consistent length
                stim_pupil_interp = interp1(1:length(pupil_avg(idx_stimulus:idx_end)), ...
                    pupil_avg(idx_stimulus:idx_end), linspace(1, length(pupil_avg(idx_stimulus:idx_end)), num_stimulus_points), 'linear', 'extrap');
                
                % Append resampled data for the current condition
                stim_pupil_avg = [stim_pupil_avg; stim_pupil_interp];
            end
        end
        
        % Calculate the mean pupil diameter for the stimulus segment
        mean_stim = nanmean(stim_pupil_avg, 1);
        
        % Plot the mean pupil diameter for the stimulus segment
        plot(linspace(0, 1, num_stimulus_points), mean_stim, 'Color', color_map(condition), 'LineWidth', 1.5, 'DisplayName', condition);
        
        % Store the averaged stimulus data in the structure
        avg_results.(condition) = mean_stim;
    end
    xlabel('Trial Duration Normalized to [0:1]');
    ylabel('Average Pupil Diameter (mm)');
    title(['Average Pupil Diameter During Stimulus Presentation - Participant ', num2str(id)]);
    legend('Location', 'best');
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1])
    saveas(gcf, fullfile(analysis_folder, 'avg_pupil_diam_during_stim.png'));
    save(fullfile(analysis_folder, 'avg_pupil_diam_during_stim_results.mat'), 'avg_results');

    hold off;
end
