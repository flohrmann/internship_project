function plotAvgNormalizedPupilDiameterByCondition(id, cutData, conditions, analysis_folder)
    % Close any existing figures
    close all;

    % Define colors for plotting
    color_blue = [0, 0.4470, 0.7410]; % Blue for fixation
    color_green = [0, 0.5, 0]; % Green for blank
    color_orange = [0.9, 0.4, 0.2]; % Orange for stimulus
    color_grey = [0.5, 0.5, 0.5]; % Grey for gaps

    % Define consistent lengths for fixation, blank, and stimulus segments
    num_fixation_points = 30;  % Fixation to 30 points
    num_blank_points = 12;     % Blank to 12 points
    num_stimulus_points = 100;  % Stimulus to 100 points

    % Create a new figure for the subplots
    figure;

    % tiled plot layout
    t = tiledlayout('flow','TileSpacing','compact');
    title(t, ['Average Pupil Diameter Over Trials - Participant ', num2str(id)], 'FontSize', 14, 'FontWeight', 'bold');
    % Define sliding window size for smoothing
    window_size = 3; % Adjust this value as needed

    % Loop through each condition to create subplots
    for condIdx = 1:length(conditions)
        condition = conditions{condIdx};

        % Initialize segmented data per condition
        fixation_pupil = [];
        blank_pupil = [];
        stim_pupil = [];

        % Loop through each trial in cutData and select based on condition
        for trial = 1:size(cutData, 1)
            if strcmp(cutData.Condition{trial}, condition) % Check condition match
                % Extract the current trial data
                current_data = cutData(trial, :);

                % Extract the relevant times
                fixStartTime = current_data.fixStartTime;
                blankStartTime = current_data.blankStartTime;
                StimulusOnsetTime = current_data.StimulusOnsetTime;
                trialEndTime = current_data.trialEndTime;

                % Extract system timestamps and pupil diameters for both eyes
                timestamps = double(current_data.eyeTrial.systemTimeStamp) / 1e6;
                pupil_right = current_data.eyeTrial.right.pupil.diameter;
                pupil_left = current_data.eyeTrial.left.pupil.diameter;

                % Get the indices corresponding to the key time points
                [~, idx_fixation] = min(abs(timestamps - fixStartTime));
                [~, idx_blank] = min(abs(timestamps - blankStartTime));
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

                % Resample each segment to ensure consistent lengths

                % Fixation: resample to 30 points
                if idx_blank > idx_fixation
                    fixation_pupil_interp = interp1(linspace(-0.7, -0.2, idx_blank - idx_fixation), pupil_avg(idx_fixation:idx_blank-1), linspace(-0.7, -0.2, num_fixation_points), 'linear', 'extrap');
                else
                    fixation_pupil_interp = nan(1, num_fixation_points);
                end

                % Blank: resample to 12 points
                if idx_stimulus > idx_blank
                    blank_pupil_interp = interp1(linspace(-0.2, 0, idx_stimulus - idx_blank), pupil_avg(idx_blank:idx_stimulus-1), linspace(-0.2, 0, num_blank_points), 'linear', 'extrap');
                else
                    blank_pupil_interp = nan(1, num_blank_points);
                end

                % Stimulus: resample to 100 points
                stim_pupil_interp = interp1(1:length(pupil_avg(idx_stimulus:idx_end)), pupil_avg(idx_stimulus:idx_end), linspace(1, length(pupil_avg(idx_stimulus:idx_end)), num_stimulus_points), 'linear', 'extrap');

                % Append resampled data for the current condition
                fixation_pupil = [fixation_pupil; fixation_pupil_interp];
                blank_pupil = [blank_pupil; blank_pupil_interp];
                stim_pupil = [stim_pupil; stim_pupil_interp];
            end
        end

        % Calculate the mean pupil diameter for each segment
        mean_fixation = nanmean(fixation_pupil, 1);
        mean_blank = nanmean(blank_pupil, 1);
        mean_stim = nanmean(stim_pupil, 1);

        %% Create a subplot for each condition
        nexttile
        hold on;

        title(['Condition: ', condition]);
        ylabel('Average Pupil Diameter (mm)');
        xlabel('Time Normalized to -0.7:1');

        % Plot the mean pupil diameter with the appropriate colors
        plot(linspace(-0.7, -0.2, num_fixation_points), mean_fixation, 'Color', color_blue, 'LineWidth', 1.5);
        plot(linspace(-0.2, 0, num_blank_points), mean_blank, 'Color', color_green, 'LineWidth', 1.5);
        plot(linspace(0, 1, num_stimulus_points), mean_stim, 'Color', color_orange, 'LineWidth', 1.5);

        % Add a vertical line at x = 0
        xline(0, '--', 'LineWidth', 1.1, 'Color', color_grey); % Dashed black line at x=0
        % hidden plot just for legend
        h = plot([0 0], ylim, '--', 'Color', color_grey, 'LineWidth', 1.5, 'Visible', 'off');
    end
    hold off;
    %% Add the legend and global title outside the subplots
    lgd = legend({'Fixation Screen', 'Blank Screen', 'Stimulation Screen', 'Stimulation Onset'}, ...
        'Orientation', 'horizontal');
    lgd.NumColumns = 1;
    lgd.Layout.Tile = 5;
    lgd.Location = 'bestoutside';

    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1])
    saveas(gcf, fullfile(analysis_folder, 'avg_pupil_diam_over_trials_per_condition.png'));
end
