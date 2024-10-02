function plotNormalizedPupilDiameterByCondition(id, cutData, conditions, analysis_folder)
    % Close any existing figures
    % sliding average statt interpolating
    % average both eyes
    close all;

    % Define colors for plotting
    color_blue = [0, 0.4470, 0.7410]; % Blue for fixation
    color_green = [0, 0.5, 0]; % Green for blank
    color_orange = [0.9, 0.4, 0.2]; % Orange for stimulus
    color_grey = [0.5, 0.5, 0.5]; % Grey for gaps

    % Define consistent lengths for fixation, blank, and stimulus segments
    num_fixation_points = 31;  % Normalize fixation to 30 points
    num_blank_points = 13;     % Normalize blank to 12 points
    num_stimulus_points = 100;  % Interpolate stimulus to 30 points

    % Create a new figure for the subplots
    figure;

    % Initialize arrays to hold segmented data for right and left pupil for each condition
    fixation_pupil_right = [];
    blank_pupil_right = [];
    stim_pupil_right = [];

    fixation_pupil_left = [];
    blank_pupil_left = [];
    stim_pupil_left = [];

    % tiled plot layout
    t = tiledlayout('flow','TileSpacing','compact');
            
    % Loop through each condition to create subplots
    for condIdx = 1:length(conditions)
        condition = conditions{condIdx};
        
        % Initialize segmented data per condition
        fixation_pupil_right = [];
        blank_pupil_right = [];
        stim_pupil_right = [];
    
        fixation_pupil_left = [];
        blank_pupil_left = [];
        stim_pupil_left = [];
        
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
        
                % Interpolate and normalize the fixation, blank, and stimulus periods
        
                % Fixation (normalize to 30 points, from -0.7s to -0.2s)
                fixation_time = linspace(-0.7, -0.2, idx_blank - idx_fixation +1); % Fixation happens from -0.7s to -0.2s
                if length(fixation_time) == (idx_blank - idx_fixation +1)
                    fixation_pupil_right_interp = interp1(fixation_time, pupil_right(idx_fixation:idx_blank), linspace(-0.7, -0.2, num_fixation_points), 'linear', 'extrap');
                    fixation_pupil_left_interp  = interp1(fixation_time, pupil_left(idx_fixation:idx_blank), linspace(-0.7, -0.2, num_fixation_points), 'linear', 'extrap');
                else
                    fixation_pupil_right_interp = nan(1, num_fixation_points);
                    fixation_pupil_left_interp = nan(1, num_fixation_points);
                end
        
                % Blank (normalize to 12 points, from -0.2s to 0s)
                blank_time = linspace(-0.2, 0, idx_stimulus - idx_blank +1); % Blank happens from -0.2s to 0s
                if length(blank_time) == (idx_stimulus - idx_blank +1)
                    blank_pupil_right_interp = interp1(blank_time, pupil_right(idx_blank:idx_stimulus), linspace(-0.2, 0, num_blank_points), 'linear', 'extrap');
                    blank_pupil_left_interp = interp1(blank_time, pupil_left(idx_blank:idx_stimulus), linspace(-0.2, 0, num_blank_points), 'linear', 'extrap');
                else
                    blank_pupil_right_interp = nan(1, num_blank_points);
                    blank_pupil_left_interp = nan(1, num_blank_points);
                end
        
                % Stimulus (interpolate to 100 points)
                stim_pupil_right_interp = interp1(1:length(pupil_right(idx_stimulus:idx_end)), pupil_right(idx_stimulus:idx_end), linspace(1, length(pupil_right(idx_stimulus:idx_end)), num_stimulus_points), 'linear', 'extrap');
                stim_pupil_left_interp = interp1(1:length(pupil_left(idx_stimulus:idx_end)), pupil_left(idx_stimulus:idx_end), linspace(1, length(pupil_left(idx_stimulus:idx_end)), num_stimulus_points), 'linear', 'extrap');
        
                % Append interpolated data for the current condition
                fixation_pupil_right = [fixation_pupil_right; fixation_pupil_right_interp];
                blank_pupil_right = [blank_pupil_right; blank_pupil_right_interp];
                stim_pupil_right = [stim_pupil_right; stim_pupil_right_interp];
        
                fixation_pupil_left = [fixation_pupil_left; fixation_pupil_left_interp];
                blank_pupil_left = [blank_pupil_left; blank_pupil_left_interp];
                stim_pupil_left = [stim_pupil_left; stim_pupil_left_interp];
            end
        end
    
        % Calculate the mean pupil diameter for each segment
        mean_fixation_right = nanmean(fixation_pupil_right, 1);
        mean_blank_right = nanmean(blank_pupil_right, 1);
        mean_stim_right = nanmean(stim_pupil_right, 1);
    
        mean_fixation_left = nanmean(fixation_pupil_left, 1);
        mean_blank_left = nanmean(blank_pupil_left, 1);
        mean_stim_left = nanmean(stim_pupil_left, 1);
    
        %% Create a subplot for each condition
        %subplot(2, 2, condIdx); % 2 rows, 2 columns
        
        nexttile
        hold on;
        
        title(['Condition: ', condition]);
        ylabel('Pupil Diameter (mm)');
        xlabel('Time Normalized to -0.7:1');
 
        % Plot the mean pupil diameter with the appropriate colors
        plot(linspace(-0.7, -0.2, num_fixation_points), mean_fixation_right, 'Color', color_blue, 'LineWidth', 1.5);
        plot(linspace(-0.2, 0, num_blank_points), mean_blank_right, 'Color', color_green, 'LineWidth', 1.5);
        plot(linspace(0, 1, num_stimulus_points), mean_stim_right, 'Color', color_orange, 'LineWidth', 1.5);
    
        
        % Add the "Raw Data" and "Normalized Data" text for the whole plot
%         xlims = get(gca, 'xlim');
%         ylims = get(gca, 'ylim');
%         annotation('textbox', [0.08, 0.6, 0.5, 0.05], 'String', 'Raw Data [s]', 'LineStyle', 'none'); % 'HorizontalAlignment', 'center',
%         annotation('textbox', [0.75, 0.05, 0.5, 0.05], 'String', 'Normalized Data [0:1]', 'HorizontalAlignment', 'center', 'LineStyle', 'none');


        plot(linspace(-0.7, -0.2, num_fixation_points), mean_fixation_left, '--', 'Color', color_blue, 'LineWidth', 1.5);
        plot(linspace(-0.2, 0, num_blank_points), mean_blank_left, '--', 'Color', color_green, 'LineWidth', 1.5);
        plot(linspace(0, 1, num_stimulus_points), mean_stim_left, '--', 'Color', color_orange, 'LineWidth', 1.5);
    
        % Add a vertical line at x = 0
        xline(0, '--', 'LineWidth', 1.1, 'Color', color_grey); % Dashed black line at x=0
        % hidden plot just for legend
        h = plot([0 0], ylim, '--', 'Color', color_grey, 'LineWidth', 1.5, 'Visible', 'off');
    end
    hold off;
    lgd = legend({'Fixation (Right)', 'Blank (Right)', 'Stimulus (Right)', ...
            'Fixation (Left)', 'Blank (Left)', 'Stimulus (Left)', 'Stimulation Onset'}, ...
            'Orientation', 'horizontal'); 
    lgd.NumColumns = 1;
    lgd.Layout.Tile = 5;
    lgd.Location = 'bestoutside';
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1])
    saveas(gcf, fullfile(analysis_folder, 'pupil_diam_combined.png'));
end
