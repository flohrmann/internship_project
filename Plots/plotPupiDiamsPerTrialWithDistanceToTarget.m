function plotPupiDiamsPerTrialWithDistanceToTarget(cutData, diam_t0, screenXpixels, screenYpixels, tolerance, num_before, num_after, baseline_length, analysis_folder)
% Function to visualize different steps of data processing in findTargetAndExtractData
for i = 1:size(cutData, 1)
    
    % Select a sample trial to visualize each step
    current_trial = cutData(i, :);
    current_diam_start  = diam_t0.DiamBothStartOutliersRemoved(i,:);
    current_diam_blank  = diam_t0.DiamBothBlankOutliersRemoved(i,:);
    
    timestamps_eyeTrial = double(current_trial.eyeTrial.systemTimeStamp) / 1e6;  % Convert eyeTrial timestamps to seconds
    timestamps_stimulusTrial = double(current_trial.stimulusTrial.systemTimeStamp) / 1e6;  % Convert stimulusTrial timestamps to seconds
    
    figure;
    
    % if all elements are nans skip trial
    %% Step 1: Plot Baseline from Blank Screen
    if  numel(current_diam_blank(~isnan(current_diam_blank))) == size(current_diam_blank, 2)
        subplot(2, 3, 1);
        %         blank_d_avg = nanmean([sample_trial.eyeTrial.left.pupil.diameter; sample_trial.eyeTrial.right.pupil.diameter], 1);
        blank_start = current_trial.blankStartTime;
        [~, idx_blank] = min(abs(timestamps_eyeTrial - blank_start));
        %baseline_data = blank_d_avg(idx_blank:idx_blank + baseline_length);
        %         plot(timestamps_eyeTrial(idx_blank:idx_blank + baseline_length), baseline_data, 'b');
        plot(timestamps_eyeTrial(idx_blank:idx_blank + baseline_length-1), current_diam_blank, 'b');
        hold on;
        yline(nanmean(current_diam_blank), 'r--', 'Mean');  % Mean line
        title('Baseline from Blank Screen');
        xlabel('Time (s)');
        ylabel('Pupil Diameter');
        hold off;
    else
        fprintf('Error processing trial %d. Baseline was previously discarded.\n', i);
        continue;
    end
    
    %% Step 2: Plot Baseline from Trial Start (Stimulation Screen)   
    if numel(current_diam_start(~isnan(current_diam_start))) == size(current_diam_start, 2)        
        subplot(2, 3, 2);
        %start_d_avg = nanmean([current_trial.eyeTrial.left.pupil.diameter; current_trial.eyeTrial.right.pupil.diameter], 1);
        %baseline_start_data = start_d_avg(1:baseline_length);
        plot(timestamps_eyeTrial(1:baseline_length), current_diam_start, 'g');
        hold on;
        yline(nanmean(current_diam_start), 'r--', 'Mean');  % Mean line
        title('Baseline from Start of Trial');
        xlabel('Time (s)');
        ylabel('Pupil Diameter');
        hold off;
    else
        fprintf('Error processing trial %d. Baseline was previously discarded.\n', i);
        continue;
    end
    try
        % Common y-axis
        ylim_values = [min([current_diam_blank, current_diam_start]), max([current_diam_blank, current_diam_start])];
        subplot(2, 3, 1); ylim(ylim_values);
        subplot(2, 3, 2); ylim(ylim_values);
    catch % if not possible bc only one plotted ignore
    end
    
    
    try
        
        %% Step 3: Target Position on Screen with Initial Gaze Position
        subplot(2, 3, 3);
        targetPos = current_trial.TargetPosition;
        targetX = current_trial.x_centers{1}(targetPos(2), targetPos(1));
        targetY = screenYpixels - current_trial.y_centers{1}(targetPos(2), targetPos(1));
        plot(targetX, targetY, 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Target');
        hold on;
        
        % All gaze positions
        x_r = current_trial.stimulusTrial.right.gazePoint.onDisplayArea(1, :) * screenXpixels;
        y_r = screenYpixels - current_trial.stimulusTrial.right.gazePoint.onDisplayArea(2, :) * screenYpixels;
        x_l = current_trial.stimulusTrial.left.gazePoint.onDisplayArea(1, :) * screenXpixels;
        y_l = screenYpixels - current_trial.stimulusTrial.left.gazePoint.onDisplayArea(2, :) * screenYpixels;
        scatter(x_r, y_r, 5, 'r', 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', 'Right Eye');
        scatter(x_l, y_l, 5, 'b', 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', 'Left Eye');
        
        % Initial gaze positions at start of trial
        initial_x_r = current_trial.stimulusTrial.right.gazePoint.onDisplayArea(1, 1) * screenXpixels;
        initial_y_r = screenYpixels - current_trial.stimulusTrial.right.gazePoint.onDisplayArea(2, 1) * screenYpixels;
        initial_x_l = current_trial.stimulusTrial.left.gazePoint.onDisplayArea(1, 1) * screenXpixels;
        initial_y_l = screenYpixels - current_trial.stimulusTrial.left.gazePoint.onDisplayArea(2, 1) * screenYpixels;
        plot(initial_x_r, initial_y_r, 'ro', 'DisplayName', 'Initial Right Eye Gaze');
        plot(initial_x_l, initial_y_l, 'bo', 'DisplayName', 'Initial Left Eye Gaze');
        legend;
        legend boxoff;
        legend('Location', 'best');
        xlim([0 screenXpixels]);
        ylim([0 screenYpixels]);
        title('Target Position and Initial Gaze');
        xlabel('X Position');
        ylabel('Y Position');
        hold off;
        
        %% Step 4: Pupil Diameter and Distance to Target Over Whole Trial
        subplot(2, 3, [4 6]);
        hold on;
        
        % Define stimulation start time for time alignment
        stimulation_start_time = timestamps_stimulusTrial(1);
        
        % Calculate target found index
        x_avg = nanmean([x_r; x_l], 1);
        y_avg = nanmean([y_r; y_l], 1);
        distance_to_target = sqrt((x_avg - targetX).^2 + (y_avg - targetY).^2);
        found_target_idx = find(distance_to_target <= tolerance, 1, 'first');
        found_target_time = timestamps_stimulusTrial(found_target_idx);
        
        % Adjust timestamps relative to stimulation start
        adjusted_timestamps_eyeTrial = timestamps_eyeTrial - stimulation_start_time;
        adjusted_timestamps_stimulusTrial = timestamps_stimulusTrial - stimulation_start_time;
        found_target_time_adjusted = found_target_time - stimulation_start_time;
        
        
        % Define the pupil diameter data for eyeTrial and stimulusTrial
        pupil_diam_avg_eyeTrial = nanmean([current_trial.eyeTrial.left.pupil.diameter; current_trial.eyeTrial.right.pupil.diameter], 1);
        pupil_diam_avg_stimulusTrial = nanmean([current_trial.stimulusTrial.left.pupil.diameter; current_trial.stimulusTrial.right.pupil.diameter], 1);
        
        % Plot EyeTrial data before stimulation (entire period)
        plot(adjusted_timestamps_eyeTrial, pupil_diam_avg_eyeTrial, 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, 'DisplayName', 'EyeTrial');
        
        % Plot StimulusTrial data with different colors for segments
        % Segment before the target was found (light green)
        plot(adjusted_timestamps_stimulusTrial(max(1, found_target_idx - num_before):found_target_idx), ...
            pupil_diam_avg_stimulusTrial(max(1, found_target_idx - num_before):found_target_idx), ...
            'Color', [0.6 1 0.6], 'LineWidth', 1.5, 'DisplayName', 'Before Target');
        
        % Segment after the target was found (dark green)
        plot(adjusted_timestamps_stimulusTrial(found_target_idx:min(found_target_idx + num_after, length(adjusted_timestamps_stimulusTrial))), ...
            pupil_diam_avg_stimulusTrial(found_target_idx:min(found_target_idx + num_after, length(adjusted_timestamps_stimulusTrial))), ...
            'Color', [0.1 0.5 0.1], 'LineWidth', 1.5, 'DisplayName', 'After Target');
        
        % Plot remaining StimulusTrial data in dark grey (outside defined processing range)
        if found_target_idx + num_after < length(adjusted_timestamps_stimulusTrial)
            plot(adjusted_timestamps_stimulusTrial(found_target_idx + num_after:end), ...
                pupil_diam_avg_stimulusTrial(found_target_idx + num_after:end), ...
                'Color', [0.3 0.3 0.3], 'LineWidth', 1.5, 'DisplayName', 'Remaining Data');
        end
        
        % Overlay Distance to Target in Orange
        yyaxis right;  % Use the right y-axis for the distance to target
        plot(adjusted_timestamps_stimulusTrial(1:length(distance_to_target)), distance_to_target, 'Color', [1 0.5 0], 'LineWidth', 1.5, 'DisplayName', 'Distance to Target');
        ylabel('Distance to Target (pixels)');
        set(gca, 'YColor', [1 0.5 0]);  % Set the right y-axis color to orange
        yyaxis left;  % Return to the left y-axis for the rest of the data
        
        % Add lines for when the target was found and blank screen start, with labels on the right, excluding from legend
        xline(blank_start - stimulation_start_time, '--', 'Blank Screen', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
        xline(current_trial.fixStartTime - stimulation_start_time, '--', 'Fixation Screen', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
        xline(current_trial.StimulusOnsetTime - stimulation_start_time, '--', 'Stimulation Screen', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
        xline(found_target_time_adjusted, '--', 'Gaze at Target', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
        xline(current_trial.rt, '--', 'Button Pressed', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
        
        title('Pupil Diameter and Distance to Target Over Whole Trial');
        xlabel('Time from Stimulation Start (s)');
        ylabel('Pupil Diameter');
        legend;
        legend boxoff;
        legend('Location', 'bestoutside');
        grid off;
        hold off;
        
        % Set figure properties and save the figure
        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf, fullfile(analysis_folder, sprintf('single_trial_baseline_pupil_diam_distance_target_%02d.png', i)));
        
        % Close the figure after saving
        close(gcf);
        
    catch
        % Close all open figures in case of error
        close all;
        fprintf('Error processing trial %d. Skipping to next.\n', i);
    end
end
end
