function plotPupilDiameterSequential(id, cutData, analysis_folder)
    % Close any existing figures
    close all;

    % Define colors for each segment of the trial
    color_grey = [0.5, 0.5, 0.5];          % Gap (Trial Start)
    color_green = [0, 0.5, 0];              % Fixation
    color_orange = [0.8500, 0.3250, 0.0980]; % Blank
    color_light_orange = [0.9290, 0.6940, 0.1250]; % Stimulus

    % Initialize cumulative time tracker for a continuous x-axis
    cumulative_time = 0; % Tracks the cumulative time across trials

    % Create figure with two subplots
    figure;

    % First subplot: Individual Eyes
    subplot(2,1,1);
    hold on;
    title(['Subject ', num2str(id), ': Sequential Pupil Diameter Across Trials (Individual Eyes)']);
    xlabel('Time (arbitrary)');
    ylabel('Pupil Diameter (mm)');

    % Second subplot: Average of Both Eyes
    subplot(2,1,2);
    hold on;
    title(['Subject ', num2str(id), ': Sequential Average Pupil Diameter Across Trials']);
    xlabel('Time (arbitrary)');
    ylabel('Average Pupil Diameter (mm)');

    % Loop through each trial in cutData
    for trial = 1:size(cutData, 1)
        % Extract the current trial data
        current_data = cutData(trial, :);

        % Extract relevant times
        trialStartTime = current_data.trialStartTime;
        fixStartTime = current_data.fixStartTime;
        blankStartTime = current_data.blankStartTime;
        StimulusOnsetTime = current_data.StimulusOnsetTime;
        trial_resp_time = current_data.trialEndTime;

        % Extract system timestamps and pupil diameters for both eyes
        timestamps = double(current_data.eyeTrial.systemTimeStamp) / 1e6;
        pupil_right = current_data.eyeTrial.right.pupil.diameter;
        pupil_left = current_data.eyeTrial.left.pupil.diameter;

        % Resample the data to a common time base for interpolation
        resampled_time = linspace(-0.7, 1, 1000); % Example common time base
        resampled_pupil_right = interp1(timestamps - StimulusOnsetTime, pupil_right, resampled_time, 'linear', 'extrap');
        resampled_pupil_left = interp1(timestamps - StimulusOnsetTime, pupil_left, resampled_time, 'linear', 'extrap');

        % Calculate average of both eyes
        resampled_pupil_both = (resampled_pupil_right + resampled_pupil_left) / 2;

        % Calculate adjusted times for each trial segment based on cumulative time
        adjusted_trialStart = cumulative_time + (trialStartTime - StimulusOnsetTime);
        adjusted_fixStart = cumulative_time + (fixStartTime - StimulusOnsetTime);
        adjusted_blankStart = cumulative_time + (blankStartTime - StimulusOnsetTime);
        adjusted_stimOnset = cumulative_time;
        adjusted_trialEnd = cumulative_time + (trial_resp_time - StimulusOnsetTime);

        % Plot each segment for both eyes in the first subplot
        subplot(2,1,1);
        % Right Eye Segments
        plot_segment(cumulative_time + resampled_time, resampled_pupil_right, adjusted_trialStart, adjusted_fixStart, color_grey);
        plot_segment(cumulative_time + resampled_time, resampled_pupil_right, adjusted_fixStart, adjusted_blankStart, color_green);
        plot_segment(cumulative_time + resampled_time, resampled_pupil_right, adjusted_blankStart, adjusted_stimOnset, color_orange);
        plot_segment(cumulative_time + resampled_time, resampled_pupil_right, adjusted_stimOnset, adjusted_trialEnd, color_light_orange);

        % Left Eye Segments
        plot_segment(cumulative_time + resampled_time, resampled_pupil_left, adjusted_trialStart, adjusted_fixStart, color_grey);
        plot_segment(cumulative_time + resampled_time, resampled_pupil_left, adjusted_fixStart, adjusted_blankStart, color_green);
        plot_segment(cumulative_time + resampled_time, resampled_pupil_left, adjusted_blankStart, adjusted_stimOnset, color_orange);
        plot_segment(cumulative_time + resampled_time, resampled_pupil_left, adjusted_stimOnset, adjusted_trialEnd, color_light_orange);

        % Plot average in the second subplot
        subplot(2,1,2);
        plot_segment(cumulative_time + resampled_time, resampled_pupil_both, adjusted_trialStart, adjusted_fixStart, color_grey);
        plot_segment(cumulative_time + resampled_time, resampled_pupil_both, adjusted_fixStart, adjusted_blankStart, color_green);
        plot_segment(cumulative_time + resampled_time, resampled_pupil_both, adjusted_blankStart, adjusted_stimOnset, color_orange);
        plot_segment(cumulative_time + resampled_time, resampled_pupil_both, adjusted_stimOnset, adjusted_trialEnd, color_light_orange);

        % Update cumulative time by adding the duration of the current trial
        cumulative_time = adjusted_trialEnd + 0.5; % Add a small gap between trials
    end

    % Add legend manually to the first subplot
    subplot(2,1,1);
    ylims = get(gca, 'ylim');
    xlims = get(gca, 'xlim');
    legend_y = ylims(2) - (ylims(2) - ylims(1)) * 0.05;
    legend_x_start = xlims(1) + (xlims(2) - xlims(1)) * 0.01;
    legend_x_end = legend_x_start + (xlims(2) - xlims(1)) * 0.05;

    % Plot legend lines
    plot([legend_x_start, legend_x_end], [legend_y, legend_y], 'Color', color_grey, 'LineWidth', 2);
    text(legend_x_end + (xlims(2) - xlims(1)) * 0.01, legend_y, 'Gap', 'VerticalAlignment', 'middle');

    legend_y = legend_y - (ylims(2) - ylims(1)) * 0.05;
    plot([legend_x_start, legend_x_end], [legend_y, legend_y], 'Color', color_green, 'LineWidth', 2);
    text(legend_x_end + (xlims(2) - xlims(1)) * 0.01, legend_y, 'Fixation', 'VerticalAlignment', 'middle');

    legend_y = legend_y - (ylims(2) - ylims(1)) * 0.05;
    plot([legend_x_start, legend_x_end], [legend_y, legend_y], 'Color', color_orange, 'LineWidth', 2);
    text(legend_x_end + (xlims(2) - xlims(1)) * 0.01, legend_y, 'Blank', 'VerticalAlignment', 'middle');

    legend_y = legend_y - (ylims(2) - ylims(1)) * 0.05;
    plot([legend_x_start, legend_x_end], [legend_y, legend_y], 'Color', color_light_orange, 'LineWidth', 2);
    text(legend_x_end + (xlims(2) - xlims(1)) * 0.01, legend_y, 'Stimulus', 'VerticalAlignment', 'middle');

    % Add legend manually to the second subplot
    subplot(2,1,2);
    ylims = get(gca, 'ylim');
    xlims = get(gca, 'xlim');
    legend_y = ylims(2) - (ylims(2) - ylims(1)) * 0.05;
    legend_x_start = xlims(1) + (xlims(2) - xlims(1)) * 0.01;
    legend_x_end = legend_x_start + (xlims(2) - xlims(1)) * 0.05;

    % Plot legend lines
    plot([legend_x_start, legend_x_end], [legend_y, legend_y], 'Color', color_grey, 'LineWidth', 2);
    text(legend_x_end + (xlims(2) - xlims(1)) * 0.01, legend_y, 'Gap', 'VerticalAlignment', 'middle');

    legend_y = legend_y - (ylims(2) - ylims(1)) * 0.05;
    plot([legend_x_start, legend_x_end], [legend_y, legend_y], 'Color', color_green, 'LineWidth', 2);
    text(legend_x_end + (xlims(2) - xlims(1)) * 0.01, legend_y, 'Fixation', 'VerticalAlignment', 'middle');

    legend_y = legend_y - (ylims(2) - ylims(1)) * 0.05;
    plot([legend_x_start, legend_x_end], [legend_y, legend_y], 'Color', color_orange, 'LineWidth', 2);
    text(legend_x_end + (xlims(2) - xlims(1)) * 0.01, legend_y, 'Blank', 'VerticalAlignment', 'middle');

    legend_y = legend_y - (ylims(2) - ylims(1)) * 0.05;
    plot([legend_x_start, legend_x_end], [legend_y, legend_y], 'Color', color_light_orange, 'LineWidth', 2);
    text(legend_x_end + (xlims(2) - xlims(1)) * 0.01, legend_y, 'Stimulus', 'VerticalAlignment', 'middle');

    hold off;
    saveas(gcf, fullfile(analysis_folder, 'pupil_diam_sequential.png'));

    %% Nested function to plot segments with gaps filled in grey
    function plot_segment(timestamps, data, start_time, end_time, color)
        idx_segment = timestamps >= start_time & timestamps < end_time;
        idx_gaps = timestamps >= start_time & timestamps < end_time & isnan(data);

        % Plot the data segment
        plot(timestamps(idx_segment), data(idx_segment), 'Color', color, 'LineWidth', 1.5);

        % Plot gaps in grey
        if any(idx_gaps)
            plot(timestamps(idx_gaps), data(idx_gaps), 'Color', color_grey, 'LineWidth', 1.5);
        end
    end
end
