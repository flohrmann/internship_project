function plotPupilDiameterOverTime(id, cutData, analysis_folder)
    % Close any existing figures
    close all;

    % Define colors for each segment of the trial
    color_green = [0, 0.5, 0];        % Fixation
    color_orange = [0.8500, 0.3250, 0.0980]; % Blank
    color_light_orange = [0.9290, 0.6940, 0.1250]; % Stimulus
    color_grey = [0.5, 0.5, 0.5];     % Gap (missing data)

    % Initialize arrays for a common time base
    resampled_time = linspace(-3, 1, 1000); % Example common time base

    % Create figure with two subplots
    figure;

    % First subplot: Individual Eyes
    subplot(2,1,1);
    hold on;
    title(['Subject ', num2str(id), ': Pupil Diameter Across Trials (Individual Eyes)']);
    xlabel('Time (s)');
    ylabel('Pupil Diameter (mm)');

    % Second subplot: Average of Both Eyes
    subplot(2,1,2);
    hold on;
    title(['Subject ', num2str(id), ': Average Pupil Diameter Across Trials']);
    xlabel('Time (s)');
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

        % Resample the data to a common time base
        resampled_pupil_right = interp1(timestamps - StimulusOnsetTime, pupil_right, resampled_time, 'linear', 'extrap');
        resampled_pupil_left = interp1(timestamps - StimulusOnsetTime, pupil_left, resampled_time, 'linear', 'extrap');

        % Calculate average of both eyes
        resampled_pupil_both = (resampled_pupil_right + resampled_pupil_left) / 2;

        % Plot each segment for both eyes in the first subplot
        subplot(2,1,1);
        plot_segment(resampled_time, resampled_pupil_right, trialStartTime - StimulusOnsetTime, fixStartTime - StimulusOnsetTime, color_grey);
        plot_segment(resampled_time, resampled_pupil_right, fixStartTime - StimulusOnsetTime, blankStartTime - StimulusOnsetTime, color_green);
        plot_segment(resampled_time, resampled_pupil_right, blankStartTime - StimulusOnsetTime, StimulusOnsetTime - StimulusOnsetTime, color_orange);
        plot_segment(resampled_time, resampled_pupil_right, StimulusOnsetTime - StimulusOnsetTime, trial_resp_time - StimulusOnsetTime, color_light_orange);

        plot_segment(resampled_time, resampled_pupil_left, trialStartTime - StimulusOnsetTime, fixStartTime - StimulusOnsetTime, color_grey);
        plot_segment(resampled_time, resampled_pupil_left, fixStartTime - StimulusOnsetTime, blankStartTime - StimulusOnsetTime, color_green);
        plot_segment(resampled_time, resampled_pupil_left, blankStartTime - StimulusOnsetTime, StimulusOnsetTime - StimulusOnsetTime, color_orange);
        plot_segment(resampled_time, resampled_pupil_left, StimulusOnsetTime - StimulusOnsetTime, trial_resp_time - StimulusOnsetTime, color_light_orange);

        % Plot average in the second subplot
        subplot(2,1,2);
        plot_segment(resampled_time, resampled_pupil_both, trialStartTime - StimulusOnsetTime, fixStartTime - StimulusOnsetTime, color_grey);
        plot_segment(resampled_time, resampled_pupil_both, fixStartTime - StimulusOnsetTime, blankStartTime - StimulusOnsetTime, color_green);
        plot_segment(resampled_time, resampled_pupil_both, blankStartTime - StimulusOnsetTime, StimulusOnsetTime - StimulusOnsetTime, color_orange);
        plot_segment(resampled_time, resampled_pupil_both, StimulusOnsetTime - StimulusOnsetTime, trial_resp_time - StimulusOnsetTime, color_light_orange);
    end

    % Add legend manually
    subplot(2,1,1);
    add_legend(gca, color_grey, color_green, color_orange, color_light_orange);

    subplot(2,1,2);
    add_legend(gca, color_grey, color_green, color_orange, color_light_orange);

    hold off;
    saveas(gcf, fullfile(analysis_folder, 'pupil_diam_over_time.png'));

    %% Function to plot segments with gaps filled in grey
    function plot_segment(timestamps, data, start_time, end_time, color)
        idx_segment = timestamps >= start_time & timestamps < end_time;
        idx_gaps = timestamps >= start_time & timestamps < end_time & isnan(data);

        % Plot the data segment
        plot(timestamps(idx_segment), data(idx_segment), 'Color', color, 'LineWidth', 0.5);

        % Plot gaps in grey
        if any(idx_gaps)
            plot(timestamps(idx_gaps), data(idx_gaps), 'Color', color_grey, 'LineWidth', 0.5);
        end
    end

    %% Function to add a legend with custom colors for each segment
    function add_legend(axis_handle, color_grey, color_green, color_orange, color_light_orange)
        ylims = get(axis_handle, 'ylim');
        xlims = get(axis_handle, 'xlim');
        legend_y = ylims(2) - (ylims(2) - ylims(1)) * 0.05;
        legend_x_start = xlims(1) + (xlims(2) - xlims(1)) * 0.01;
        legend_x_end = legend_x_start + (xlims(2) - xlims(1)) * 0.05;

        plot([legend_x_start, legend_x_end], [legend_y, legend_y], 'Color', color_green, 'LineWidth', 1);
        text(legend_x_end + (xlims(2) - xlims(1)) * 0.01, legend_y, 'Fixation', 'VerticalAlignment', 'middle');

        legend_y = legend_y - (ylims(2) - ylims(1)) * 0.05;
        plot([legend_x_start, legend_x_end], [legend_y, legend_y], 'Color', color_orange, 'LineWidth', 1);
        text(legend_x_end + (xlims(2) - xlims(1)) * 0.01, legend_y, 'Blank', 'VerticalAlignment', 'middle');

        legend_y = legend_y - (ylims(2) - ylims(1)) * 0.05;
        plot([legend_x_start, legend_x_end], [legend_y, legend_y], 'Color', color_light_orange, 'LineWidth', 1);
        text(legend_x_end + (xlims(2) - xlims(1)) * 0.01, legend_y, 'Stimulus', 'VerticalAlignment', 'middle');

        legend_y = legend_y - (ylims(2) - ylims(1)) * 0.05;
        plot([legend_x_start, legend_x_end], [legend_y, legend_y], 'Color', color_grey, 'LineWidth', 1);
        text(legend_x_end + (xlims(2) - xlims(1)) * 0.01, legend_y, 'Gap', 'VerticalAlignment', 'middle');
    end
end
