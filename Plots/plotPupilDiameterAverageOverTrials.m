function plotPupilDiameterAverageOverTrials(id, cutData, analysis_folder)
% Close any existing figures
close all;

color_blue = [0, 0.4470, 0.7410]; % Blue
color_green = [0, 0.5, 0]; % Green
color_orange = [0.8500, 0.3250, 0.0980]; % Orange
color_light_orange = [0.9290, 0.6940, 0.1250];
color_grey = [0.5, 0.5, 0.5]; % Grey

% Initialize arrays to hold resampled data
resampled_time = linspace(-0.7, 1, 1000); % Example common time base
resampled_pupil_right = nan(size(cutData, 1), length(resampled_time));
resampled_pupil_left = nan(size(cutData, 1), length(resampled_time));

% Loop through each trial in cutData
for trial = 1:size(cutData, 1)
    % Extract the current trial data
    current_data = cutData(trial, :);
    
    % Extract the relevant times
    %trialStartTime = current_data.StimulusOnsetTime; % todo get real data
    fixStartTime = current_data.fixStartTime;
    blankStartTime = current_data.blankStartTime;
    StimulusOnsetTime = current_data.StimulusOnsetTime;
    trial_resp_time = current_data.trialEndTime;
    
    % Extract system timestamps and pupil diameters for both eyes
    timestamps = double(current_data.eyeTrial.systemTimeStamp) / 1e6;
    pupil_right = current_data.eyeTrial.right.pupil.diameter;
    pupil_left = current_data.eyeTrial.left.pupil.diameter;
    
    % Resample the data to the common time base
    resampled_pupil_right(trial, :) = interp1(timestamps - StimulusOnsetTime, pupil_right, resampled_time, 'linear', 'extrap');
    resampled_pupil_left(trial, :) = interp1(timestamps - StimulusOnsetTime, pupil_left, resampled_time, 'linear', 'extrap');
end

% Calculate the mean pupil diameter across trials
mean_pupil_right = nanmean(resampled_pupil_right, 1);
mean_pupil_left = nanmean(resampled_pupil_left, 1);

% Create a new figure
figure;
hold on;
title(['Subject ', num2str(id),': Mean Pupil Diameter Over Time']);
%xlabel('Time (s)');
ylabel('Left and Right Eye: Pupil Diameter (mm)');

% Plot the mean pupil diameter with the appropriate colors
%plot_segment(resampled_time, mean_pupil_right, -1, fixStartTime - StimulusOnsetTime, color_blue);
plot_segment(resampled_time, mean_pupil_right, fixStartTime - StimulusOnsetTime, blankStartTime - StimulusOnsetTime, color_green);
plot_segment(resampled_time, mean_pupil_right, blankStartTime - StimulusOnsetTime, StimulusOnsetTime - StimulusOnsetTime, color_orange);
plot_segment(resampled_time, mean_pupil_right, StimulusOnsetTime - StimulusOnsetTime, trial_resp_time - StimulusOnsetTime, color_light_orange);

%plot_segment(resampled_time, mean_pupil_left, -1, fixStartTime - StimulusOnsetTime, color_blue);
plot_segment(resampled_time, mean_pupil_left, fixStartTime - StimulusOnsetTime, blankStartTime - StimulusOnsetTime, color_green);
plot_segment(resampled_time, mean_pupil_left, blankStartTime - StimulusOnsetTime, StimulusOnsetTime - StimulusOnsetTime, color_orange);
plot_segment(resampled_time, mean_pupil_left, StimulusOnsetTime - StimulusOnsetTime, trial_resp_time - StimulusOnsetTime, color_light_orange);

% handmade legend
ylims = get(gca, 'ylim');
xlims = get(gca, 'xlim');
legend_y = ylims(2) - (ylims(2) - ylims(1)) * 0.05;
legend_x_start = xlims(1) + (xlims(2) - xlims(1)) * 0.01;
legend_x_end = legend_x_start + (xlims(2) - xlims(1)) * 0.05;

% plot([legend_x_start, legend_x_end], [legend_y, legend_y], 'Color', color_blue, 'LineWidth', 2);
% text(legend_x_end + (xlims(2) - xlims(1)) * 0.01, legend_y, 'Press Any Button', 'VerticalAlignment', 'middle');

legend_y = legend_y - (ylims(2) - ylims(1)) * 0.05;
plot([legend_x_start, legend_x_end], [legend_y, legend_y], 'Color', color_green, 'LineWidth', 2);
text(legend_x_end + (xlims(2) - xlims(1)) * 0.01, legend_y, 'Fixation', 'VerticalAlignment', 'middle');

legend_y = legend_y - (ylims(2) - ylims(1)) * 0.05;
plot([legend_x_start, legend_x_end], [legend_y, legend_y], 'Color', color_orange, 'LineWidth', 2);
text(legend_x_end + (xlims(2) - xlims(1)) * 0.01, legend_y, 'Blank', 'VerticalAlignment', 'middle');

legend_y = legend_y - (ylims(2) - ylims(1)) * 0.05;
plot([legend_x_start, legend_x_end], [legend_y, legend_y], 'Color', color_light_orange, 'LineWidth', 2);
text(legend_x_end + (xlims(2) - xlims(1)) * 0.01, legend_y, 'Stimulus', 'VerticalAlignment', 'middle');

legend_y = legend_y - (ylims(2) - ylims(1)) * 0.05;
plot([legend_x_start, legend_x_end], [legend_y, legend_y], 'Color', color_grey, 'LineWidth', 2);
text(legend_x_end + (xlims(2) - xlims(1)) * 0.01, legend_y, 'Gap', 'VerticalAlignment', 'middle');

%% change labels (unused)
% % Modify the x-axis: before 0 is "Time (s)", after 0 is normalized to [0, 1]
% x_ticks = linspace(-0.7, 1, 17); % Define x-axis ticks (adjust as needed)
% x_labels = arrayfun(@(x) num2str(x, '%.1f'), x_ticks, 'UniformOutput', false);
% % Change labels after 0 to [0, 1]
% for i = 1:length(x_ticks)
%     if x_ticks(i) >= 0
%         % Convert the time after 0 to normalized values
%         x_labels{i} = num2str((x_ticks(i) - 0), '%.1f');
%     end
% end
% set(gca, 'XTick', x_ticks, 'XTickLabel', x_labels);

% Add the "Raw Data" and "Normalized Data" text below the x-axis
% Get axis limits for dynamic placement
xlims = get(gca, 'xlim');
ylims = get(gca, 'ylim');
xline(0, '--', 'LineWidth', 1.1, 'Color', 'k'); % 'k' for black color

% Position the text right below the x-axis based on the axis limits
text(-0.4, ylims(1) - 0.07 * (ylims(2) - ylims(1)), 'Raw Data [s]', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(0.5, ylims(1) - 0.07 * (ylims(2) - ylims(1)), 'Normalized Data [0:1]', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');




hold off;
saveas(gcf,strcat(analysis_folder, '\pupil_diam_averaged.png'));


%% Function to plot segments with gaps filled with grey
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
