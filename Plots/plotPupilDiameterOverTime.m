function plotPupilDiameterOverTime(id, cutData, eye_tracking_data, trial_results, analysis_folder)
% Close any existing figures
close all;

% Create a new figure
figure;
hold on;
title(['Subject ', num2str(id),': Pupil Diameter Over Time']);
%title('Pupil Diameter Over Time');
xlabel('Time (s)');
ylabel('Left and Right Eye: Pupil Diameter (mm)');

% Define colors for plotting
color_blue = [0, 0.4470, 0.7410]; % Blue
color_green = [0, 0.5, 0]; % Green
color_orange = [0.8500, 0.3250, 0.0980]; % Orange
color_light_orange = [0.9290, 0.6940, 0.1250];
color_grey = [0.5, 0.5, 0.5]; % Grey

% Create arrays to store the indices that are already plotted
plotted_indices = false(size(double(eye_tracking_data.systemTimeStamp) / 1e6));

% Loop through each trial in cutData
for trial = 1:size(cutData, 1)
    % Extract the current trial data
    current_data = trial_results(trial, :);
    
    % Extract the relevant times
    trialStartTime = current_data.trialStartTime;
    fixStartTime = current_data.fixStartTime;
    blankStartTime = current_data.blankStartTime;
    StimulusOnsetTime = current_data.StimulusOnsetTime;
    trial_resp_time = current_data.trialEndTime;
    
    % Extract system timestamps and pupil diameters for both eyes
    timestamps = double(eye_tracking_data.systemTimeStamp) / 1e6;
    
    % Plot continuous eyetracking data in different colors for each time segment
    plotted_indices = plot_segment(timestamps, eye_tracking_data.right.pupil.diameter, trialStartTime, fixStartTime, color_blue, plotted_indices);
    plotted_indices = plot_segment(timestamps, eye_tracking_data.right.pupil.diameter, fixStartTime, blankStartTime, color_green, plotted_indices);
    plotted_indices = plot_segment(timestamps, eye_tracking_data.right.pupil.diameter, blankStartTime, StimulusOnsetTime, color_orange, plotted_indices);
    plotted_indices = plot_segment(timestamps, eye_tracking_data.right.pupil.diameter, StimulusOnsetTime, trial_resp_time, color_light_orange, plotted_indices);
    
    plotted_indices = plot_segment(timestamps, eye_tracking_data.left.pupil.diameter, trialStartTime, fixStartTime, color_blue, plotted_indices);
    plotted_indices = plot_segment(timestamps, eye_tracking_data.left.pupil.diameter, fixStartTime, blankStartTime, color_green, plotted_indices);
    plotted_indices = plot_segment(timestamps, eye_tracking_data.left.pupil.diameter, blankStartTime, StimulusOnsetTime, color_orange, plotted_indices);
    plotted_indices = plot_segment(timestamps, eye_tracking_data.left.pupil.diameter, StimulusOnsetTime, trial_resp_time, color_light_orange, plotted_indices);
end

% Plot the data points not already plotted in grey
plot_gaps(timestamps, eye_tracking_data.right.pupil.diameter, ~plotted_indices, color_grey);
plot_gaps(timestamps, eye_tracking_data.left.pupil.diameter, ~plotted_indices, color_grey);

% Create a handmade legend
ylims = get(gca, 'ylim');
xlims = get(gca, 'xlim');
legend_y = ylims(2) - (ylims(2) - ylims(1)) * 0.05;
legend_x_start = xlims(1) + (xlims(2) - xlims(1)) * 0.01;
legend_x_end = legend_x_start + (xlims(2) - xlims(1)) * 0.05;

plot([legend_x_start, legend_x_end], [legend_y, legend_y], 'Color', color_blue, 'LineWidth', 2);
text(legend_x_end + (xlims(2) - xlims(1)) * 0.01, legend_y, 'Press Any Button', 'VerticalAlignment', 'middle');

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

set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1])
saveas(gcf,strcat(analysis_folder, '\pupil_diam_over_time.png'));



%% Function to plot segments with gaps filled with grey
    function plotted_indices = plot_segment(timestamps, data, start_time, end_time, color, plotted_indices)
        idx_segment = timestamps >= start_time & timestamps < end_time;
        idx_gaps = timestamps >= start_time & timestamps < end_time & isnan(data);
        
        % Plot the data segment
        plot(timestamps(idx_segment), data(idx_segment), 'Color', color, 'LineWidth', 1.5);
        plotted_indices = plotted_indices | idx_segment;
        
        % Plot gaps in grey
        if any(idx_gaps)
            plot(timestamps(idx_gaps), data(idx_gaps), 'Color', color_grey, 'LineWidth', 1.5);
        end
    end

%% Function to plot gaps in grey
    function plot_gaps(timestamps, data, idx_gaps, color)
        segments = find(idx_gaps);
        if isempty(segments)
            return; % If all are NaNs, skip
        end
        
        for i = 1:length(segments)-1
            if segments(i+1) - segments(i) == 1  % Adjacent points are valid
                line(timestamps(segments(i):segments(i+1)), data(segments(i):segments(i+1)), 'Color', color, 'LineWidth', 1.5);
            else  % There's a gap
                %                 xi = [timestamps(segments(i)), timestamps(segments(i+1))];
                %                 yi = [data(segments(i)), data(segments(i+1))];
                %                 plot(xi, yi, '--', 'Color', color);  % Grey dashed line
            end
        end
    end

hold off;
end
