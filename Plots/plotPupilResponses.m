function plotPupilResponses(processed_data, id, conditions, analysis_folder, color_map, sr, base)
% Define the different methods to process: 'Mean' and 'Median'
methods = {'Mean', 'Median'};

% Loop through both 'Mean' and 'Median' to plot both cases
for m = 1:length(methods)
    method = methods{m};  % Get the current method (Mean/Median)
    
    figure;
    hold on;
    
    % Define the number of time points before and after based on the first condition
    num_before = size(processed_data.a.(['before' method]), 2);
    num_after = size(processed_data.a.(['after' method]), 2);
    
    % Define the custom x-axis for the time points
    x_before = linspace(-sr * (num_before - 1), 0, num_before);
    x_after = linspace(sr, sr * num_after, num_after);
    x_values = [x_before, x_after];  % Time axis with 0 being the target found
    
    % Loop through conditions and plot the processed data
    for i = 1:length(conditions)
        condition = conditions{i};
        
        % Extract the processed data for this condition (mean or median)
        data_before = processed_data.(condition).(['before' method]);
        data_after = processed_data.(condition).(['after' method]);
        
        % Combine before and after data
        combined_data = [data_before, data_after];
        
        % Compute the mean and standard error
        if strcmp(method, 'Mean')
            avg = mean(combined_data, 1, 'omitnan');
        else % median
            avg = median(combined_data, 1, 'omitnan');
        end
        std_error = std(combined_data, 0, 1, 'omitnan') / sqrt(height(combined_data));
        
        % Remove any NaNs from the shading
        nan_mask = ~isnan(avg) & ~isnan(std_error);
        x_shaded = x_values(nan_mask);
        y_upper = avg(nan_mask) + std_error(nan_mask);
        y_lower = avg(nan_mask) - std_error(nan_mask);
        
        % Get the color for the current condition
        color = color_map(condition);
        
        % Plot the shaded region (mean ± standard error)
        fill([x_shaded, fliplr(x_shaded)], [y_upper, fliplr(y_lower)], ...
            color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        % Plot the mean/median line
        plot(x_values, avg, 'Color', color, 'LineWidth', 2);
    end
    xline(0, '--', 'LineWidth', 1.5, 'Color', [0.5, 0.5, 0.5]);
    plot([0 0], ylim, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5, 'Visible', 'off');  % Hidden plot for legend
    xlabel('Time (s)');
    xlim([-sr * (num_before - 1), sr * num_after]);
    ylabel(['Pupil Dilation (z-scored) [', method, ']']);
    title(['Average Pupil Change during a Trial - Participant ', num2str(id)]);
    
    legend_objects = findobj(gca, 'Type', 'Line');  % Find all line objects in the plot
    legend(legend_objects(end:-1:1), 'a', 'a simple', 'b', 'b simple', 'Gaze Reached Target', ...
        'Location', 'northwest', 'Orientation', 'vertical');
    legend('boxoff');
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf, fullfile(analysis_folder, ['norm_', base, '_', method, '_processed_pupil_diam_participant_', num2str(id), '.png']));
    hold off;
end
end
