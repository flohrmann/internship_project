
function plotRTscatterConditionAcrossGroups(groups, conditions, rt_eye_conditions, rt_button_press_conditions, color_map, safe, comparison_results_folder)
    % Calculate uniform axis limits across all conditions and groups
    all_eye_rt = cell2mat(rt_eye_conditions(:));
    all_button_rt = cell2mat(rt_button_press_conditions(:));
    x_limits = [min(all_eye_rt), max(all_eye_rt)];
    y_limits = [min(all_button_rt), max(all_button_rt)];

    % Set up the tiled layout with one extra tile for the legend
    figure;
    t = tiledlayout(1, length(conditions) + 1, 'TileSpacing', 'compact', 'Padding', 'compact');  % Layout for conditions + legend

    % Initialize arrays for legend handles and labels
    legend_handles = []; 
    legend_labels = {};

    for c = 1:length(conditions)
        nexttile;
        hold on;
        
        for g = 1:length(groups)
            % Remove NaNs and align data
            eye_rt_data = rt_eye_conditions{g, c};
            button_rt_data = rt_button_press_conditions{g, c};
            valid_idx = ~isnan(eye_rt_data) & ~isnan(button_rt_data);
            eye_rt_clean = eye_rt_data(valid_idx);
            button_rt_clean = button_rt_data(valid_idx);

            % Transform data to log scale
            eye_rt_log = log10(eye_rt_clean);
            button_rt_log = log10(button_rt_clean);

            % Plot scatter for each group with different colors
            scatter_handle = scatter(eye_rt_log, button_rt_log, 36, ...
                    'MarkerEdgeColor', color_map(groups{g}), 'DisplayName', groups{g}, 'MarkerFaceColor', 'none');
            
            % Add trend line for each group in log scale if enough points
            if length(eye_rt_log) > 10
                coeffs = polyfit(eye_rt_log, button_rt_log, 1);
                
                % Define x_fit based on the range of eye_rt_log
                x_fit = linspace(min(eye_rt_log), max(eye_rt_log), 100);
                
                y_fit = polyval(coeffs, x_fit);
                line_handle = plot(x_fit, y_fit, 'Color', color_map(groups{g}), 'LineWidth', 1.5, 'DisplayName', [groups{g}, ' Trend']);
            end
            
            % Add scatter and line handles to legend arrays if not already added
            if c == 1  % Add legend for only the first subplot
                legend_handles = [legend_handles, scatter_handle, line_handle];
                legend_labels = [legend_labels, groups{g}, [groups{g}, ' Trend']];
            end
        end

        % Set axis limits, log scale, and add labels
        xlim(log10(x_limits));
        ylim(log10(y_limits));
        xlabel('Eye RT (log scale)');
        ylabel('Button Press RT (log scale)');
        title(['Condition: ', conditions{c}]);
        grid off;
        hold off;
    end

    % Create a dedicated tile for the legend outside the main plots
    nexttile(length(conditions) + 1);
    axis off;
    legend(legend_handles, legend_labels);

    % Save the second figure if specified
    if safe == 1
        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf, fullfile(comparison_results_folder, 'scatter_RTeye_RTbutton_per_condition.png'));
    end
end
