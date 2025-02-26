function plotRTScatterGroupCondition(groups, conditions, rt_eye_conditions, rt_button_press_conditions, color_map, safe, comparison_results_folder)
    % Calculate uniform axis limits across all conditions and groups
    all_eye_rt = cell2mat(rt_eye_conditions(:));
    all_button_rt = cell2mat(rt_button_press_conditions(:));
    x_limits = [min(all_eye_rt), max(all_eye_rt)];
    y_limits = [min(all_button_rt), max(all_button_rt)];

    figure;
    for g = 1:length(groups)
        for c = 1:length(conditions)
            % Calculate subplot index
            subplot(length(groups), length(conditions), (g-1)*length(conditions) + c);
            condition_color = color_map(conditions{c});

            % Remove NaNs and align data
            eye_rt_data = rt_eye_conditions{g, c};
            button_rt_data = rt_button_press_conditions{g, c};
            valid_idx = ~isnan(eye_rt_data) & ~isnan(button_rt_data);
            eye_rt_clean = eye_rt_data(valid_idx);
            button_rt_clean = button_rt_data(valid_idx);

            % Transform data to log scale
            eye_rt_log = log10(eye_rt_clean);
            button_rt_log = log10(button_rt_clean);
            
            % Plot scatter for current group and condition
            scatter(eye_rt_log, button_rt_log, 36, ...
                    'MarkerEdgeColor', condition_color, 'MarkerFaceColor', 'none');
            hold on;
            
            % Add trend line only if enough data points are available
            if length(eye_rt_log) > 10
                coeffs = polyfit(eye_rt_log, button_rt_log, 1);
                x_fit = linspace(min(eye_rt_log), max(eye_rt_log), 100);
                y_fit = polyval(coeffs, x_fit);
                plot(x_fit, y_fit, 'Color', condition_color, 'LineWidth', 1.5);
            end

            % Set axis limits and log scale for plotting
            xlim(log10(x_limits));
            ylim(log10(y_limits));
            xlabel('Eye RT (log scale)');
            ylabel('Button Press RT (log scale)');
            title([groups{g}, ' - ', conditions{c}]);
            grid off;
        end
    end

    % Save first figure if specified
    if safe == 1
        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf, fullfile(comparison_results_folder, 'scatter_RTeye_RTbutton_group_condition.png'));
    end
end