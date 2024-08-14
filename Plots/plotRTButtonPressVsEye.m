function plotRTButtonPressVsEye(data, comparison_results_folder, safe)
    % Initialize variables for storing reaction times and eye movement times
    rt_button_press = [];
    rt_eye = [];

    % Loop through the struct to gather data for RT comparison
    for i = 1:length(data)
        rt_button_press = [rt_button_press; data(i).rt];
        rt_eye = [rt_eye; data(i).rt_eye];
    end

    % Create figure
    figure;

    % Scatter plot of RT button press vs RT eye movement
    scatter(rt_button_press, rt_eye, 50, 'b', 'filled');
    hold on;

    % Plot the line of equality (y = x)
    max_val = max(max(rt_button_press), max(rt_eye));
    plot([0, max_val], [0, max_val], 'k--', 'LineWidth', 1.5);

    % Customize the plot
    xlabel('Button Press Reaction Time (s)');
    ylabel('Eye Movement Reaction Time (s)');
    title('Comparison of Reaction Times: Button Press vs Eye Movement');
    set(gca, 'XScale', 'log', 'YScale', 'log');
    grid on;
    axis equal;

    % Save the figure if safe is set to 1
    if safe == 1
        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf, fullfile(comparison_results_folder, 'rt_comparison_button_press_vs_eye_movement.png'));
    end
end
