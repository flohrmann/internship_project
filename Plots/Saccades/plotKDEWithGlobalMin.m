function plotKDEWithGlobalMin(dni_f, dni_xi, threshold_dist_non_interp, dni_global_min_value, ...
                              vni_f, vni_xi, threshold_vel_non_interp, vni_global_min_value, ...
                              di_f, di_xi, threshold_dist_interp, di_global_min_value, ...
                              vi_f, vi_xi, threshold_vel_interp, vi_global_min_value, ...
                              safe, comparison_results_folder, trial)
    % Function to plot KDEs with global minimum thresholds in a 2x2 subplot layout
    % Inputs:
    %   dni_f, dni_xi - KDE results for non-interpolated distance
    %   threshold_dist_non_interp - Global minimum threshold for non-interpolated distance
    %   dni_global_min_value - KDE value at global minimum for non-interpolated distance
    %   vni_f, vni_xi - KDE results for non-interpolated velocity
    %   threshold_vel_non_interp - Global minimum threshold for non-interpolated velocity
    %   vni_global_min_value - KDE value at global minimum for non-interpolated velocity
    %   di_f, di_xi - KDE results for interpolated distance
    %   threshold_dist_interp - Global minimum threshold for interpolated distance
    %   di_global_min_value - KDE value at global minimum for interpolated distance
    %   vi_f, vi_xi - KDE results for interpolated velocity
    %   threshold_vel_interp - Global minimum threshold for interpolated velocity
    %   vi_global_min_value - KDE value at global minimum for interpolated velocity
    %   safe - if 1 safe png if 0 dont
    %   comparison_results_folder - folder to save the figure

    figure;

    % Subplot 1: Non-Interpolated Data - Distance KDE
    subplot(2, 2, 1);
    plot(dni_xi, dni_f, 'LineWidth', 2);
    hold on;
    plot(threshold_dist_non_interp, dni_global_min_value, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');  % Plot the threshold
    xlabel('Distance (pixels)');
    ylabel('Density');
    title('Non-Interpolated Data: Distance KDE');
    grid on;
    legend('KDE', 'Global Minimum');

    % Subplot 2: Non-Interpolated Data - Velocity KDE
    subplot(2, 2, 2);
    plot(vni_xi, vni_f, 'LineWidth', 2);
    hold on;
    plot(threshold_vel_non_interp, vni_global_min_value, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');  % Plot the threshold
    xlabel('Velocity (pixels/second)');
    ylabel('Density');
    title('Non-Interpolated Data: Velocity KDE');
    grid on;
    legend('KDE', 'Global Minimum');

    % Subplot 3: Interpolated Data - Distance KDE
    subplot(2, 2, 3);
    plot(di_xi, di_f, 'LineWidth', 2);
    hold on;
    plot(threshold_dist_interp, di_global_min_value, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');  % Plot the threshold
    xlabel('Distance (pixels)');
    ylabel('Density');
    title('Interpolated Data: Distance KDE');
    grid on;
    legend('KDE', 'Global Minimum');

    % Subplot 4: Interpolated Data - Velocity KDE
    subplot(2, 2, 4);
    plot(vi_xi, vi_f, 'LineWidth', 2);
    hold on;
    plot(threshold_vel_interp, vi_global_min_value, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');  % Plot the threshold
    xlabel('Velocity (pixels/second)');
    ylabel('Density');
    title('Interpolated Data: Velocity KDE');
    grid on;
    legend('KDE', 'Global Minimum');

    % Finalize and Save Plot
    sgtitle(strcat('Trial ', num2str(trial),': Kernel Density Estimates with Global Minimum Thresholds'));
    
    if safe == 1
        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf, fullfile(comparison_results_folder, strcat('kde_global_minimum_plot_', num2str(trial), '.png')));
    else
    end
end
