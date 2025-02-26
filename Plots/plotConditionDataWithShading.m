function [y_upper,y_lower] = plotConditionDataWithShading(data, x_values, color, condition)
    % Plot mean ± SEM with shading for each condition
    avg = mean(data, 1, 'omitnan');
    std_error = std(data, 0, 1, 'omitnan') / sqrt(size(data, 1));
    nan_mask = ~isnan(avg) & ~isnan(std_error);
    x_shaded = x_values(nan_mask);
    y_upper = avg(nan_mask) + std_error(nan_mask);
    y_lower = avg(nan_mask) - std_error(nan_mask);
    fill([x_shaded, fliplr(x_shaded)], [y_upper, fliplr(y_lower)], color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(x_values, avg, 'Color', color, 'LineWidth', 2, 'DisplayName', condition);
end