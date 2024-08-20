function plotSaccades(x, y, saccade_onsets, saccade_offsets, title_text)
    hold on;
    plot(x, y, 'k.-', 'MarkerSize', 10, 'DisplayName', 'Gaze Path');
    colors = lines(length(saccade_onsets));  % Generate different colors
    for i = 1:length(saccade_onsets)
        plot(x(saccade_onsets(i):saccade_offsets(i)), y(saccade_onsets(i):saccade_offsets(i)), ...
             'Color', colors(i,:), 'LineWidth', 2);
    end
    xlabel('Horizontal Gaze Position (pixels)');
    ylabel('Vertical Gaze Position (pixels)');
    title(title_text);
    grid on;
    hold off;
end
