function plotStimFound30Points(avg_stim_30, id, analysis_folder, color_map)
    figure;
    hold on;
    
    % new x axis 
    num_points = 30;
    x_values = linspace(-0.016 * (num_points - 1), 0, num_points); 

    for i = 1:height(avg_stim_30)
        condition = avg_stim_30.Condition{i};  
        data_points = avg_stim_30.DataPoints{i};  
        color = color_map(condition);
        plot(x_values, data_points, 'Color', color, 'LineWidth', 1.5);
    end

    xlabel('Time before Button Press (s))');
    xlim([-0.016 * (num_points - 1), 0]);
    ylabel('Pupil Diameter (mm)');
    title(['Pupil Diameters Before Button Press - Participant ', num2str(id)]);

    legend('a', 'a simple', 'b', 'b simple', 'Location', 'northwest');
    legend('boxoff');
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1])
    saveas(gcf, fullfile(analysis_folder, 'button_press_30_points_by_condition.png'));
    hold off;
end
