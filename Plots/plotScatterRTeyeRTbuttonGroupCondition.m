function plotScatterRTeyeRTbuttonGroupCondition(groups, conditions, rt_eye_conditions, rt_button_press_conditions, color_map, safe, comparison_results_folder)

figure;
for g = 1:length(groups)
    for c = 1:length(conditions)
        subplot(length(groups), length(conditions), (g-1)*length(conditions) + c);
        
        % Get the color for the current condition
        condition_color = color_map(conditions{c});
        
        % Scatter plot with the respective color outline
        scatter(rt_eye_conditions{g, c}, rt_button_press_conditions{g, c}, 36, 'MarkerEdgeColor', condition_color, 'MarkerFaceColor', 'none');
        
        xlabel('Eye RT');
        ylabel('Button Press RT');
        title([groups{g}, ' - ', conditions{c}]);
        grid on;
    end
end
if safe == 1
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf, fullfile(comparison_results_folder, 'scatter_RTeye_RTbutton_group_condition.png'));
end