function plotPupilDiameterResponses4subplotsAveragePerCondition( ...
    data1, data_type_1, title1, t01, ...
    data2, data_type_2, title2, t02,...
    data3, data_type_3, title3, t03,...
    data4, data_type_4, title4, t04,...
    big_title, conditions, condition_labels, time_vector, color_map, ...
    fullscreen, new_x_limits, new_y_limits, dimensions)
    % Initialize y-limits for unified scaling across tiles
    min_y = [];
    max_y = [];

    % Create tiled layout
    tiledlayout(2, 2);  % 2x2 grid for the four subplots
    sgtitle(big_title);

    % Prepare a cell array for data and titles for easier looping
    dataArray = {data1, data2, data3, data4};
    titleArray = {title1, title2, title3, title4};
    typeArray = {data_type_1, data_type_2, data_type_3, data_type_4};
    t0Array = {t01, t02, t03, t04};
    % Loop through each tile for each dataset and plot accordingly
    for idx = 1:4
        % Select next tile
        nexttile; hold on;
        title(titleArray{idx});
        data = dataArray{idx};
        type = typeArray{idx};
        % Plot each condition
        for i = 1:length(conditions)
            condition = conditions{i};
            condition_trials = data(strcmp(data.Condition, condition), :);
            [y_upper, y_lower] = plotConditionDataWithShading(condition_trials.(type), time_vector, color_map(condition), condition);
            min_y = [min_y, y_lower];
            max_y = [max_y, y_upper];
        end

        % Add line at t0 and labels
        xline(0, '--', t0Array{idx}, 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3], 'LabelOrientation', 'horizontal', 'HandleVisibility', 'off');
        xlabel('Time from t0 (s)');
        ylabel('Pupil Diameter (mm)');
        
        try
            xlim([new_x_limits(1), new_x_limits(2)]);
        catch
            xlim([min(time_vector), max(time_vector)]);
            disp('No extra x limits defined')
        end
        
        hold off;
    end

    
    try % if y limits explicitly stated
        for idx = 1:4
            nexttile(idx);
            ylim(new_y_limits);
        end
    catch % Set optimal unified y-axis limits across all tiles
        y_limits = [min(min_y), max(max_y)];
        for idx = 1:4
            nexttile(idx);
            ylim(y_limits);
        end
    end
    % Add a single legend for all tiles
    lgd = legend(condition_labels); 
    set(lgd, 'Orientation', 'horizontal', 'Location', 'southoutside', 'Box', 'off');


%     legend(condition_labels); legend('boxoff');
%     legend('Location', 'bestoutside');
    try
       set(gcf, 'Units', 'normalized', 'OuterPosition', dimensions); % Example: make it wider
    catch
    end
    if fullscreen == 1
        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    else
    end
end
