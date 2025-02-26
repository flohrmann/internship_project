function plotGroupConditionCorrelations(R_group_condition, groups, conditions)

    % Loop through each group
    for g = 1:length(groups)
        group = groups{g};
        
        % Create a new figure for this group
        figure('Name', ['Correlation Heatmaps - Group: ', group], 'NumberTitle', 'off');
        sgtitle(['Correlation Heatmaps for ', group, ' Group']);
        
        % Loop through each condition and plot the heatmap
        for c = 1:length(conditions)
            condition = conditions{c};
            
            % Get the correlation matrix for this group and condition
            R_matrix = R_group_condition.(group){c};
            
            % Create a subplot for each condition
            subplot(2, 2, c);
            heatmap(R_matrix, 'Colormap', parula, 'ColorbarVisible', 'on');
            
            % Set axis labels and title
            title(['Condition: ', condition]);
            xlabel('Variables (RT_ButtonPress, RT_Eye, Accuracy)');
            ylabel('Variables (RT_ButtonPress, RT_Eye, Accuracy)');
            ax = gca;
            ax.XDisplayLabels = {'RT_ButtonPress', 'RT_Eye', 'Accuracy'};
            ax.YDisplayLabels = {'RT_ButtonPress', 'RT_Eye', 'Accuracy'};
        end
        
        % Adjust the layout
        set(gcf, 'Position', [100, 100, 1000, 800]);
    end
end
