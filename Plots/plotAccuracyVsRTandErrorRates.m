function plotAccuracyVsRTandErrorRates(data, color_map)
    % Plot 1: Scatter plot of accuracy vs. reaction time
    figure;
    hold on;
    % Iterate through each participant in the data
    for i = 1:length(data)
        scatter(data(i).rt, data(i).accuracy, 'DisplayName', data(i).group, ...
                'MarkerFaceColor', color_map(data(i).group), 'MarkerEdgeColor', color_map(data(i).group));
    end
    xlabel('Reaction Time (RT)');
    ylabel('Accuracy');
    title('Accuracy vs. Reaction Time');
    legend('Location', 'northeastoutside');
    grid on;
    hold off;

    % Plot 2: Line plot of wrong trials for simple vs. nonsimple conditions per group
    figure;
    hold on;
    % Initialize arrays to collect wrong trial counts for plotting
    groups = unique({data.group});  % Get unique group labels

    % Define conditions for simple and nonsimple trials
    simple_conditions = {'a_simple', 'b_simple'};
    nonsimple_conditions = {'a', 'b'};

    % Iterate through each group
    for g = 1:length(groups)
        group_name = groups{g};

        % Get data only for participants in the current group
        group_data = data(strcmp({data.group}, group_name));
        
        % Initialize arrays to store counts for each participant in the group
        simple_wrong_counts = zeros(1, length(group_data));
        nonsimple_wrong_counts = zeros(1, length(group_data));

        % Calculate wrong trials for simple and nonsimple conditions per participant
        for i = 1:length(group_data)
            simple_wrong_counts(i) = sum(~group_data(i).accuracy & ismember(group_data(i).condition, simple_conditions));
            nonsimple_wrong_counts(i) = sum(~group_data(i).accuracy & ismember(group_data(i).condition, nonsimple_conditions));
        end

        % Plot line plot for the current group with simple and nonsimple points connected
        plot(ones(1, length(simple_wrong_counts)), simple_wrong_counts, 'o-', 'MarkerSize', 8, ...
            'Color', color_map(group_name), 'DisplayName', [group_name, ' Simple']);
        plot(2 * ones(1, length(nonsimple_wrong_counts)), nonsimple_wrong_counts, 'o-', 'MarkerSize', 8, ...
            'Color', color_map(group_name), 'DisplayName', [group_name, ' Nonsimple']);
        
        % Connect points for each participant in the group
        for i = 1:length(simple_wrong_counts)
            plot([1, 2], [simple_wrong_counts(i), nonsimple_wrong_counts(i)], '-', ...
                 'Color', color_map(group_name), 'LineWidth', 1);
        end
    end

    % Set x-axis labels and legend
    xticks([1, 2]);
    xticklabels({'Simple', 'Nonsimple'});
    ylabel('Wrong Trials Count');
    title('Wrong Trials: Simple vs. Nonsimple Conditions by Group');
    legend('Location', 'northeastoutside');
    grid on;
    hold off;
end
