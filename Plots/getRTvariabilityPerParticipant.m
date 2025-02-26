function rt_variability = getRTvariabilityPerParticipant(data, conditions, color_map, comparison_results_folder)
    % This function computes the variability (std dev or variance) of reaction times (RTs)
    % for each condition for a single participant, and plots the results.

    % Initialize a struct to store the results
    rt_variability = struct();

    % Loop through each condition
    for c = 1:length(conditions)
        condition = conditions{c};

        % Find the trials that match the current condition
        condition_trials = strcmp(data.Condition, condition);
        rt_button_all = data.rt(condition_trials);
        rt_button_var = std(rt_button_all, 'omitnan');

        % Store the RT variability for the condition
        rt_variability.(condition).RT_ButtonPress_var = rt_button_var;
    end

    % Plot RT variability per condition with error bars
    % Initialize data for plotting
    button_rtv = zeros(1, length(conditions));

    % Extract the data for each condition
    for c = 1:length(conditions)
        condition = conditions{c};
        button_rtv(c) = rt_variability.(condition).RT_ButtonPress_var;
    end

    % Create the plot
    figure;
    hold on;

    % Plot each condition separately with its color from color_map
    for c = 1:length(conditions)
        condition = conditions{c};
        plot(c, button_rtv(c), 'o-', 'Color', color_map(condition), ...
            'MarkerFaceColor', color_map(condition), 'LineWidth', 2, 'MarkerSize', 8);
    end

    % Set plot properties
    set(gca, 'XTick', 1:length(conditions), 'XTickLabel', conditions);
    xlabel('Condition');
    ylabel('RT ButtonPress Variability (Std Dev)');
    title('Reaction Time Variability (ButtonPress) by Condition');
    grid on;

    hold off;

    % Save the figure
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);    
    saveas(gcf, fullfile(comparison_results_folder, 'RTV_single_participant.png'));    

    % Display the RT variability data
    disp('Reaction Time Variability (Standard Deviation) by Condition:');
    disp(rt_variability);
end
