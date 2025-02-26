function plotAccuracyPerCondition(trial_results, color_map, safe, analysis_folder)
    % Initialize arrays
    accuracy_all = [];
    condition_all = {};

    % Populate arrays with accuracy and condition information
    for i = 1:size(trial_results,1)
        accuracy_all = [accuracy_all; trial_results.correct(i)];
        condition_all = [condition_all; trial_results.Condition(i)];
    end

    % Convert condition labels to categorical for easier comparison
    condition_all = categorical(condition_all);

    % Calculate mean accuracy and SEM per condition
    [unique_conditions, ~, cond_idx] = unique(condition_all);
    mean_accuracy = zeros(length(unique_conditions), 1);
    sem_accuracy = zeros(length(unique_conditions), 1);

    for i = 1:length(unique_conditions)
        % Extract accuracy for each condition
        condition_accuracy = accuracy_all(cond_idx == i);

        % Calculate mean accuracy and SEM for each condition
        mean_accuracy(i) = mean(condition_accuracy);
        sem_accuracy(i) = std(condition_accuracy) / sqrt(length(condition_accuracy)); % Standard Error of Mean
    end

    % Create figure
    figure;
    hold on;

    % Plot mean points with error bars for each condition
    x_data = 1:length(unique_conditions);
    y_data = mean_accuracy;
    error_data = sem_accuracy;
    for i = 1:length(unique_conditions)
        condition = char(unique_conditions(i));
        errorbar(x_data(i), y_data(i), error_data(i), 'o', 'Color', color_map(condition), ...
                 'MarkerFaceColor', color_map(condition), 'LineWidth', 2, 'MarkerSize', 8);
    end

    % Customize the plot
    xticks(1:length(unique_conditions));
    xticklabels(unique_conditions);
    xlabel('Condition');
    ylabel('Mean Accuracy');
    title('Mean Accuracy per Condition');
    legend(unique_conditions, 'Location', 'northeastoutside');
    hold off;

    % Save the figure
    if safe == 1
        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf, fullfile(analysis_folder, 'mean_accuracy_per_condition.png'));
    end
end
