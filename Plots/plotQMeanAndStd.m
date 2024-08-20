function plotQMeanAndStd(adhd_numeric, non_adhd_numeric, color_map, safe, comparison_results_folder)
    % Function to plot the mean and standard deviation of responses by question
    %
    % Inputs:
    %   adhd_numeric: Numeric responses for the ADHD group
    %   non_adhd_numeric: Numeric responses for the non-ADHD group
    %   color_map: Color mapping for groups

    num_questions = size(adhd_numeric, 2);

    % Calculate the mean and standard deviation for each question
    adhd_means = mean(adhd_numeric, 1, 'omitnan');
    adhd_stds = std(adhd_numeric, 0, 1, 'omitnan');
    non_adhd_means = mean(non_adhd_numeric, 1, 'omitnan');
    non_adhd_stds = std(non_adhd_numeric, 0, 1, 'omitnan');

    % Create the plot
    figure;
    hold on;

    % Plot the mean and standard deviation for non-ADHD group
    errorbar(1:num_questions, non_adhd_means, non_adhd_stds, '-o', 'Color', color_map('nonADHD'), 'LineWidth', 1.5, ...
        'MarkerFaceColor', color_map('nonADHD'), 'MarkerSize', 8);

    % Plot the mean and standard deviation for ADHD group
    errorbar(1:num_questions, adhd_means, adhd_stds, '-o', 'Color', color_map('ADHD'), 'LineWidth', 1.5, ...
        'MarkerFaceColor', color_map('ADHD'), 'MarkerSize', 8);

    xlabel('Question Number');
    ylabel('Response');
    title('Mean and Standard Deviation of Responses by Question');
    xticks(1:num_questions);
    xlim([1 num_questions]);
    yticks(1:5);
    yticklabels({'Never', 'Rarely', 'Sometimes', 'Often', 'Very Often'});
    legend({'Non-ADHD', 'ADHD'}, 'Location', 'eastoutside');
    grid on;
    hold off;
    
    if safe == 1
        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf, fullfile(comparison_results_folder, 'quest_line_all_questions_group_mean_std.png'));
end
end
