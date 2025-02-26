function plotViolinNumFixations(fixationStats, group_labels, conditions, color_map, comparison_results_folder, safe)
    % Initialize containers for ADHD and non-ADHD groups
    mean_adhd_conditions = [];
    mean_nonadhd_conditions = [];
    
    % Loop through participants to calculate mean number of fixations per condition
    for i = 1:length(fixationStats)
        trials = fixationStats(i).trials;
        numFixations = [trials.numFixations]; % Extract number of fixations
        trialConditions = {trials.conditions}; % Extract conditions for each trial
        
        % Aggregate data by condition
        participant_means = zeros(1, length(conditions));
        for c = 1:length(conditions)
            % Extract the string content from the nested cells
            trialConditionsStrings = cellfun(@(x) x{1}, trialConditions, 'UniformOutput', false);
            condition_idx = strcmp(trialConditionsStrings, conditions{c});
            participant_means(c) = mean(numFixations(condition_idx)); % Mean for this condition
        end
        
        % Assign to ADHD or non-ADHD group
        if strcmp(group_labels{i}, 'ADHD')
            mean_adhd_conditions = [mean_adhd_conditions; participant_means];
        elseif strcmp(group_labels{i}, 'nonADHD')
            mean_nonadhd_conditions = [mean_nonadhd_conditions; participant_means];
        end
    end

    % Prepare data for violin plots
    data_for_violin_mean_adhd = cell(1, length(conditions));
    data_for_violin_mean_nonadhd = cell(1, length(conditions));
    facecolors_adhd = zeros(length(conditions), 3);
    facecolors_nonadhd = zeros(length(conditions), 3);

    for i = 1:length(conditions)
        % Collect per-condition data
        data_for_violin_mean_adhd{i} = mean_adhd_conditions(:, i);
        data_for_violin_mean_nonadhd{i} = mean_nonadhd_conditions(:, i);
        
        % Set face colors
        facecolors_adhd(i, :) = color_map(conditions{i});
        facecolors_nonadhd(i, :) = color_map(conditions{i});
    end

    % Create figure with two subplots for number of fixations
    figure;

    % Subplot for number of fixations - non-ADHD group
    ax1 = subplot(1, 2, 1);
    violin(data_for_violin_mean_nonadhd, 'xlabel', conditions, 'facecolor', facecolors_nonadhd, ...
           'edgecolor', 'none', 'facealpha', 0.5, 'mc', 'k', 'medc', 'k--');
    ylabel('Mean Number of Fixations');
    xlabel('Condition');
    title('non-ADHD');

    % Subplot for number of fixations - ADHD group
    ax2 = subplot(1, 2, 2);
    violin(data_for_violin_mean_adhd, 'xlabel', conditions, 'facecolor', facecolors_adhd, ...
           'edgecolor', 'none', 'facealpha', 0.5, 'mc', 'k', 'medc', 'k--');
    ylabel('Mean Number of Fixations');
    xlabel('Condition');
    title('ADHD');

    % Link y-axes for consistency
    linkaxes([ax1, ax2], 'y');

    % Set a global title for the plot
    sgtitle('Comparison of Mean Number of Fixations per Condition');

    % Save the plot if required
    if safe == 1
        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf, fullfile(comparison_results_folder, '01_fixation_violin_mean_num_fixations_condition_group.png'));
    end
end
