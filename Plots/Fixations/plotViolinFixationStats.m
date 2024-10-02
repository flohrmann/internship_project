function plotViolinFixationStats(fixationStats, group_labels, conditions, color_map, comparison_results_folder, safe)
    % Initialize containers for ADHD and nonADHD data
    duration_adhd_conditions = [];
    duration_nonadhd_conditions = [];
    
    condition_adhd = {};
    condition_nonadhd = {};

    % Loop through the struct to gather data for ADHD and non-ADHD groups
    for i = 1:length(fixationStats)
        for c = 1:length(conditions)
            conditionIdx = strcmp({fixationStats(i).conditions.name}, conditions{c});
            
            if any(conditionIdx)
                conditionData = fixationStats(i).conditions(conditionIdx);
                
                if strcmp(group_labels{i}, 'ADHD')
                    % Append the fixation durations for each condition
                    duration_adhd_conditions = [duration_adhd_conditions; conditionData.fixationDurations(:)];
                    condition_adhd = [condition_adhd; repmat({conditions{c}}, length(conditionData.fixationDurations), 1)];
                elseif strcmp(group_labels{i}, 'nonADHD')
                    % Append the fixation durations for each condition
                    duration_nonadhd_conditions = [duration_nonadhd_conditions; conditionData.fixationDurations(:)];
                    condition_nonadhd = [condition_nonadhd; repmat({conditions{c}}, length(conditionData.fixationDurations), 1)];
                end
            end
        end
    end

    % Convert the condition labels to categorical and extract unique conditions
    condition_adhd = categorical(condition_adhd);
    condition_nonadhd = categorical(condition_nonadhd);
    unique_conditions = categories(condition_adhd); % Assuming both groups have the same conditions

    % Prepare data for ADHD and non-ADHD violin plots
    data_for_violin_duration_adhd = cell(1, length(unique_conditions));
    data_for_violin_duration_nonadhd = cell(1, length(unique_conditions));
    facecolors_adhd = zeros(length(unique_conditions), 3);
    facecolors_nonadhd = zeros(length(unique_conditions), 3);

    for i = 1:length(unique_conditions)
        % Log-transform the data before storing it in the cell arrays
        data_for_violin_duration_adhd{i} = log10(duration_adhd_conditions(condition_adhd == unique_conditions{i}));
        data_for_violin_duration_nonadhd{i} = log10(duration_nonadhd_conditions(condition_nonadhd == unique_conditions{i}));
        
        % Set face colors
        facecolors_adhd(i, :) = color_map(char(unique_conditions{i}));
        facecolors_nonadhd(i, :) = color_map(char(unique_conditions{i}));
    end

    % Create figure with two subplots for fixation duration
    figure;

    % Subplot for Fixation Duration - non-ADHD group
    subplot(1, 2, 1);
    [h1, L1, MX1, MED1, bw1] = violin(data_for_violin_duration_nonadhd, 'xlabel', unique_conditions, 'facecolor', facecolors_nonadhd, 'edgecolor', 'none', 'facealpha', 0.5, 'mc', 'k', 'medc', 'k--');
    ylabel('log(Fixation Duration (ms))');
    xlabel('Condition');
    title('non-ADHD');

    % Capture the y-axis limits
    ylim_nonadhd_duration = ylim();

    % Subplot for Fixation Duration - ADHD group
    subplot(1, 2, 2);
    [h2, L2, MX2, MED2, bw2] = violin(data_for_violin_duration_adhd, 'xlabel', unique_conditions, 'facecolor', facecolors_adhd, 'edgecolor', 'none', 'facealpha', 0.5, 'mc', 'k', 'medc', 'k--');
    ylabel('log(Fixation Duration (ms))');
    xlabel('Condition');
    title('ADHD');

    % Capture the y-axis limits
    ylim_adhd_duration = ylim();

    % Determine the global y-axis limits
    global_ylim_duration = [min(ylim_nonadhd_duration(1), ylim_adhd_duration(1)), max(ylim_nonadhd_duration(2), ylim_adhd_duration(2))];

    % Apply the global y-axis limits to both subplots
    subplot(1, 2, 1);
    ylim(global_ylim_duration);

    subplot(1, 2, 2);
    ylim(global_ylim_duration);

    sgtitle('Comparison of Log-Transformed Fixation Durations per Condition');

    if safe == 1
        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf, fullfile(comparison_results_folder, 'violin_log_fixation_duration_condition_group.png'));
    end
end
