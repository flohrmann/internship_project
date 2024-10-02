function plotSaccadesViolin(saccadeStats, group_labels, conditions, color_map, comparison_results_folder, safe)
    % Initialize containers for ADHD and nonADHD data
    amplitude_adhd_conditions = [];
    amplitude_nonadhd_conditions = [];
    velocity_adhd_conditions = [];
    velocity_nonadhd_conditions = [];
    
    condition_adhd = {};
    condition_nonadhd = {};

    % Loop through the struct to gather data for ADHD and non-ADHD groups
    for i = 1:length(saccadeStats)
        for c = 1:length(conditions)
            conditionIdx = strcmp({saccadeStats(i).conditions.name}, conditions{c});
            
            if any(conditionIdx)
                conditionData = saccadeStats(i).conditions(conditionIdx);
                
                if strcmp(group_labels{i}, 'ADHD')
                    % Append the saccade amplitudes and velocities for each condition
                    amplitude_adhd_conditions = [amplitude_adhd_conditions; conditionData.saccadeAmplitudes(:)];
                    velocity_adhd_conditions = [velocity_adhd_conditions; conditionData.saccadeVelocities(:)];
                    condition_adhd = [condition_adhd; repmat({conditions{c}}, length(conditionData.saccadeAmplitudes), 1)];
                elseif strcmp(group_labels{i}, 'nonADHD')
                    % Append the saccade amplitudes and velocities for each condition
                    amplitude_nonadhd_conditions = [amplitude_nonadhd_conditions; conditionData.saccadeAmplitudes(:)];
                    velocity_nonadhd_conditions = [velocity_nonadhd_conditions; conditionData.saccadeVelocities(:)];
                    condition_nonadhd = [condition_nonadhd; repmat({conditions{c}}, length(conditionData.saccadeAmplitudes), 1)];
                end
            end
        end
    end

    % Convert the condition labels to categorical and extract unique conditions
    condition_adhd = categorical(condition_adhd);
    condition_nonadhd = categorical(condition_nonadhd);
    unique_conditions = categories(condition_adhd); % Assuming both groups have the same conditions

    % Prepare data for ADHD and non-ADHD violin plots
    data_for_violin_amplitude_adhd = cell(1, length(unique_conditions));
    data_for_violin_amplitude_nonadhd = cell(1, length(unique_conditions));
    data_for_violin_velocity_adhd = cell(1, length(unique_conditions));
    data_for_violin_velocity_nonadhd = cell(1, length(unique_conditions));
    facecolors_adhd = zeros(length(unique_conditions), 3);
    facecolors_nonadhd = zeros(length(unique_conditions), 3);


    for i = 1:length(unique_conditions)
        % Log-transform the data before storing it in the cell arrays
        data_for_violin_amplitude_adhd{i} = log10(amplitude_adhd_conditions(condition_adhd == unique_conditions{i}));
        data_for_violin_amplitude_nonadhd{i} = log10(amplitude_nonadhd_conditions(condition_nonadhd == unique_conditions{i}));
        data_for_violin_velocity_adhd{i} = log10(velocity_adhd_conditions(condition_adhd == unique_conditions{i}));
        data_for_violin_velocity_nonadhd{i} = log10(velocity_nonadhd_conditions(condition_nonadhd == unique_conditions{i}));
        
        % Set face colors
        facecolors_adhd(i, :) = color_map(char(unique_conditions{i}));
        facecolors_nonadhd(i, :) = color_map(char(unique_conditions{i}));
    end

    %% Create figure with two subplots for amplitude
    figure;
    % Subplot for Saccade Amplitude - non-ADHD group
    subplot(1, 2, 1);
    [h1, L1, MX1, MED1, bw1] = violin(data_for_violin_amplitude_nonadhd, 'xlabel', unique_conditions, 'facecolor', facecolors_nonadhd, 'edgecolor', 'none', 'facealpha', 0.5, 'mc', 'k', 'medc', 'k--');
    ylabel('log Saccade Amplitude (pixels)');
    xlabel('Condition');
    title('non-ADHD');
    % Capture the y-axis limits
    ylim_nonadhd_amplitude = ylim();

    % Subplot for Saccade Amplitude - ADHD group
    subplot(1, 2, 2);
    [h2, L2, MX2, MED2, bw2] = violin(data_for_violin_amplitude_adhd, 'xlabel', unique_conditions, 'facecolor', facecolors_adhd, 'edgecolor', 'none', 'facealpha', 0.5, 'mc', 'k', 'medc', 'k--');
    ylabel('log Saccade Amplitude (pixels)');
    xlabel('Condition');
    title('ADHD');
    % Capture the y-axis limits
    ylim_adhd_amplitude = ylim();

    % Determine the global y-axis limits
    global_ylim_amplitude = [min(ylim_nonadhd_amplitude(1), ylim_adhd_amplitude(1)), max(ylim_nonadhd_amplitude(2), ylim_adhd_amplitude(2))];
    % Apply the global y-axis limits to both subplots
    subplot(1, 2, 1);
    ylim(global_ylim_amplitude);
    xticks(1:length(conditions));
    xticklabels(conditions);
    xlabel('Condition');
    
    subplot(1, 2, 2);
    ylim(global_ylim_amplitude);
    xticks(1:length(conditions));
    xticklabels(conditions);
    xlabel('Condition');
    sgtitle('Comparison of Saccade Amplitudes per Condition');

    if safe == 1
        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf, fullfile(comparison_results_folder, 'violin_saccade_amplitude_condition_group.png'));
    end

    
    
    
    %% Create figure with two subplots for velocity
    figure;

    % Subplot for Saccade Velocity - non-ADHD group
    subplot(1, 2, 1);
    [h1, L1, MX1, MED1, bw1] = violin(data_for_violin_velocity_nonadhd, 'xlabel', unique_conditions, 'facecolor', facecolors_nonadhd, 'edgecolor', 'none', 'facealpha', 0.5, 'mc', 'k', 'medc', 'k--');
    ylabel('log Saccade Velocity (pixels/ms)');
    xlabel('Condition');
    title('non-ADHD');
    % Capture the y-axis limits
    ylim_nonadhd_velocity = ylim();

    % Subplot for Saccade Velocity - ADHD group
    subplot(1, 2, 2);
    [h2, L2, MX2, MED2, bw2] = violin(data_for_violin_velocity_adhd, 'xlabel', unique_conditions, 'facecolor', facecolors_adhd, 'edgecolor', 'none', 'facealpha', 0.5, 'mc', 'k', 'medc', 'k--');
    ylabel('log Saccade Velocity (pixels/ms)');
    xlabel('Condition');
    title('ADHD');
    % Capture the y-axis limits
    ylim_adhd_velocity = ylim();

    % Determine the global y-axis limits
    global_ylim_velocity = [min(ylim_nonadhd_velocity(1), ylim_adhd_velocity(1)), max(ylim_nonadhd_velocity(2), ylim_adhd_velocity(2))];
    % Apply the global y-axis limits to both subplots
    subplot(1, 2, 1);
    ylim(global_ylim_velocity);
    xticks(1:length(conditions));
    xticklabels(conditions);
    xlabel('Condition');
    
    subplot(1, 2, 2);
    ylim(global_ylim_velocity);
    xticks(1:length(conditions));
    xticklabels(conditions);
    xlabel('Condition');
    sgtitle('Comparison of Saccade Velocities per Condition');

    if safe == 1
        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf, fullfile(comparison_results_folder, 'violin_saccade_velocity_condition_group.png'));
    end
end

