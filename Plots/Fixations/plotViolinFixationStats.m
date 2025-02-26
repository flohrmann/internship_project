function plotViolinFixationStats(fixationStats, group_labels, conditions, condition_labels, color_map, color_map_individual, comparison_results_folder, safe)
    % Initialize containers for ADHD and non-ADHD groups
    mean_adhd_conditions = [];
    mean_nonadhd_conditions = [];
    sem_adhd_conditions = [];
    sem_nonadhd_conditions = [];
    all_means_condition = [];

    % Loop through the struct to gather mean and SEM data for each group
    for i = 1:size(fixationStats, 2)
        conditionData = fixationStats(i);
        all_means_condition = [conditionData.a_mean, conditionData.as_mean, ...
                               conditionData.b_mean, conditionData.bs_mean];
        if strcmp(group_labels{i}, 'ADHD')
            % Append means and SEMs for ADHD group
            mean_values = [conditionData.a_mean, conditionData.as_mean, ...
                           conditionData.b_mean, conditionData.bs_mean];
            sem_values = [conditionData.a_sem, conditionData.as_sem, ...
                          conditionData.b_sem, conditionData.bs_mean];
            mean_adhd_conditions = [mean_adhd_conditions; mean_values];
            sem_adhd_conditions = [sem_adhd_conditions; sem_values];
        elseif strcmp(group_labels{i}, 'nonADHD')
            % Append means and SEMs for non-ADHD group
            mean_values = [conditionData.a_mean, conditionData.as_mean, ...
                           conditionData.b_mean, conditionData.bs_mean];
            sem_values = [conditionData.a_sem, conditionData.as_sem, ...
                          conditionData.b_sem, conditionData.bs_mean];
            mean_nonadhd_conditions = [mean_nonadhd_conditions; mean_values];
            sem_nonadhd_conditions = [sem_nonadhd_conditions; sem_values];
        end
    end

    % Prepare data for violin plots
    num_conditions = size(mean_adhd_conditions, 2); % Number of conditions
    data_for_violin_mean_adhd = cell(1, num_conditions);
    data_for_violin_mean_nonadhd = cell(1, num_conditions);

    for i = 1:num_conditions
        % Extract mean data for each condition
        data_for_violin_mean_adhd{i} = mean_adhd_conditions(:, i);
        data_for_violin_mean_nonadhd{i} = mean_nonadhd_conditions(:, i);
    end

    %% Create figure with two subplots for fixation duration
    figure;
   
    % Subplot for Fixation Duration - non-ADHD group
    ax1 = subplot(1, 2, 2);
    violin(data_for_violin_mean_nonadhd, 'xlabel', condition_labels, 'facecolor', repmat(color_map('nonADHD'), num_conditions, 1), ...
           'edgecolor', 'none', 'facealpha', 0.5, 'mc', 'k', 'medc', 'k--');
    ylabel('Fixation Duration (s)');
    title('nonADHD');
    hold on;

    % Add SEM as error bars for non-ADHD group
    for i = 1:num_conditions
        x = repmat(i, size(mean_nonadhd_conditions, 1), 1);
        hError = errorbar(x, mean_nonadhd_conditions(:, i), sem_nonadhd_conditions(:, i), 'k.', 'CapSize', 0); % Example of creating an error bar
        set(get(get(hError, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end

    % Add individual markers for non-ADHD group using color_map_individual
    %% TODO fix individual markers
    jitterAmount = 0.1; % Adjust for marker spread
    for participant = 1:size(all_means_condition, 1)
        for condition = 1:num_conditions
            x = condition + (rand(1) - 0.5) * jitterAmount; % Jitter x-coordinates
            scatter(x, all_means_condition(participant, condition), 60, ...
                'MarkerFaceColor', color_map_individual(participant).color, ... % Use individual color for nonADHD
                'MarkerEdgeColor', color_map_individual(participant).color, ...
                'Marker', color_map_individual(participant).marker, ...
                'HandleVisibility', 'off'); % Exclude individual markers from legend
        end
    end
    hold off;

    
    
    
    % Subplot for Fixation Duration - ADHD group
    ax2 = subplot(1, 2, 1);
    violin(data_for_violin_mean_adhd, 'xlabel', condition_labels, 'facecolor', repmat(color_map('ADHD'), num_conditions, 1), ...
           'edgecolor', 'none', 'facealpha', 0.5, 'mc', 'k', 'medc', 'k--');
    ylabel('Fixation Duration (s)');
    title('ADHD');
    hold on;

    % Add SEM as error bars for ADHD group
    for i = 1:num_conditions
        x = repmat(i, size(mean_adhd_conditions, 1), 1);
        hError = errorbar(x, mean_adhd_conditions(:, i), sem_adhd_conditions(:, i), 'k.', 'CapSize', 0); % Example of creating an error bar
        set(get(get(hError, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end

    % Add individual markers for ADHD group using color_map_individual
    for participant = 1:size(mean_adhd_conditions, 1)
        
        for condition = 1:num_conditions
            x = condition + (rand(1) - 0.5) * jitterAmount; % Jitter x-coordinates
            scatter(x, mean_adhd_conditions(participant, condition), 60, ...
                'MarkerFaceColor', color_map_individual(participant), ... % Use individual color for ADHD
                'MarkerEdgeColor', 'k', ...
                'Marker', 'o', ...
                'HandleVisibility', 'off'); % Exclude individual markers from legend
        end
    end
    hold off;   
    
    linkaxes([ax1, ax2], 'y');    % Link y-axes of the two subplots
    sgtitle('Comparison of Fixation Durations per Condition');
    if safe == 1
        set(gcf, 'Position', [0 0 1 1]);
        saveas(gcf, fullfile(comparison_results_folder, '01_fixation_violin_fixation_duration_condition_group.png'));
    end
end



% function plotViolinFixationStats(fixationStats, group_labels, conditions, condition_labels, color_map, comparison_results_folder, safe)
%     % Initialize containers for ADHD and non-ADHD groups
%     mean_adhd_conditions = [];
%     mean_nonadhd_conditions = [];
%     sem_adhd_conditions = [];
%     sem_nonadhd_conditions = [];
%     
%     % Loop through the struct to gather mean and SEM data for each group
%     for i = 1:size(fixationStats, 2)
%         conditionData = fixationStats(i);
%         if strcmp(group_labels{i}, 'ADHD')
%             % Append means and SEMs for ADHD group
%             mean_values = [conditionData.a_mean, conditionData.as_mean, ...
%                            conditionData.b_mean, conditionData.bs_mean];
%             sem_values = [conditionData.a_sem, conditionData.as_sem, ...
%                           conditionData.b_sem, conditionData.bs_sem];
%             mean_adhd_conditions = [mean_adhd_conditions; mean_values];
%             sem_adhd_conditions = [sem_adhd_conditions; sem_values];
%         elseif strcmp(group_labels{i}, 'nonADHD')
%             % Append means and SEMs for non-ADHD group
%             mean_values = [conditionData.a_mean, conditionData.as_mean, ...
%                            conditionData.b_mean, conditionData.bs_mean];
%             sem_values = [conditionData.a_sem, conditionData.as_sem, ...
%                           conditionData.b_sem, conditionData.bs_sem];
%             mean_nonadhd_conditions = [mean_nonadhd_conditions; mean_values];
%             sem_nonadhd_conditions = [sem_nonadhd_conditions; sem_values];
%         end
%     end
% 
%     % Prepare data for violin plots
%     num_conditions = size(mean_adhd_conditions, 2); % Number of conditions
%     data_for_violin_mean_adhd = cell(1, num_conditions);
%     data_for_violin_mean_nonadhd = cell(1, num_conditions);
%     facecolors_adhd = zeros(num_conditions, 3);
%     facecolors_nonadhd = zeros(num_conditions, 3);
% 
%     for i = 1:num_conditions
%         % Extract mean data for each condition
%         data_for_violin_mean_adhd{i} = mean_adhd_conditions(:, i);
%         data_for_violin_mean_nonadhd{i} = mean_nonadhd_conditions(:, i);
%         
%         % Set face colors
%         facecolors_adhd(i, :) = color_map(conditions{i});
%         facecolors_nonadhd(i, :) = color_map(conditions{i});
%     end
% 
%     %% Create figure with two subplots for fixation duration
%     figure;
%    
%     % Subplot for Fixation Duration - non-ADHD group
%     ax1 = subplot(1, 2, 2);
%     
%     violin(data_for_violin_mean_nonadhd, 'xlabel', {'a', 'a simple', 'b', 'b simple'}, 'facecolor', facecolors_nonadhd, ...
%            'edgecolor', 'none', 'facealpha', 0.5, 'mc', 'k', 'medc', 'k--');
%     ylabel('Fixation Duration (s)');
%     %xlabel(condition_labels);
%     title('nonADHD');
%     % Add SEM as error bars for non-ADHD group
%     hold on;
%     for i = 1:num_conditions
%         x = repmat(i, size(mean_nonadhd_conditions, 1), 1);
%         hError = errorbar(x, mean_nonadhd_conditions(:, i), sem_nonadhd_conditions(:, i), 'k.', 'CapSize', 0); % Example of creating an error bar
%         set(get(get(hError, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
%     end
%     hold off;
% 
%     % Subplot for Fixation Duration - ADHD group
%     ax2 = subplot(1, 2, 1);
%     violin(data_for_violin_mean_adhd, 'xlabel', {'a', 'a simple', 'b', 'b simple'}, 'facecolor', facecolors_adhd, ...
%            'edgecolor', 'none', 'facealpha', 0.5, 'mc', 'k', 'medc', 'k--');
%     ylabel('Fixation Duration (s)');
%     title('ADHD');
%     hold on;
%     for i = 1:num_conditions    % Add SEM as error bars for ADHD group
%         x = repmat(i, size(mean_adhd_conditions, 1), 1);
%         hError = errorbar(x, mean_adhd_conditions(:, i), sem_adhd_conditions(:, i), 'k.', 'CapSize', 0); % Example of creating an error bar
%         set(get(get(hError, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
%     end
%     hold off;   
%     
%     linkaxes([ax1, ax2], 'y');    % Link y-axes of the two subplots
%     sgtitle('Comparison of Fixation Durations per Condition');
%     if safe == 1
%         set(gcf, 'Position', [0 0 1 1]);
%         saveas(gcf, fullfile(comparison_results_folder, '01_fixation_violin_fixation_duration_condition_group.png'));
%     end
% end