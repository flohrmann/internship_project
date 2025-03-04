function plotBarNumFixations(fixationStats, group_labels, ids,  conditions, condition_labels, color_map, color_map_individual, comparison_results_folder, safe)
    % Initialize containers for ADHD and non-ADHD groups
    avg_adhd_conditions = [];
    avg_nonadhd_conditions = [];
    adhd_subject_ids = {};
    nonadhd_subject_ids = {};
    
    % Loop through participants to calculate avg number of fixations per condition
    for i = 1:length(fixationStats)
        trials = fixationStats(i).trials;
        numFixations = [trials.numFixations]; % Extract number of fixations
        trialConditions = {trials.conditions}; % Extract conditions for each trial
        
        % Aggregate data by condition
        participant_avgs = zeros(1, length(conditions));
        for c = 1:length(conditions)
            % Extract strings from the nested cells
            trialConditionsStrings = cellfun(@(x) x{1}, trialConditions, 'UniformOutput', false);
            
            % Find indices matching the current condition
            condition_idx = strcmp(trialConditionsStrings, conditions{c});
            
            % Compute the mean number of fixations for this condition
            participant_avgs(c) = median(numFixations(condition_idx)); % median for this condition
        end
        
        % Assign to ADHD or non-ADHD group
        if strcmp(group_labels{i}, 'ADHD')
            avg_adhd_conditions = [avg_adhd_conditions; participant_avgs];
            adhd_subject_ids = [adhd_subject_ids; fixationStats(i).id]; % Store subject IDs
        elseif strcmp(group_labels{i}, 'nonADHD')
            avg_nonadhd_conditions = [avg_nonadhd_conditions; participant_avgs];
            nonadhd_subject_ids = [nonadhd_subject_ids; fixationStats(i).id]; % Store subject IDs
        end
    end

    % group-level avg and SEM
    avg_adhd = median(avg_adhd_conditions, 1);
    sem_adhd = std(avg_adhd_conditions, [], 1) ./ sqrt(size(avg_adhd_conditions, 1));
    
    avg_nonadhd = median(avg_nonadhd_conditions, 1);
    sem_nonadhd = std(avg_nonadhd_conditions, [], 1) ./ sqrt(size(avg_nonadhd_conditions, 1));

    plotADHDnonADHDandDiff('Number of Fixations per Trial',... % sgtitle
                            avg_adhd_conditions, avg_adhd, sem_adhd, 'ADHD', 'northeast', ...  % adhd data
                            avg_nonadhd_conditions, avg_nonadhd, sem_nonadhd, 'nonADHD', 'northeast', ...  % nonadhd data
                            ids,condition_labels, 'Median Fixations per Trial', ...% x, y axis labels
                            group_labels, conditions, color_map, color_map_individual,...
                            fullfile(comparison_results_folder, '01_fixation_bar_mean_num_fixations_allinone.png'))

end

