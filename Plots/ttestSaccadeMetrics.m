function ttestSaccadeMetrics(merged_data, condition, groups, adhd_indices, nonadhd_indices)

% Initialize arrays to hold data for each metric for condition A only
mean_first_saccade_start_adhd_a = [];
mean_first_saccade_start_nonadhd_a = [];
n_sac_adhd_a = [];
n_sac_nonadhd_a = [];
direction_difference_adhd_a = [];
direction_difference_nonadhd_a = [];
total_dist_sac_adhd_a = [];
total_dist_sac_nonadhd_a = [];
mean_dist_sac_adhd_a = [];
mean_dist_sac_nonadhd_a = [];
saccades_until_target_adhd_a = [];
saccades_until_target_nonadhd_a = [];
saccades_after_target_adhd_a = [];
saccades_after_target_nonadhd_a = [];
relevant_saccade_start_adhd_a = [];
relevant_saccade_start_nonadhd_a = [];

% Loop through each participant and extract data for condition
for i = 1:length(adhd_indices)
    if adhd_indices(i)
        % Condition indices for first group
        condition_a_indices = strcmp(merged_data.Condition{i}, condition);
        
        % Filter and calculate means for each metric for condition
        mean_first_saccade_start_adhd_a = [mean_first_saccade_start_adhd_a; ...
            nanmean(merged_data.first_saccade_start{i}(condition_a_indices))];
        n_sac_adhd_a = [n_sac_adhd_a; ...
            nanmean(merged_data.n_sac{i}(condition_a_indices))];
        direction_difference_adhd_a = [direction_difference_adhd_a; ...
            nanmean(merged_data.direction_difference{i}(condition_a_indices))];
        total_dist_sac_adhd_a = [total_dist_sac_adhd_a; ...
            nanmean(merged_data.total_dist_sac{i}(condition_a_indices))];
        mean_dist_sac_adhd_a = [mean_dist_sac_adhd_a; ...
            nanmean(merged_data.mean_dist_sac{i}(condition_a_indices))];
        saccades_until_target_adhd_a = [saccades_until_target_adhd_a; ...
            nanmean(merged_data.saccades_until_target{i}(condition_a_indices))];
        saccades_after_target_adhd_a = [saccades_after_target_adhd_a; ...
            nanmean(merged_data.saccades_after_target{i}(condition_a_indices))];
        relevant_saccade_start_adhd_a = [relevant_saccade_start_adhd_a; ...
            nanmean(merged_data.relevant_saccade_start{i}(condition_a_indices))];
        
    elseif nonadhd_indices(i)
        % Condition indices for second group participants
        condition_a_indices = strcmp(merged_data.Condition{i}, condition);
        
        % Filter and calculate means for each metric for condition
        mean_first_saccade_start_nonadhd_a = [mean_first_saccade_start_nonadhd_a; ...
            nanmean(merged_data.first_saccade_start{i}(condition_a_indices))];
        n_sac_nonadhd_a = [n_sac_nonadhd_a; ...
            nanmean(merged_data.n_sac{i}(condition_a_indices))];
        direction_difference_nonadhd_a = [direction_difference_nonadhd_a; ...
            nanmean(merged_data.direction_difference{i}(condition_a_indices))];
        total_dist_sac_nonadhd_a = [total_dist_sac_nonadhd_a; ...
            nanmean(merged_data.total_dist_sac{i}(condition_a_indices))];
        mean_dist_sac_nonadhd_a = [mean_dist_sac_nonadhd_a; ...
            nanmean(merged_data.mean_dist_sac{i}(condition_a_indices))];
        saccades_until_target_nonadhd_a = [saccades_until_target_nonadhd_a; ...
            nanmean(merged_data.saccades_until_target{i}(condition_a_indices))];
        saccades_after_target_nonadhd_a = [saccades_after_target_nonadhd_a; ...
            nanmean(merged_data.saccades_after_target{i}(condition_a_indices))];
        relevant_saccade_start_nonadhd_a = [relevant_saccade_start_nonadhd_a; ...
            nanmean(merged_data.relevant_saccade_start{i}(condition_a_indices))];
    end
end

% Perform t-tests for condition A on all metrics
[~, p_first_sacc_start_a] = ttest2(mean_first_saccade_start_adhd_a, mean_first_saccade_start_nonadhd_a);
fprintf('T-test p-values for %s Condition %s: first_saccade_start: p = %.3f\n', groups, condition, p_first_sacc_start_a);

[~, p_n_sacc_a] = ttest2(n_sac_adhd_a, n_sac_nonadhd_a);
fprintf('T-test p-values for %s Condition %s: number of saccades: p = %.3f\n', groups, condition, p_n_sacc_a);

[~, p_direction_difference_a] = ttest2(direction_difference_adhd_a, direction_difference_nonadhd_a);
fprintf('T-test p-values for %s Condition %s: Cosine Similarity: p = %.3f\n', groups, condition,p_direction_difference_a);

[~, p_total_dist_sac_a] = ttest2(total_dist_sac_adhd_a, total_dist_sac_nonadhd_a);
fprintf('T-test p-values for %s Condition %s: total_dist_sac: p = %.3f\n', groups, condition,p_total_dist_sac_a);

[~, p_mean_dist_sac_a] = ttest2(mean_dist_sac_adhd_a, mean_dist_sac_nonadhd_a);
fprintf('T-test p-values for %s Condition %s: mean_dist_sac: p = %.3f\n', groups, condition,p_mean_dist_sac_a);

[~, p_saccades_until_target_a] = ttest2(saccades_until_target_adhd_a, saccades_until_target_nonadhd_a);
fprintf('T-test p-values for %s Condition %s: saccades_until_target: p = %.3f\n', groups, condition,p_saccades_until_target_a);

[~, p_saccades_after_target_a] = ttest2(saccades_after_target_adhd_a, saccades_after_target_nonadhd_a);
fprintf('T-test p-values for %s Condition %s: saccades_after_target: p = %.3f\n', groups, condition,p_saccades_after_target_a);

[~, p_relevant_saccade_start_a] = ttest2(relevant_saccade_start_adhd_a, relevant_saccade_start_nonadhd_a);
fprintf('T-test p-values for %s Condition %s: relevant_saccade_start: p = %.3f\n',groups, condition, p_relevant_saccade_start_a);
