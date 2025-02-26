function repeatedMeasuresANOVA(avg_adhd, avg_nonadhd, comparison_results_folder)
    % Group labels: 1 for ADHD, 2 for Non-ADHD
    % Define conditions and pupil measures
    conditions = {'a', 'a_simple', 'b', 'b_simple'};
    measures = {'min', 'max', 'mean'};  % Pupil measures to analyze

    combined_data = [avg_adhd; avg_nonadhd];
    group_labels = [ones(size(avg_adhd, 1), 1); 2 * ones(size(avg_nonadhd, 1), 1)];  % 1 = ADHD, 2 = Non-ADHD
    % Extract data for each measure (min, max, mean)
    pupil_min = [];  
    pupil_max = [];  
    pupil_mean = [];
    
    % Loop through subjects to get data across conditions
    for subj_idx = 1:height(combined_data)
        % Get the data for the subject and append to matrices
        pupil_min(subj_idx, :) = [combined_data.a_min(subj_idx), combined_data.as_min(subj_idx), ...
                                  combined_data.b_min(subj_idx), combined_data.bs_min(subj_idx)];
        pupil_max(subj_idx, :) = [combined_data.a_max(subj_idx), combined_data.as_max(subj_idx), ...
                                  combined_data.b_max(subj_idx), combined_data.bs_max(subj_idx)];
        pupil_mean(subj_idx, :) = [combined_data.a_mean(subj_idx), combined_data.as_mean(subj_idx), ...
                                   combined_data.b_mean(subj_idx), combined_data.bs_mean(subj_idx)];
    end

    % Group variable (ADHD = 1, Non-ADHD = 2)
    group_var = group_labels;

    % Define condition labels for repeated measures
    condition_factor = table({'a'; 'a_simple'; 'b'; 'b_simple'}, 'VariableNames', {'Condition'});

    % Create a table that combines all data and group labels
    pupil_min_table = array2table(pupil_min, 'VariableNames', {'a', 'a_simple', 'b', 'b_simple'});
    pupil_max_table = array2table(pupil_max, 'VariableNames', {'a', 'a_simple', 'b', 'b_simple'});
    pupil_mean_table = array2table(pupil_mean, 'VariableNames', {'a', 'a_simple', 'b', 'b_simple'});

    % Add group labels to the table
    pupil_min_table.Group = group_var;
    pupil_max_table.Group = group_var;
    pupil_mean_table.Group = group_var;

    % Run Repeated Measures ANOVA for each measure
    fprintf('Running Repeated Measures ANOVA for Min Pupil Diameter...\n');
    rm_anova_min = fitrm(pupil_min_table, 'a-b_simple ~ Group', 'WithinDesign', condition_factor);
    ranova_min_results = ranova(rm_anova_min);
    disp(ranova_min_results);

    fprintf('Running Repeated Measures ANOVA for Max Pupil Diameter...\n');
    rm_anova_max = fitrm(pupil_max_table, 'a-b_simple ~ Group', 'WithinDesign', condition_factor);
    ranova_max_results = ranova(rm_anova_max);
    disp(ranova_max_results);

    fprintf('Running Repeated Measures ANOVA for Mean Pupil Diameter...\n');
    rm_anova_mean = fitrm(pupil_mean_table, 'a-b_simple ~ Group', 'WithinDesign', condition_factor);
    ranova_mean_results = ranova(rm_anova_mean);
    disp(ranova_mean_results);

    save(fullfile(comparison_results_folder, 'anova_min_results.mat'), 'ranova_min_results');
    save(fullfile(comparison_results_folder, 'anova_max_results.mat'), 'ranova_max_results');
    save(fullfile(comparison_results_folder, 'anova_mean_results.mat'), 'ranova_mean_results');
end