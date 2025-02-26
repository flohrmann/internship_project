function [R_group_condition, P_group_condition] = correlationByGroupCondition(data, groups, conditions, comparison_results_folder)

    % Initialize to store results for each group and condition
    R_group_condition = struct();
    P_group_condition = struct();
    
    % Loop through groups (ADHD, non-ADHD)
    for g = 1:length(groups)
        group = groups{g};
        
        % Initialize group-specific results
        R_group_condition.(group) = cell(1, length(conditions));
        P_group_condition.(group) = cell(1, length(conditions));

        % Loop through conditions (a, a simple, b, b simple)
        for c = 1:length(conditions)
            condition = conditions{c};
            
            % Initialize arrays for the current group and condition
            rt_button = [];
            rt_eye = [];
            accuracy = [];
            lapse = [];
            % Loop through participants and trials
            for participant = 1:length(data)
                if strcmp(data(participant).group, group)
                    for trial = 1:length(data(participant).Condition)
                        if strcmp(data(participant).Condition{trial}, condition)
                            % Collect reaction times and accuracy for this group and condition
                            rt_button = [rt_button; data(participant).rt(trial)];
                            rt_eye = [rt_eye; data(participant).rt_eye(trial)];
                            lapse = [lapse; data(participant).rt(trial) -data(participant).rt_eye(trial)];
                            accuracy = [accuracy; data(participant).accuracy(trial)];
                        end
                    end
                end
            end
            
            % Handle NaNs: remove rows where any data is NaN
            valid_idx = ~isnan(rt_button) & ~isnan(rt_eye) & ~isnan(accuracy);
            rt_button = rt_button(valid_idx);
            rt_eye = rt_eye(valid_idx);
            accuracy = accuracy(valid_idx);
            lapse = lapse(valid_idx);
            
            % Perform correlation analysis for the current group and condition
            [R, P] = corr([rt_button, rt_eye, accuracy, lapse], 'Type', 'Pearson');
            
            % Store results for the current group and condition
            R_group_condition.(group){c} = R;
            P_group_condition.(group){c} = P;
            
            % Display results
            disp(['Group: ', group, ', Condition: ', condition]);
            disp('Correlation Matrix (R):');
            disp(R);
            disp('P-values Matrix (P):');
            disp(P);
            disp('------------------------------------------');
        end
    end
    
    % Save results
    save(fullfile(comparison_results_folder, 'correlation_by_group_condition_results.mat'), 'R_group_condition', 'P_group_condition');
    disp('Group-specific and condition-specific correlation results saved.');
end
