function correlationsSpeedAccuracyTradeOff(data_median, accuracy,rt_median_eye,rt_median_button, group_labels, groups, conditions, comparison_results_folder)

% averages over participants
% rtbutton, rteye, acc, lapse
[R_group_condition, P_group_condition] = correlationByGroupConditionAvg(accuracy,rt_median_eye,rt_median_button, group_labels, groups, conditions, comparison_results_folder);

% single trials (not as robust?)
% rtbutton, rteye, acc, lapse
%[R_group_condition, P_group_condition] = correlationByGroupCondition(data_median,  groups, conditions, comparison_results_folder);
% these return slightly different values since trials with no GazeRT have to be
% removed since NaN values cant be used for correlation analysis

end



function [R_group_condition, P_group_condition] = correlationByGroupConditionAvg(accuracy,rt_median_eye,rt_median_button, group_labels, groups, conditions, comparison_results_folder);
% adhd, nonadhd
which_index = [strcmp(group_labels, groups{1}), strcmp(group_labels, groups{2})];
% Loop through groups (ADHD, non-ADHD)
for g = 1:length(groups)
    group = groups{g};
    idxs = which_index(:,g);
    % Initialize group-specific results
    R_group_condition.(group) = cell(1, length(conditions));
    a= [];b= [];e= [];l = [];
    for c = 1:length(conditions)
        a = accuracy(idxs,c);
        b = rt_median_button(idxs,c);
        e = rt_median_eye(idxs,c);
        l = b-e;
        
        % Perform correlation analysis for the current group and condition
        [R, P] = corr([b, e, a, l], 'Type', 'Pearson');
        
        % Store results for the current group and condition
        R_group_condition.(group){c} = R;
        P_group_condition.(group){c} = P;
        
        % Display results
        disp(['Group: ', group, ', Condition: ', conditions{c}]);
        %disp(['gaze at target in: ', num2str(sum(valid_idx)/size(valid_idx, 1) *100) , ' % of trials'])
        disp('Correlation Matrix (R) rtbutton, rteye, acc, lapse:');
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
