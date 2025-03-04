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
    
    % Loop through conditions (a, b, a simple, b simple)
    for c = 1:length(conditions)
        
        % Initialize arrays for the current group and condition
        rt_button = [];
        rt_eye = [];
        accuracy = [];
        lapse = [];
        avg_rt_button = [];
        avg_rt_eye = [];
        avg_accuracy = [];
        avg_lapse = [];
        
        % Loop through participants and trials
        for participant = 1:length(data)
            
            if strcmp(data(participant).group, group)
                current_data = data(participant);
                condition_trials = strcmp(current_data.Condition, conditions{c});
                
                % Collect reaction times and accuracy for this group and condition
                rt_button = [rt_button; current_data.rt(condition_trials)];
                rt_eye = [rt_eye; current_data.rt_eye(condition_trials)];
                lapse = [lapse; current_data.rt(condition_trials) - current_data.rt_eye(condition_trials)];
                accuracy = [accuracy; current_data.accuracy(condition_trials)];
                
                
                % medians per participant
                a1 = current_data.accuracy(condition_trials);
                b1 = current_data.rt(condition_trials);
                e1 = current_data.rt_eye(condition_trials);
                
                valid_idx = ~isnan(a1) & ~isnan(b1) & ~isnan(e1);

                avg_rt_button = [avg_rt_button; nanmedian(b1(valid_idx))];
                avg_rt_eye = [avg_rt_eye; nanmedian(e1(valid_idx))];
                avg_lapse = [avg_lapse; nanmedian(b1(valid_idx)) - nanmedian(e1(valid_idx))];
                avg_accuracy = [avg_accuracy; nanmedian(a1(valid_idx))];
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
        [R, P] = corr([avg_rt_button, avg_rt_eye, avg_accuracy, avg_lapse], 'Type', 'Pearson');
     
        % Store results for the current group and condition
        R_group_condition.(group){c} = R;
        P_group_condition.(group){c} = P;
        
        % Display results
        disp(['Group: ', group, ', Condition: ', conditions{c}]);
        disp(['gaze at target in: ', num2str(sum(valid_idx)/size(valid_idx, 1) *100) , ' % of trials'])
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






% function [R_group_condition, P_group_condition] = correlationByGroupCondition(accuracy,rt_median_eye,rt_median_button, group_labels, groups, conditions, comparison_results_folder)
% % adhd, nonadhd
% which_index = [strcmp(group_labels, groups{1}), strcmp(group_labels, groups{2})];
% % Loop through groups (ADHD, non-ADHD)
% for g = 1:length(groups)
%     group = groups{g};
%     idxs = which_index(:,g);
%     % Initialize group-specific results
%     R_group_condition.(group) = cell(1, length(conditions));
%     a= [];b= [];e= [];l = [];
%     for c = 1:length(conditions)
%         a = accuracy(idxs,c);
%         b = rt_median_button(idxs,c);
%         e = rt_median_eye(idxs,c);
%         l = b-e;
%         
%         % Perform correlation analysis for the current group and condition
%         [R, P] = corr([b, e, a, l], 'Type', 'Pearson');
%         
%         % Store results for the current group and condition
%         R_group_condition.(group){c} = R;
%         P_group_condition.(group){c} = P;
%         
%         % Display results
%         disp(['Group: ', group, ', Condition: ', conditions{c}]);
%         %disp(['gaze at target in: ', num2str(sum(valid_idx)/size(valid_idx, 1) *100) , ' % of trials'])
%         disp('Correlation Matrix (R) rtbutton, rteye, acc, lapse:');
%         disp(R);
%         disp('P-values Matrix (P):');
%         disp(P);
%         disp('------------------------------------------');
%         
%     end
%     
%     
% end
% 
% % Save results
% save(fullfile(comparison_results_folder, 'correlation_by_group_condition_results.mat'), 'R_group_condition', 'P_group_condition');
% disp('Group-specific and condition-specific correlation results saved.');
% end
