function plotAccuracyVsButtonPressRT(data, mistakes, accuracy, rt_median_eye, rt_median_button, ...
                                    group_labels, ids, unique_conditions, condition_labels, color_map,color_map_individual, comparison_results_folder, safe)

%% data by group                                
adhd_eye    = rt_median_eye(strcmp(group_labels, 'ADHD'), :);
nonadhd_eye = rt_median_eye(strcmp(group_labels, 'nonADHD'), :);
                               
adhd_button    = rt_median_button(strcmp(group_labels, 'ADHD'), :);
nonadhd_button = rt_median_button(strcmp(group_labels, 'nonADHD'), :);

adhd_mistakes    = mistakes(strcmp(group_labels, 'ADHD'), :);
nonadhd_mistakes = mistakes(strcmp(group_labels, 'nonADHD'), :);

adhd_accuracy    = accuracy(strcmp(group_labels, 'ADHD'), :);
nonadhd_accuracy = accuracy(strcmp(group_labels, 'nonADHD'), :);

num_participants = length(group_labels);
num_conditions = length(unique_conditions);


%% accuracy anova
% % Define the Within-Subjects Design for 4 conditions
% conditions = {'a', 'b', 'a_simple', 'b_simple'}; % Match column order in accuracy
% within = table(conditions', 'VariableNames', {'Condition'});
% % Create RM-ANOVA Table (Using Separate Columns for Each Condition)
% tbl_rm = table(group_labels, accuracy(:,1), accuracy(:,2), accuracy(:,3), accuracy(:,4), ...
%                'VariableNames', {'Group', 'a', 'b', 'a_simple', 'b_simple'});
% % Convert categorical variables
% tbl_rm.Group = categorical(tbl_rm.Group);
% % Fit Repeated-Measures Model for Accuracy
% rm_accuracy = fitrm(tbl_rm, 'a,b,a_simple,b_simple ~ Group', 'WithinDesign', within);
% ranova_accuracy = ranova(rm_accuracy);
% disp('Repeated-Measures ANOVA (Accuracy)');
% disp(ranova_accuracy);
% % Run post-hoc comparisons for Condition
% multcompare(rm_accuracy, 'Condition', 'ComparisonType', 'bonferroni')

%% Linear Mixed-Effects Model (LME)
ids = repelem(1:num_participants, num_conditions, 1)';
groups = repelem(group_labels, 1, num_conditions);
cond = repelem(1:num_conditions,num_participants, 1);
% Create table for analysis
tbl = table(ids(:), groups(:), cond(:), ...
                accuracy(:), rt_median_eye(:), rt_median_button(:), ...
                'VariableNames', {'Participant', 'Group', 'Condition', 'Accuracy', 'RT_eye', 'RT_button'});

% Convert categorical variables
tbl.Group = categorical(tbl.Group);
tbl.Condition = categorical(tbl.Condition);

% Fit LME model: Does RT predict Accuracy?
lme_eye = fitlme(tbl, 'Accuracy ~ RT_eye * Group * Condition + (1|Participant)');
disp('Linear Mixed-Effects Model (RT Eye ~ Accuracy)');
disp(lme_eye);

lme_button = fitlme(tbl, 'Accuracy ~ RT_button * Group * Condition + (1|Participant)');
disp('Linear Mixed-Effects Model (RT Button ~ Accuracy)');
disp(lme_button);







%% spearman correlation
% eye
% [r_eye, p_eye] = corr(rt_median_eye(:), accuracy(:), 'Type', 'Spearman');
% fprintf('Spearman correlation (All RT_eye vs Accuracy): r = %.3f, p = %.3f\n', r_eye, p_eye);
% [r_eye_adhd, p_eye_adhd] = corr(adhd_eye(:), adhd_accuracy(:), 'Type', 'Spearman');
% fprintf('Spearman correlation (ADHD RT_eye vs Accuracy): r = %.3f, p = %.3f\n', r_eye_adhd, p_eye_adhd);
% [r_eye_non, p_eye_non] = corr(nonadhd_eye(:), nonadhd_accuracy(:), 'Type', 'Spearman');
% fprintf('Spearman correlation (nonDHD RT_eye vs Accuracy): r = %.3f, p = %.3f\n', r_eye_non, p_eye_non);
% 
% % button 
% [r_button, p_button] = corr(rt_median_button(:), accuracy(:), 'Type', 'Spearman');
% fprintf('Spearman correlation (All RT_button vs Accuracy): r = %.3f, p = %.3f\n', r_button, p_button);
% [r_button_adhd, p_button_adhd] = corr(adhd_button(:), adhd_accuracy(:), 'Type', 'Spearman');
% fprintf('Spearman correlation (ADHD RT_button vs Accuracy): r = %.3f, p = %.3f\n', r_button_adhd, p_button_adhd);
% [r_button_non, p_button_non] = corr(nonadhd_button(:), nonadhd_accuracy(:), 'Type', 'Spearman');
% fprintf('Spearman correlation (nonADHD RT_button vs Accuracy): r = %.3f, p = %.3f\n', r_button_non, p_button_non);





       
                                
                                
end
