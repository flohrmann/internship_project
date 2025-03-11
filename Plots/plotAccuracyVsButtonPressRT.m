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


%% anova
%[p,tbl,stats] = anova1(rt_median_button(:), mistakes(:));
%disp(tbl)

%% eye
ymax = max(accuracy(:)); 
ymin = min(accuracy(:));
xmax = max(rt_median_eye(:)); 
xmin = min(rt_median_eye(:)); 
sharedYLimits = [ymin, ymax];
sharedXLimits = [xmin, xmax];

figure;
sgtitle('Gaze-RT vs Accuracy by Condition [Participant Medians and Regression Line per Condition]');

% ADHD Subplot
[r_1, p_1] = corr(adhd_eye(:), adhd_accuracy(:), 'Type', 'Spearman');
fprintf('adhd Spearman Correlation eye x acc: r = %.2f, p = %.3f\n', r_1, p_1)

subplot(1, 2, 1); hold on;
legendEntries_ADHD = cell(1, num_conditions);
for c = 1:num_conditions
    % Extract RT and Accuracy values for ADHD group
    rt_data = adhd_eye(:, c); 
    acc_data = adhd_accuracy(:, c); 
    % Scatter plot for ADHD group
    scatterHandles_ADHD(c) = scatter(rt_data, acc_data, 50, color_map(unique_conditions{c}), 'filled', 'DisplayName', condition_labels{c});
    % Store legend entry
    legendEntries_ADHD{c} = sprintf('%s', condition_labels{c});
    % Fit regression line
    coeffs = polyfit(rt_data, acc_data, 1);
    x_fit = linspace(min(rt_data), max(rt_data), 100);
    y_fit = polyval(coeffs, x_fit);
    plot(x_fit, y_fit, 'Color', color_map(unique_conditions{c}), 'LineWidth', 1.5);
end
xlabel('Reaction Time (s)');ylabel('Accuracy (%)');
%title(sprintf('ADHD - Spearman Correlation: r = %.2f, p = %.3f\n', r_1, p_1));
title('ADHD');
legend(scatterHandles_ADHD, legendEntries_ADHD, 'Location', 'southwest', 'Interpreter', 'none');
ylim(sharedYLimits);xlim(sharedXLimits);
hold off;

% Non-ADHD Subplot
[r_2, p_2] = corr(nonadhd_eye(:), nonadhd_accuracy(:), 'Type', 'Spearman');
fprintf('nonadhd Spearman Correlation eye x acc: r = %.2f, p = %.3f\n', r_2, p_2)
subplot(1, 2, 2); hold on;
legendEntries_NonADHD = cell(1, num_conditions);

for c = 1:num_conditions
    % Extract RT and Accuracy values for Non-ADHD group
    rt_data = nonadhd_eye(:, c); 
    acc_data = nonadhd_accuracy(:, c); 
    % Scatter plot for Non-ADHD group
    scatterHandles_nonADHD(c) = scatter(rt_data, acc_data, 50, color_map(unique_conditions{c}), 'filled', 'DisplayName', condition_labels{c});
    % Store legend entry
    legendEntries_NonADHD{c} = sprintf('%s', condition_labels{c}); 
    % Fit regression line
    coeffs = polyfit(rt_data, acc_data, 1);
    x_fit = linspace(min(rt_data), max(rt_data), 100);
    y_fit = polyval(coeffs, x_fit);
    plot(x_fit, y_fit, 'Color', color_map(unique_conditions{c}), 'LineWidth', 1.5);
end
xlabel('Reaction Time (s)');ylabel('Accuracy (%)');
%title(sprintf('nonADHD - Spearman Correlation: r = %.2f, p = %.3f\n', r_2, p_2));
title('nonADHD');
legend(scatterHandles_nonADHD, legendEntries_NonADHD, 'Location', 'southwest', 'Interpreter', 'none');
ylim(sharedYLimits);xlim(sharedXLimits);
hold off;

set(gcf, 'Position', [50, 50, 1200, 600]); 
saveas(gcf, strcat(comparison_results_folder, "\11_accuracy_rt_eye.png"));


%% button press
xmax = max(adhd_button(:)); 
xmin = min(nonadhd_button(:)); 
sharedYLimits = [ymin, ymax];
sharedXLimits = [xmin, xmax];

figure;
sgtitle('Button-Press-RT vs Accuracy by Condition [Participant Medians and Regression Line per Condition]');

% ADHD Subplot
[r_3, p_3] = corr(adhd_button(:), adhd_accuracy(:), 'Type', 'Spearman');
fprintf('adhd Spearman Correlation button x acc: r = %.2f, p = %.3f\n', r_3, p_3)

subplot(1, 2, 1); hold on;
legendEntries_ADHD = cell(1, num_conditions);
for c = 1:num_conditions
    % Extract RT and Accuracy values for ADHD group
    rt_data = adhd_button(:, c); 
    acc_data = adhd_accuracy(:, c); 
    % Scatter plot for ADHD group
    scatterHandles_ADHD(c) = scatter(rt_data, acc_data, 50, color_map(unique_conditions{c}), 'filled', 'DisplayName', condition_labels{c});
    % Store legend entry
    legendEntries_ADHD{c} = sprintf('%s', condition_labels{c});
    % Fit regression line
    coeffs = polyfit(rt_data, acc_data, 1);
    x_fit = linspace(min(rt_data), max(rt_data), 100);
    y_fit = polyval(coeffs, x_fit);
    plot(x_fit, y_fit, 'Color', color_map(unique_conditions{c}), 'LineWidth', 1.5);
end
xlabel('Reaction Time (s)');ylabel('Accuracy (%)');
%title(sprintf('ADHD - Spearman Correlation: r = %.2f, p = %.3f\n', r_3, p_3));
title('ADHD');
legend(scatterHandles_ADHD, legendEntries_ADHD, 'Location', 'best', 'Interpreter', 'none');
ylim(sharedYLimits);xlim(sharedXLimits);
hold off;

% Non-ADHD Subplot
[r_spearman, p_spearman] = corr(nonadhd_button(:), nonadhd_accuracy(:), 'Type', 'Spearman');
fprintf('nonadhd Spearman Correlation button x acc: r = %.2f, p = %.3f\n', r_spearman, p_spearman)
subplot(1, 2, 2); hold on;
legendEntries_NonADHD = cell(1, num_conditions);

for c = 1:num_conditions
    % Extract RT and Accuracy values for Non-ADHD group
    rt_data = nonadhd_button(:, c); 
    acc_data = nonadhd_accuracy(:, c); 
    % Scatter plot for Non-ADHD group
    scatterHandles_nonADHD(c) = scatter(rt_data, acc_data, 50, color_map(unique_conditions{c}), 'filled', 'DisplayName', condition_labels{c});
    % Store legend entry
    legendEntries_NonADHD{c} = sprintf('%s', condition_labels{c}); 
    % Fit regression line
    coeffs = polyfit(rt_data, acc_data, 1);
    x_fit = linspace(min(rt_data), max(rt_data), 100);
    y_fit = polyval(coeffs, x_fit);
    plot(x_fit, y_fit, 'Color', color_map(unique_conditions{c}), 'LineWidth', 1.5);
end
xlabel('Reaction Time (s)');ylabel('Accuracy (%)');
%title(sprintf('nonADHD - Spearman Correlation: r = %.2f, p = %.3f\n', r_spearman, p_spearman));
title('nonADHD');
legend(scatterHandles_nonADHD, legendEntries_NonADHD, 'Location', 'best', 'Interpreter', 'none');
ylim(sharedYLimits);xlim(sharedXLimits);
hold off;

set(gcf, 'Position', [50, 50, 1200, 600]); 
saveas(gcf, strcat(comparison_results_folder, "\11_accuracy_rt_buttonpress.png"));


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

       
                                
                                
end
