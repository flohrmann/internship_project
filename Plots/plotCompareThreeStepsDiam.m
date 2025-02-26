function plotCompareTwoStepsDiam(data1, data_type_1, subtitle1, data2, data_type_2, subtitle2, bigtitle1, color_map, conditions, time_vector, x_label_text)
%% first data
    min_y = [];
    max_y = [];

num_rows = 1; num_col = 3;
figure; sgtitle(bigtitle1);

subplot(num_rows, num_col, 1); hold on; title(strcat('Average per Condition (',subtitle1,')'));
for i = 1:4
    condition = conditions{i};
    condition_trials = data1(strcmp(data1.Condition, condition), :);
    [y_upper, y_lower] =  plotConditionDataWithShading(condition_trials.(data_type_1), time_vector, color_map(condition), condition);
    min_y = [min_y, y_lower];
    max_y = [max_y, y_upper];
end
xline(0, '--', 't0', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
xlabel(x_label_text); ylabel('Pupil Diameter (mm)'); legend('show'); legend('boxoff');hold off;

%% 2nd data
subplot(num_rows, num_col, 2); hold on; title(strcat('Average per Condition Trials (',subtitle2,')'));
for i = 1:4
    condition = conditions{i}; 
    condition_trials = data2(strcmp(data2.Condition, condition), :);
    [y_upper, y_lower] =  plotConditionDataWithShading(condition_trials.(data_type_2), time_vector, color_map(condition), condition);
    min_y = [min_y, y_lower];
    max_y = [max_y, y_upper];
end
xline(0, '--', 't0', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
xlabel(x_label_text); ylabel('Pupil Diameter (mm)'); legend('show'); legend('boxoff');hold off;
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); %saveas(gcf, fullfile(analysis_folder, 'diam_t0_saccade_0_raw_0aligned_each_condition.png'));

%% 3rd data
subplot(num_rows, num_col, 3); hold on; title(strcat('Average per Condition Trials (',subtitle3,')'));
for i = 1:4
    condition = conditions{i}; 
    condition_trials = data3(strcmp(data3.Condition, condition), :);
    [y_upper, y_lower] = plotConditionDataWithShading(condition_trials.(data_type_3), time_vector, color_map(condition), condition);
    min_y = [min_y, y_lower];
    max_y = [max_y, y_upper];
end
xline(0, '--', 't0', 'LabelHorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
xlabel(x_label_text); ylabel('Pupil Diameter (mm)'); legend('show'); legend('boxoff');hold off;
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); %saveas(gcf, fullfile(analysis_folder, 'diam_t0_saccade_0_raw_0aligned_each_condition.png'));



y_limits = [min(min_y), max(max_y)];
% Apply unified y-limits and x-limits across individual trial subplots
for i = 1:4
    subplot(num_rows, num_col, i);
    ylim(y_limits);  % Set uniform y-axis
    xlim([min(time_vector), max(time_vector)]);  % Tight x-axis limits
end


