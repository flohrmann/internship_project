function plotCompareTwoStepsDiam(data1, data_type_1, subtitle1, ...
                                 data2, data_type_2, subtitle2, ...
                                 bigtitle1, color_map, conditions, condition_labels, time_vector, x_label_text, t0_label, ...
                                 analysis_folder, compare_folder, t0, id)
    %% Initialize variables
    all_y_values = [];
    num_rows = 2; 
    num_col = 2;

    %% Create the figure
    figure;
    sgtitle(bigtitle1);

    %% First data - Individual trials
    subplot(num_rows, num_col, 1); 
    hold on; 
    title(strcat('Individual Trials (', subtitle1, ')'));
    for i = 1:4
        condition = conditions{i};
        condition_trials = data1(strcmp(data1.Condition, condition), :);
        for j = 1:height(condition_trials)
            y_data = condition_trials.(data_type_1)(j, :);     
            plot(time_vector, y_data, 'Color', color_map(condition));
            all_y_values = [all_y_values; [max(y_data), min(y_data)]];
        end
    end
    xline(0, '--', t0_label, 'LabelHorizontalAlignment', 'right', 'LabelOrientation', 'horizontal', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    xlabel(x_label_text); 
    ylabel('Pupil Diameter (mm)');
    hold off;

    %% First data - Condition averages
    subplot(num_rows, num_col, 2); 
    hold on; 
    title(strcat('Condition Average (', subtitle1, ')'));
    for i = 1:4
        condition = conditions{i};
        condition_trials = data1(strcmp(data1.Condition, condition), :);
        plotConditionDataWithShading(condition_trials.(data_type_1), time_vector, color_map(condition), condition);
    end
    xline(0, '--', t0_label, 'LabelHorizontalAlignment', 'right', 'LabelOrientation', 'horizontal', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    xlabel(x_label_text); 
    ylabel('Pupil Diameter (mm)');
    hold off;

    %% Second data - Individual trials
    subplot(num_rows, num_col, 3); 
    hold on; 
    title(strcat('Individual Trials (', subtitle2, ')'));
    for i = 1:4
        condition = conditions{i};
        condition_trials = data2(strcmp(data2.Condition, condition), :);
        for j = 1:height(condition_trials)
            y_data = condition_trials.(data_type_2)(j, :);     
            plot(time_vector, y_data, 'Color', color_map(condition));
            all_y_values = [all_y_values; [max(y_data), min(y_data)]];
        end
    end
    xline(0, '--', t0_label, 'LabelHorizontalAlignment', 'right', 'LabelOrientation', 'horizontal', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    xlabel(x_label_text); 
    ylabel('Pupil Diameter (mm)');
    hold off;

    %% Second data - Condition averages
    subplot(num_rows, num_col, 4); 
    hold on; 
    title(strcat('Condition Average (', subtitle2, ')'));
    for i = 1:4
        condition = conditions{i};
        condition_trials = data2(strcmp(data2.Condition, condition), :);
        plotConditionDataWithShading(condition_trials.(data_type_2), time_vector, color_map(condition), condition);
    end
    xline(0, '--', t0_label, 'LabelHorizontalAlignment', 'right', 'LabelOrientation', 'horizontal', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    xlabel(x_label_text); 
    ylabel('Pupil Diameter (mm)');
    hold off;

    %% Apply unified y-limits and x-limits
    y_limits = [min(all_y_values, [], 'all'), max(all_y_values, [], 'all')];
    for i = 1:4
        subplot(num_rows, num_col, i);
        ylim(y_limits);  % Set uniform y-axis
        xlim([min(time_vector), max(time_vector)]);  % Tight x-axis limits
    end

    %% Add a single shared legend at the bottom
    lgd = legend(condition_labels, 'Orientation', 'horizontal', 'Location', 'southoutside', 'Box', 'off');
    lgd.Position = [0.35, -0.01, 0.3, 0.05]; % Adjust this to change position (normalized units)

    %% Save the figure
    saveas(gcf, fullfile(analysis_folder, strcat('diam_t0_', t0, '.png')));
    print(gcf, fullfile(compare_folder, strcat('diam_t0_', t0, '_id', num2str(id), '.svg')), '-dsvg');
end














% function plotCompareTwoStepsDiam(data1, data_type_1, subtitle1, ...
%                                  data2, data_type_2, subtitle2, ...
%                                  bigtitle1, color_map, conditions, condition_labels, time_vector, x_label_text, t0_label, ...
%                                  analysis_folder, compare_folder, t0, id)
% %% first data
% all_y_values = [];
% avg_y_values = [];
% 
% num_rows = 2; num_col = 2;
% figure; sgtitle(bigtitle1);
% subplot(num_rows, num_col, 1); hold on; title(strcat('Individual Trials (',subtitle1,')'));
% for i = 1:4
%     condition = conditions{i};
%     condition_trials = data1(strcmp(data1.Condition, condition), :);
%     for j = 1:height(condition_trials)
%         y_data = condition_trials.(data_type_1)(j,:);     
%         plot(time_vector, y_data, 'Color', color_map(condition));
%         all_y_values = [all_y_values; [max(y_data), min(y_data)]];  % Collect y-values for unified axis
%     end
% end
% xline(0, '--', t0_label, 'LabelHorizontalAlignment', 'right','LabelOrientation', 'horizontal',  'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
% xlabel(x_label_text); ylabel('Pupil Diameter (mm)'); hold off;
% %
% subplot(num_rows, num_col, 2); hold on; title(strcat('Condition Average (',subtitle1,')'));
% for i = 1:4
%     condition = conditions{i};
%     condition_trials = data1(strcmp(data1.Condition, condition), :);
%     plotConditionDataWithShading(condition_trials.(data_type_1), time_vector, color_map(condition), condition);
% end
% % xline_handle = 
% xline(0, '--', t0_label, 'LabelHorizontalAlignment', 'right', 'LabelOrientation', 'horizontal', 'Color', [0.3 0.3 0.3]);
% xlabel(x_label_text); ylabel('Pupil Diameter (mm)'); legend('show'); legend('boxoff');hold off;
% %set(xline_handle, 'DisplayName', t0_label);
% legend('off');
% %% 2nd data
% for i = 1:4
%     subplot(num_rows, num_col, 3); hold on;
%     condition = conditions{i}; title(strcat('Individual Trials (',subtitle2,')'));
%     condition_trials = data2(strcmp(data2.Condition, condition), :);
%     for j = 1:height(condition_trials)
%         y_data = condition_trials.(data_type_2)(j,:);     
%         plot(time_vector, y_data, 'Color', color_map(condition));
%         all_y_values = [all_y_values; [max(y_data), min(y_data)]];  % Collect y-values for unified axis
%     end
% end
% xline(0, '--', t0_label, 'LabelHorizontalAlignment', 'right', 'LabelOrientation', 'horizontal', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
% xlabel(x_label_text); ylabel('Pupil Diameter Change (mm)');
% hold off;
% %
% subplot(num_rows, num_col, 4); hold on; title(strcat('Condition Average(',subtitle2,')'));
% for i = 1:4
%     condition = conditions{i}; 
%     condition_trials = data2(strcmp(data2.Condition, condition), :);
%     plotConditionDataWithShading(condition_trials.(data_type_2), time_vector, color_map(condition), condition);
% end
% xline(0, '--', t0_label, 'LabelHorizontalAlignment', 'right', 'LabelOrientation', 'horizontal', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
% xlabel(x_label_text); ylabel('Pupil Diameter (mm)');% legend('show'); legend('boxoff'); 
% hold off;
% %set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); %saveas(gcf, fullfile(analysis_folder, 'diam_t0_saccade_0_raw_0aligned_each_condition.png'));
% 
% y_limits = [min(all_y_values, [], 'all'), max(all_y_values, [], 'all')];
% % Apply unified y-limits and x-limits across individual trial subplots
% for i = 1:4
%     subplot(num_rows, num_col, i);
%     ylim(y_limits);  % Set uniform y-axis
%     xlim([min(time_vector), max(time_vector)]);  % Tight x-axis limits
% end
% 
%  
% lgd = legend(condition_labels); 
% set(lgd, 'Orientation', 'horizontal', 'Location', 'southoutside', 'Box', 'off');
% 
% saveas(gcf, fullfile(analysis_folder, strcat('diam_t0_',t0,'.png')));
% print(gcf, fullfile(compare_folder, strcat('diam_t0_',t0,'_id',num2str(id),'.svg')), '-dsvg');
