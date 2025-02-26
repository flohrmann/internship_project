function plotDiamsPerConditionAndAverage(PLOTDATA, data_type, conditions, condition_labels, time_vector, x_label_text,color_map, bigtitle, ...
                                         xlinelabel, analysis_folder, compare_folder, t0,id)

num_rows = 3; 
num_col = 2; 
%figure; 
%sgtitle(bigtitle);
all_y_values = [];

t = tiledlayout(num_rows, num_col, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Add shared title and axis labels
title(t,bigtitle)
xlabel(t,x_label_text)
ylabel(t,'Pupil Diameter (mm)')

for i = 1:4
    %subplot(num_rows, num_col, i); 
    ax = nexttile;  % Move to the next tile
    hold on;
    
    condition = conditions{i}; 
    title(strcat('Individual Trials: ', condition_labels{i}));
    condition_trials = PLOTDATA(strcmp(PLOTDATA.Condition, condition), :);
    for j = 1:height(condition_trials)
        try
            y_data = condition_trials.(data_type)(j,:);
            plot(time_vector, y_data, 'Color', color_map(condition));
            all_y_values = [all_y_values; [max(y_data), min(y_data)]];  % Collect y-values for unified axis
        catch 
        end
    end
    xline(0, '--', xlinelabel,  'Color', [0.3 0.3 0.3], 'LabelOrientation', 'horizontal', 'HandleVisibility', 'off');
    %xlabel(x_label_text); ylabel('Pupil Diameter (mm)');
    hold off;
end


y_limits = [min(min(all_y_values)), max(max(all_y_values))];
% Apply unified y-limits and x-limits across individual trial subplots
for i = 1:4
    ax = nexttile(i);
    ylim(ax, y_limits);
    xlim(ax, [min(time_vector), max(time_vector)]);
%     subplot(num_rows, num_col, i);
%     ylim(y_limits);  % Set uniform y-axis
%     xlim([min(time_vector), max(time_vector)]);  % Tight x-axis limits
end

% Plot the condition average in a larger tile (spanning width of layout)
ax_avg = nexttile([1, 2]);  % Span across 2 columns
%subplot(num_rows, num_col, 5); 
hold on; title('Condition Average');
for i = 1:4
    condition = conditions{i}; 
    condition_trials = PLOTDATA(strcmp(PLOTDATA.Condition, condition), :);
    plotConditionDataWithShading(condition_trials.(data_type), time_vector, color_map(condition), condition);
end
xline(0, '--', xlinelabel, 'LabelHorizontalAlignment', 'right', 'LabelOrientation', 'horizontal', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
%xlabel(x_label_text); ylabel('Pupil Diameter (mm)'); 

lgd = legend(condition_labels); 
set(lgd, 'Orientation', 'horizontal', 'Location', 'southoutside', 'Box', 'off');

%ylim(y_limits);  % Apply same y-limits to the average plot
%xlim([min(time_vector), max(time_vector)]);  % Tight x-axis limits


% Apply consistent y-limits and x-limits to the average plot
%ylim(ax_avg, y_limits);
xlim(ax_avg, [min(time_vector), max(time_vector)]);

hold off;


%set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf, fullfile(analysis_folder, strcat('diam_',t0,'_0aligned_each_condition.png')));
print(gcf, fullfile(compare_folder, strcat('diam_',t0,'_zero_aligned_subj',num2str(id),'.svg')), '-dsvg');

