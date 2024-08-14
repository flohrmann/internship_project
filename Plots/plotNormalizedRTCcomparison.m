function plotNormalizedRTCcomparison(data, color_map, comparison_results_folder, safe)

rt_adhd = [];
rt_nonadhd = [];
for i = 1:length(data)
    if strcmp(data(i).group, 'ADHD')
        rt_adhd = [rt_adhd; data(i).normalized_rt];
    elseif strcmp(data(i).group, 'nonADHD')
        rt_nonadhd = [rt_nonadhd; data(i).normalized_rt];
    end
end

figure;
% Plot boxcharts for ADHD group
boxchart(ones(size(rt_adhd)), rt_adhd, 'BoxFaceColor', color_map('ADHD'), 'BoxWidth', 0.5);
hold on;

% Plot boxcharts for non-ADHD group with x-offset
boxchart(2 * ones(size(rt_nonadhd)), rt_nonadhd, 'BoxFaceColor', color_map('nonADHD'), 'BoxWidth', 0.5);
hold off;

set(gca, 'XTick', [1, 2], 'XTickLabel', {'ADHD', 'nonADHD'});
xlabel('Group');
set(gca, 'YScale', 'log');
ylabel('Reaction Time (s)');
title('Normalized Reaction Time Comparison Between Groups');
if safe == 1
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf, fullfile(comparison_results_folder, 'boxchart_group_normalized_rt.png'));
else
end