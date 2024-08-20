function plotViolinMeanRTGroupCondition(data, color_map, comparison_results_folder, safe)

rt_adhd_conditions = [];
condition_adhd = {};
rt_nonadhd_conditions = [];
condition_nonadhd = {};

% Loop through the struct to gather data for ADHD and non-ADHD groups
for i = 1:length(data)
    if strcmp(data(i).group, 'ADHD')
        % Append the normalized RTs for each condition
        rt_adhd_conditions = [rt_adhd_conditions; data(i).nRTa; data(i).nRTb; data(i).nRTasimple; data(i).nRTbsimple];
        condition_adhd = [condition_adhd; repmat({'a'}, length(data(i).nRTa), 1); ...
                                      repmat({'b'}, length(data(i).nRTb), 1); ...
                                      repmat({'a_simple'}, length(data(i).nRTasimple), 1); ...
                                      repmat({'b_simple'}, length(data(i).nRTbsimple), 1)];
    elseif strcmp(data(i).group, 'nonADHD')
        % Append the normalized RTs for each condition
        rt_nonadhd_conditions = [rt_nonadhd_conditions; data(i).nRTa; data(i).nRTb; data(i).nRTasimple; data(i).nRTbsimple];
        condition_nonadhd = [condition_nonadhd; repmat({'a'}, length(data(i).nRTa), 1); ...
                                           repmat({'b'}, length(data(i).nRTb), 1); ...
                                           repmat({'a_simple'}, length(data(i).nRTasimple), 1); ...
                                           repmat({'b_simple'}, length(data(i).nRTbsimple), 1)];
    end
end

% Convert the condition labels to categorical and extract unique conditions
condition_adhd = categorical(condition_adhd);
condition_nonadhd = categorical(condition_nonadhd);
unique_conditions = categories(condition_adhd); % Assuming both groups have the same conditions

% Prepare data for ADHD and non-ADHD violin plots
data_for_violin_adhd = cell(1, length(unique_conditions));
data_for_violin_nonadhd = cell(1, length(unique_conditions));
facecolors_adhd = zeros(length(unique_conditions), 3);
facecolors_nonadhd = zeros(length(unique_conditions), 3);

for i = 1:length(unique_conditions)
    data_for_violin_adhd{i} = log10(rt_adhd_conditions(condition_adhd == unique_conditions{i}));
    %data_for_violin_adhd{i} = log10(rt_adhd_conditions(condition_adhd == unique_conditions{i}));
    data_for_violin_nonadhd{i} = log10(rt_nonadhd_conditions(condition_nonadhd == unique_conditions{i}));
    %data_for_violin_nonadhd{i} = log10(rt_nonadhd_conditions(condition_nonadhd == unique_conditions{i}));
    facecolors_adhd(i, :) = color_map(char(unique_conditions{i}));
    facecolors_nonadhd(i, :) = color_map(char(unique_conditions{i}));
end

% Create figure with two subplots
figure;

% Subplot for non-ADHD group
subplot(1, 2, 1);
[h1, L1, MX1, MED1, bw1] = violin(data_for_violin_nonadhd, 'xlabel', unique_conditions, 'facecolor', facecolors_nonadhd, 'edgecolor', 'none', 'facealpha', 0.5, 'mc', 'k', 'medc', 'k--');
ylabel('log(Reaction Time (s))');
xlabel('Condition');
title('non-ADHD');

% Capture the y-axis limits
ylim_nonadhd = ylim();

% Subplot for ADHD group
subplot(1, 2, 2);
[h2, L2, MX2, MED2, bw2] = violin(data_for_violin_adhd, 'xlabel', unique_conditions, 'facecolor', facecolors_adhd, 'edgecolor', 'none', 'facealpha', 0.5, 'mc', 'k', 'medc', 'k--');
ylabel('log(Reaction Time (s))');
xlabel('Condition');
title('ADHD');

% Capture the y-axis limits
ylim_adhd = ylim();

% Determine the global y-axis limits
global_ylim = [min(ylim_nonadhd(1), ylim_adhd(1)), max(ylim_nonadhd(2), ylim_adhd(2))];

% Apply the global y-axis limits to both subplots
subplot(1, 2, 1);
ylim(global_ylim);

subplot(1, 2, 2);
ylim(global_ylim);


sgtitle('Comparison of Normalized Mean RTs per Condition');
if safe == 1
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf, fullfile(comparison_results_folder, 'violin_condition_group_normalized_rt.png'));
end

end
