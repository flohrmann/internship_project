function [mistakes, accuracy] = plotMeanAccuracyPerGroupCondition(data, group_labels, ids, unique_conditions, condition_labels, color_map,color_map_individual, comp_results_fix, safe)

numParticipants = length(data); % Number of participants
accuracy = zeros(numParticipants, length(unique_conditions));
mistakes = zeros(numParticipants, length(unique_conditions));
% Loop through each participant
for p = 1:numParticipants
    % Loop through each condition
    for c = 1:length(unique_conditions)
        condition = unique_conditions{c};
        condition_trials = strcmp(data(p).Condition, condition);

        acc_cond = data(p).accuracy(condition_trials);
        
        
        mistakes(p, c) = height(acc_cond) - sum(acc_cond);
        accuracy(p, c) = sum(acc_cond)/height(acc_cond);
    end
end


%% num mistakes
adhd_condition_avgs    = mistakes(strcmp(group_labels, 'ADHD'), :);
nonadhd_condition_avgs = mistakes(strcmp(group_labels, 'nonADHD'), :);
adhd_median    = median(adhd_condition_avgs, 1, 'omitnan');
adhd_sem       = std(adhd_condition_avgs, 1, 'omitnan')./ sqrt(size(adhd_condition_avgs, 1));
nonAdhd_median = median(nonadhd_condition_avgs, 1, 'omitnan');
nonAdhd_sem    = std(nonadhd_condition_avgs, 1, 'omitnan')./ sqrt(size(nonadhd_condition_avgs, 1));

plotADHDnonADHDVariance('Mistakes Per Condition',...
    adhd_condition_avgs, adhd_median, adhd_sem, 'ADHD', 'northeast', ...
    nonadhd_condition_avgs, nonAdhd_median, nonAdhd_sem, 'nonADHD', 'northeast', ...
    ids, condition_labels, 'Number of Wrong Trials', group_labels, unique_conditions, color_map, color_map_individual, ...
    fullfile(comp_results_fix, '11_mistakes_allinone_median.png'));

%% % accuracy
adhd_condition_avgs    = accuracy(strcmp(group_labels, 'ADHD'), :);
nonadhd_condition_avgs = accuracy(strcmp(group_labels, 'nonADHD'), :);
adhd_median    = median(adhd_condition_avgs, 1, 'omitnan');
adhd_sem       = std(adhd_condition_avgs, 1, 'omitnan')./ sqrt(size(adhd_condition_avgs, 1));
nonAdhd_median = median(nonadhd_condition_avgs, 1, 'omitnan');
nonAdhd_sem    = std(nonadhd_condition_avgs, 1, 'omitnan')./ sqrt(size(nonadhd_condition_avgs, 1));

plotADHDnonADHDVariance('Accuracy Per Condition',...
    adhd_condition_avgs, adhd_median, adhd_sem, 'ADHD', 'southeast', ...
    nonadhd_condition_avgs, nonAdhd_median, nonAdhd_sem, 'nonADHD', 'southeast', ...
    ids, condition_labels, '% of correct trials', group_labels, unique_conditions, color_map, color_map_individual, ...
    fullfile(comp_results_fix, '11_accuracy_allinone_median.png'));

end
