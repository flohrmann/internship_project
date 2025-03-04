function [mistakes, accuracy] = plotMeanAccuracyPerGroupCondition(data, group_labels, ids, unique_conditions, condition_labels, color_map,color_map_individual, comp_results_fix, safe)



%% same but only for simple vs nonsimple condition
numParticipants = length(data); % Number of participants
accuracy = NaN(numParticipants, length(unique_conditions)/2);
mistakes = NaN(numParticipants, length(unique_conditions)/2);

for p = 1:numParticipants
        
        %nonnsimple
        condition_trials_a = strcmp(data(p).Condition, unique_conditions{1});
        acc_cond_a = data(p).accuracy(condition_trials_a);
        
        condition_trials_b = strcmp(data(p).Condition, unique_conditions{2});
        acc_cond_b = data(p).accuracy(condition_trials_b);
        
        accs_nonsimple = [acc_cond_a;acc_cond_b];
        mistakes(p, 1) = size(accs_nonsimple,1) - sum(accs_nonsimple);
        accuracy(p, 1) = sum(accs_nonsimple)/size(accs_nonsimple,1);
        
        %simple
        condition_trials_as = strcmp(data(p).Condition, unique_conditions{3});
        acc_cond_as = data(p).accuracy(condition_trials_as);
        
        condition_trials_bs = strcmp(data(p).Condition, unique_conditions{4});
        acc_cond_bs = data(p).accuracy(condition_trials_bs);
        
        accs_simple = [acc_cond_as;acc_cond_bs];
        mistakes(p, 2) = size(accs_simple,1) - sum(accs_simple);
        accuracy(p, 2) = sum(accs_simple)/size(accs_simple,1);
end


% num mistakes
adhd_condition_avgs    = mistakes(strcmp(group_labels, 'ADHD'), :);
nonadhd_condition_avgs = mistakes(strcmp(group_labels, 'nonADHD'), :);
adhd_median    = median(adhd_condition_avgs, 1, 'omitnan');
adhd_sem       = std(adhd_condition_avgs, 1, 'omitnan')./ sqrt(size(adhd_condition_avgs, 1));
nonAdhd_median = median(nonadhd_condition_avgs, 1, 'omitnan');
nonAdhd_sem    = std(nonadhd_condition_avgs, 1, 'omitnan')./ sqrt(size(nonadhd_condition_avgs, 1));

plotHEREADHDnonADHDVariance('Mistakes Per Condition',...
    adhd_condition_avgs, adhd_median, adhd_sem, 'ADHD', 'northeast', ...
    nonadhd_condition_avgs, nonAdhd_median, nonAdhd_sem, 'nonADHD', 'east', ...
    ids, {"nonSimple", "simple"}, 'Number of Wrong Trials', group_labels, unique_conditions, color_map, color_map_individual, ...
    fullfile(comp_results_fix, '11_mistakes_allinone_median_simplevsnonsimple.png'));

% accuracy
adhd_condition_avgs    = accuracy(strcmp(group_labels, 'ADHD'), :);
nonadhd_condition_avgs = accuracy(strcmp(group_labels, 'nonADHD'), :);
adhd_median    = median(adhd_condition_avgs, 1, 'omitnan');
adhd_sem       = std(adhd_condition_avgs, 1, 'omitnan')./ sqrt(size(adhd_condition_avgs, 1));
nonAdhd_median = median(nonadhd_condition_avgs, 1, 'omitnan');
nonAdhd_sem    = std(nonadhd_condition_avgs, 1, 'omitnan')./ sqrt(size(nonadhd_condition_avgs, 1));

plotHEREADHDnonADHDVariance('Accuracy Per Condition',...
    adhd_condition_avgs, adhd_median, adhd_sem, 'ADHD', 'southeast', ...
    nonadhd_condition_avgs, nonAdhd_median, nonAdhd_sem, 'nonADHD', 'east', ...
    ids, {"nonSimple", "simple"}, '% of correct trials', group_labels, unique_conditions, color_map, color_map_individual, ...
    fullfile(comp_results_fix, '11_accuracy_allinone_median_simplevsnonsimple.png'));



%% a vs rest
numParticipants = length(data); % Number of participants
accuracy = NaN(numParticipants, length(unique_conditions)/2);
mistakes = NaN(numParticipants, length(unique_conditions)/2);
accs_nonsimple = [];
accs_simple = [];
for p = 1:numParticipants
        %nonnsimple
        condition_trials_a = strcmp(data(p).Condition, unique_conditions{1});
        acc_cond_a = data(p).accuracy(condition_trials_a);
        
        accs_nonsimple = [acc_cond_a];
        mistakes(p, 1) = size(accs_nonsimple,1) - sum(accs_nonsimple);
        accuracy(p, 1) = sum(accs_nonsimple)/size(accs_nonsimple,1);
        
        %simple
        condition_trials_b = strcmp(data(p).Condition, unique_conditions{2});
        acc_cond_b = data(p).accuracy(condition_trials_b);
        
        condition_trials_as = strcmp(data(p).Condition, unique_conditions{3});
        acc_cond_as = data(p).accuracy(condition_trials_as);
        
        condition_trials_bs = strcmp(data(p).Condition, unique_conditions{4});
        acc_cond_bs = data(p).accuracy(condition_trials_bs);
        
        accs_simple = [acc_cond_as;acc_cond_bs; acc_cond_b];
        mistakes(p, 2) = size(accs_simple,1) - sum(accs_simple);
        accuracy(p, 2) = sum(accs_simple)/size(accs_simple,1);
end
% num mistakes
adhd_condition_avgs    = mistakes(strcmp(group_labels, 'ADHD'), :);
nonadhd_condition_avgs = mistakes(strcmp(group_labels, 'nonADHD'), :);
adhd_median    = median(adhd_condition_avgs, 1, 'omitnan');
adhd_sem       = std(adhd_condition_avgs, 1, 'omitnan')./ sqrt(size(adhd_condition_avgs, 1));
nonAdhd_median = median(nonadhd_condition_avgs, 1, 'omitnan');
nonAdhd_sem    = std(nonadhd_condition_avgs, 1, 'omitnan')./ sqrt(size(nonadhd_condition_avgs, 1));

plotHEREADHDnonADHDVariance('Mistakes Per Condition',...
    adhd_condition_avgs, adhd_median, adhd_sem, 'ADHD', 'northeast', ...
    nonadhd_condition_avgs, nonAdhd_median, nonAdhd_sem, 'nonADHD', 'east', ...
    ids, {"a", "b, a simple, b simple"}, 'Number of Wrong Trials', group_labels, unique_conditions, color_map, color_map_individual, ...
    fullfile(comp_results_fix, '11_mistakes_allinone_median_avsrest.png'));

% accuracy
adhd_condition_avgs    = accuracy(strcmp(group_labels, 'ADHD'), :);
nonadhd_condition_avgs = accuracy(strcmp(group_labels, 'nonADHD'), :);
adhd_median    = median(adhd_condition_avgs, 1, 'omitnan');
adhd_sem       = std(adhd_condition_avgs, 1, 'omitnan')./ sqrt(size(adhd_condition_avgs, 1));
nonAdhd_median = median(nonadhd_condition_avgs, 1, 'omitnan');
nonAdhd_sem    = std(nonadhd_condition_avgs, 1, 'omitnan')./ sqrt(size(nonadhd_condition_avgs, 1));

plotHEREADHDnonADHDVariance('Accuracy Per Condition',...
    adhd_condition_avgs, adhd_median, adhd_sem, 'ADHD', 'southeast', ...
    nonadhd_condition_avgs, nonAdhd_median, nonAdhd_sem, 'nonADHD', 'east', ...
    ids, {"a", "b, a simple, b simple"}, '% of correct trials', group_labels, unique_conditions, color_map, color_map_individual, ...
    fullfile(comp_results_fix, '11_accuracy_allinone_median_avsrest.png'));


%% 4 conditions
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
        
        
        mistakes(p, c) = size(acc_cond,1) - sum(acc_cond);
        accuracy(p, c) = sum(acc_cond)/size(acc_cond,1);
    end
end


% num mistakes
adhd_condition_avgs    = mistakes(strcmp(group_labels, 'ADHD'), :);
nonadhd_condition_avgs = mistakes(strcmp(group_labels, 'nonADHD'), :);
adhd_median    = median(adhd_condition_avgs, 1, 'omitnan');
adhd_sem       = std(adhd_condition_avgs, 1, 'omitnan')./ sqrt(size(adhd_condition_avgs, 1));
nonAdhd_median = median(nonadhd_condition_avgs, 1, 'omitnan');
nonAdhd_sem    = std(nonadhd_condition_avgs, 1, 'omitnan')./ sqrt(size(nonadhd_condition_avgs, 1));

plotADHDnonADHDandDiff('Mistakes Per Condition',...
    adhd_condition_avgs, adhd_median, adhd_sem, 'ADHD', 'east', ...
    nonadhd_condition_avgs, nonAdhd_median, nonAdhd_sem, 'nonADHD', 'east', ...
    ids, condition_labels, 'Number of Wrong Trials', group_labels, unique_conditions, color_map, color_map_individual, ...
    fullfile(comp_results_fix, '11_mistakes_allinone_median.png'));

% accuracy
adhd_condition_avgs    = accuracy(strcmp(group_labels, 'ADHD'), :);
nonadhd_condition_avgs = accuracy(strcmp(group_labels, 'nonADHD'), :);
adhd_median    = median(adhd_condition_avgs, 1, 'omitnan');
adhd_sem       = std(adhd_condition_avgs, 1, 'omitnan')./ sqrt(size(adhd_condition_avgs, 1));
nonAdhd_median = median(nonadhd_condition_avgs, 1, 'omitnan');
nonAdhd_sem    = std(nonadhd_condition_avgs, 1, 'omitnan')./ sqrt(size(nonadhd_condition_avgs, 1));

plotADHDnonADHDandDiff('Accuracy Per Condition',...
    adhd_condition_avgs, adhd_median, adhd_sem, 'ADHD', 'southeast', ...
    nonadhd_condition_avgs, nonAdhd_median, nonAdhd_sem, 'nonADHD', 'east', ...
    ids, condition_labels, '% of correct trials', group_labels, unique_conditions, color_map, color_map_individual, ...
    fullfile(comp_results_fix, '11_accuracy_allinone_median.png'));

end










function plotHEREADHDnonADHDVariance(bigtitle,... % sgtitle
                        data_1, group_1_medians, group_1_sems, label_1, legend_1_place, ...  % adhd data
                        data_2, group_2_medians, group_2_sems, label_2, legend_2_place, ...  % nonadhd data
                        ids,...     % all participant ids
                        x_labels, y_label, ...% x, y axis labels
                        group_labels, conditions, color_map, participant_map, ...
                        safe_name)
% data_1&2 : [participants x conditions]
ymax = max(max(data_1(:)), max(data_2(:)))+ max(max(group_1_sems(:)), max(group_2_sems(:))); % Maximum individual median value
ymin = min(min(data_1(:)), min(data_2(:)))- max(max(group_1_sems(:)), max(group_2_sems(:)))/2; % Maximum individual median value

sharedYLimits = [ymin, ymax];

figure; sgtitle(bigtitle);
% ADHD Plot
subplot(1, 3, 1);
plotHEREStandardTwoGroupsEachData(data_1, group_1_medians, group_1_sems,label_1, y_label, legend_1_place, ids, group_labels, conditions, x_labels, color_map, participant_map, sharedYLimits);

% nonADHD Plot
subplot(1, 3, 2);
plotHEREStandardTwoGroupsEachData(data_2, group_2_medians, group_2_sems, label_2, y_label, legend_2_place, ids, group_labels, conditions, x_labels, color_map, participant_map, sharedYLimits);

% Variance Difference Plot
subplot(1, 3, 3); 
plotHEREVarianceDifference(data_1, group_1_sems, data_2, group_2_sems, conditions, x_labels, color_map, label_2);
  
set(gcf, 'Position', [50, 50, 1400, 700]); % Resize the figure window (x, y, width, height)
saveas(gcf, safe_name);

% end 
end



function plotHEREStandardTwoGroupsEachData(data, medians, sems, group, y_label, position, ids, group_labels, conditions, labels, color_map, participant_map, sharedYLimits)
    % Plot group medians and scatter individual participant data
    numConditions = length(conditions)/2;
    hold on;
    
    % Bar plot for medians with error bars
    for c = 1:numConditions
        bar(c, medians(c), 'FaceColor', color_map(conditions{c}), 'EdgeColor', 'none', ...
            'FaceAlpha', 0.3, 'DisplayName', group, 'HandleVisibility', 'off');
    end
    
    current_indices = strcmp(group_labels, group);
    group_ids = ids(current_indices);

    % Scatter individual participant data
    jitterAmount = 0.4; 
    for participant = 1:length(group_ids)
        participant_id = group_ids(participant);

        %if strcmp(group_labels{participant}, group)
            info = participant_map(participant_id); % to retrieve color and marker
            current_data = data(participant,:);
            jitteredX = (1:numConditions) + jitterAmount * (rand(1, numConditions) - 0.5); % Jitter x-coordinates
            
            scatter(jitteredX, current_data, 40, ... % Use scatter for points only
                'MarkerEdgeColor', info.color, ...
                'MarkerFaceColor', info.color, ...
                'Marker', info.marker, ...
                'DisplayName', sprintf('ID%d', participant_id));
        %end
    end
    
    for c = 1:numConditions % errorbar on top for better visibility
       errorbar(c, medians(c), sems(c), '.', 'LineWidth', 1.5, 'CapSize', 10, 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    end
    xticks(1:numConditions); xticklabels(labels);
    ylabel(y_label);  ylim(sharedYLimits);
    legend('boxoff');
    legend('Location', position, 'Orientation', 'horizontal', 'NumColumns', 2); % Two-column legend
    title(group);    hold off;
end

%%
function plotHEREVarianceDifference(adhd_medians, adhd_sem, nonadhd_medians, nonAdhd_sem, conditions, labels, color_map, y_label)
    % Calculate variance metric and perform T-tests or permutation tests
    difference = median(adhd_medians) - median(nonadhd_medians);
    % difference in standard error of the mean 
    difference_sem = sqrt(adhd_sem.^2 + nonAdhd_sem.^2) ;

    % Perform statistical tests and store p-values
    numConditions = length(conditions)/2;
    p_values = nan(1, numConditions); % Store p-values for annotation
    legendEntries = cell(1, numConditions); % Store legend labels
    
    hold on;
    for c = 1:numConditions
        % Extract medians for each condition
        adhd_data = adhd_medians(:, c); % ADHD group medians for condition c
        nonadhd_data = nonadhd_medians(:, c); % nonADHD group medians for condition c
        
        % Perform a T-test (or permutation test if needed)
        [~, p] = ttest2(adhd_data, nonadhd_data);
        p_values(c) = p; % Store p-value
        
        % Plot variance metric
        bar(c, difference(c), 'FaceColor', color_map(conditions{c}), 'EdgeColor', 'none', ...
            'FaceAlpha', 0.3, 'DisplayName', sprintf('%s Variance', conditions{c}));
        
        errorbar(c, difference(c), difference_sem(c), '.', 'LineWidth', 1.5, 'CapSize', 10, 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
        legendEntries{c} = sprintf('p=%.2f', p);
    end
    
    title('ADHD vs. nonADHD');
    xticks(1:numConditions); xticklabels(labels);
    ylabel(strcat('Difference in',{' '}, y_label))

    lgd = legend(legendEntries, 'Location', 'southoutside');
    legend('boxoff'); %lgd.ItemTextAlignment = 'left'; 
    title(lgd, 'ADHD/nonADHD Medians t-test'); 
    hold off;
end
