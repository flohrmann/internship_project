function plotQuestionnaires(ids, quest_struct, quest_table, group_labels, color_map, color_map_other, color_map_individual,comparison_results_folder, safe)
adhd_data = quest_struct.ADHD;
non_adhd_data = quest_struct.nonADHD;

% get questionnaire results/mean answers
quest_scores = questionnaireScale(quest_table, comparison_results_folder);

% Convert responses to numeric for easier analysis and plotting
adhd_numeric = adhd_data{:, 2:end};  % Skip the ID column
non_adhd_numeric = non_adhd_data{:, 2:end};  % Skip the ID column

plotQPieChart(adhd_numeric, non_adhd_numeric, safe, comparison_results_folder)

% Answers per question with a unique color for each participant
%plotQParticipantAnswers(quest_table, color_map, safe, comparison_results_folder);

% Means of answers per question per group
plotQMeanAndStd(adhd_numeric, non_adhd_numeric, color_map, safe, comparison_results_folder)

%%
disp('quest 6 : a,b,as,bs')
adhd_condition_avgs = adhd_numeric(:,1:6);
nonadhd_condition_avgs = non_adhd_numeric(:,1:6);
adhd_median    = median(adhd_condition_avgs, 1, 'omitnan')
adhd_sem       = std(adhd_condition_avgs, 1, 'omitnan')./ sqrt(size(adhd_condition_avgs, 1))
nonAdhd_median = median(nonadhd_condition_avgs, 1, 'omitnan')
nonAdhd_sem    = std(nonadhd_condition_avgs, 1, 'omitnan')./ sqrt(size(nonadhd_condition_avgs, 1))

plotADHDnonADHDandDiff('ASRS Results',...
    adhd_condition_avgs, adhd_median, adhd_sem, 'ADHD', 'northeast', ...
    nonadhd_condition_avgs, nonAdhd_median, nonAdhd_sem, 'nonADHD', 'northeast', ...
    ids, {'Q1','Q2','Q3','Q4','Q5','Q6'}, 'Answer (1-Never: 5:Very Often)', group_labels, {'1', '2', '3', '4', '5', '6'}, color_map_other, color_map_individual, ...
    fullfile(comparison_results_folder, '20_questionaire_6_allinone_median.png'));

%% 
disp('quest sum: a,b,as,bs')
adhd_condition_avgs = [quest_scores.TotalSumSymptoms(strcmp(group_labels, 'ADHD'), :), quest_scores.PartASumSymptoms(strcmp(group_labels, 'ADHD'), :)];
nonadhd_condition_avgs = [quest_scores.TotalSumSymptoms(strcmp(group_labels, 'nonADHD'), :), quest_scores.PartASumSymptoms(strcmp(group_labels, 'nonADHD'), :)];
adhd_median    = median(adhd_condition_avgs, 1, 'omitnan')
adhd_sem       = std(adhd_condition_avgs, 1, 'omitnan')./ sqrt(size(adhd_condition_avgs, 1))
nonAdhd_median = median(nonadhd_condition_avgs, 1, 'omitnan')
nonAdhd_sem    = std(nonadhd_condition_avgs, 1, 'omitnan')./ sqrt(size(nonadhd_condition_avgs, 1))

plotADHDnonADHDandDiff('Summed up ASRS Results',...
    adhd_condition_avgs, adhd_median, adhd_sem, 'ADHD', 'northeast', ...
    nonadhd_condition_avgs, nonAdhd_median, nonAdhd_sem, 'nonADHD', 'northeast', ...
    ids, {'All Questions', 'Part A'}, 'Summed up Questionaire Scores', group_labels, {'1', '2', '3', '4', '5', '6'}, color_map_other, color_map_individual, ...
    fullfile(comparison_results_folder, '20_questionaire_sums_allinone_median.png'));
