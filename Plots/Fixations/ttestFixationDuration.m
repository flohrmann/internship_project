function ttestFixationDuration(group_labels, fixationStats, saccadeStats, conditions)
%%
% Assume `fixationStats` and `saccadeStats` are structs that hold your fixation and saccade data
% Assume `group_labels` holds the ADHD status for each participant
% Assume `conditions` is an array of your experimental conditions

% Initialize containers for ADHD and non-ADHD data
adhdFixationDurations = [];
adhdSaccadeAmplitudes = [];
adhdSaccadeCounts = [];

nonAdhdFixationDurations = [];
nonAdhdSaccadeAmplitudes = [];
nonAdhdSaccadeCounts = [];

% Loop through participants and segregate data
for participant = 1:length(group_labels)
    if strcmp(group_labels{participant}, 'ADHD')
        % Aggregate ADHD data
        for condition = 1:length(conditions)
            adhdFixationDurations = [adhdFixationDurations, fixationStats(participant).conditions(condition).avgFixationDuration];
            adhdSaccadeAmplitudes = [adhdSaccadeAmplitudes, saccadeStats(participant).conditions(condition).avgSaccadeAmplitude];
            adhdSaccadeCounts = [adhdSaccadeCounts, saccadeStats(participant).conditions(condition).avgSaccadeCount];
        end
    else
        % Aggregate non-ADHD data
        for condition = 1:length(conditions)
            nonAdhdFixationDurations = [nonAdhdFixationDurations, fixationStats(participant).conditions(condition).avgFixationDuration];
            nonAdhdSaccadeAmplitudes = [nonAdhdSaccadeAmplitudes, saccadeStats(participant).conditions(condition).avgSaccadeAmplitude];
            nonAdhdSaccadeCounts = [nonAdhdSaccadeCounts, saccadeStats(participant).conditions(condition).avgSaccadeCount];
        end
    end
end

% Combine data into matrices for analysis
adhdData = [adhdFixationDurations; adhdSaccadeAmplitudes; adhdSaccadeCounts];
nonAdhdData = [nonAdhdFixationDurations; nonAdhdSaccadeAmplitudes; nonAdhdSaccadeCounts];

% Perform a t-test for each metric
[hFix, pFix] = ttest2(adhdData(1, :), nonAdhdData(1, :));  % Fixation Duration
[hAmp, pAmp] = ttest2(adhdData(2, :), nonAdhdData(2, :));  % Saccade Amplitude
[hCount, pCount] = ttest2(adhdData(3, :), nonAdhdData(3, :));  % Saccade Count

% Display the results
fprintf('Fixation Duration: p = %.3f\n', pFix);
fprintf('Saccade Amplitude: p = %.3f\n', pAmp);
fprintf('Saccade Count: p = %.3f\n', pCount);


figure;
subplot(1, 3, 1);
bar([mean(adhdData(1, :)), mean(nonAdhdData(1, :))]);
title('Fixation Duration');
xticklabels({'ADHD', 'nonADHD'});
ylabel('Duration (ms)');
grid on;

subplot(1, 3, 2);
bar([mean(adhdData(2, :)), mean(nonAdhdData(2, :))]);
title('Saccade Amplitude');
xticklabels({'ADHD', 'nonADHD'});
ylabel('Amplitude (pixels)');
grid on;

subplot(1, 3, 3);
bar([mean(adhdData(3, :)), mean(nonAdhdData(3, :))]);
title('Saccade Count');
xticklabels({'ADHD', 'nonADHD'});
ylabel('Count');
grid on;
