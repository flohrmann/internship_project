function plotGroupFixationDurations(fixationStats, ids, group_labels, labels, conditions, participant_map, color_map, comparison_results_folder, fixation_duration, safe)
% Plot fixation durations for ADHD and nonADHD groups with individual participant data.
% Arguments:
%   fixationStats: Struct array with fixation data for each participant.
%   group_labels: Cell array of group labels ('nonADHD' or 'ADHD').
%   labels: Cell array of condition names.
%   conditions: Array of condition names.
%   participant_map: Map containing color and marker data for participants.
%   color_map: Map containing predefined colors for conditions and groups.
%   comparison_results_folder: Folder to save the plots.

% Initialize variables to store fixation data by group and condition
numConditions = length(conditions);
adhd_medians = zeros(0, numConditions); % Store data for ADHD group
nonadhd_medians = zeros(0, numConditions); % Store data for non-ADHD group

% Loop over participants to group data by ADHD status
for participant = 1:length(fixationStats)
    if strcmp(group_labels{participant}, 'ADHD')
        adhd_medians = [adhd_medians; [fixationStats(participant).a_median, fixationStats(participant).b_median, ...
                                       fixationStats(participant).as_median,fixationStats(participant).bs_median]];
    else
        nonadhd_medians = [nonadhd_medians; [fixationStats(participant).a_median, fixationStats(participant).b_median, ...
                                             fixationStats(participant).as_median,fixationStats(participant).bs_median]];
    end
end

% Calculate means and standard deviations for each group
median_adhdMedians = median(adhd_medians, 1, 'omitnan');
sem_adhd = std(adhd_medians, 0, 1, 'omitnan')./ sqrt(size(adhd_medians, 1));
median_nonAdhdMedians = median(nonadhd_medians, 1, 'omitnan');
sem_nonadhd = std(nonadhd_medians, 0, 1, 'omitnan')./ sqrt(size(nonadhd_medians, 1));

plotADHDnonADHDandDiff('Average Fixation Duration',... % sgtitle
                        adhd_medians, median_adhdMedians, sem_adhd, 'ADHD', 'northwest', ...       % adhd data
                        nonadhd_medians, median_nonAdhdMedians, sem_nonadhd, 'nonADHD', 'northwest',...  % nonadhd data
                        ids, labels, 'Median Fixation Duration (s)', ... % x, y axis labels
                        group_labels, conditions, color_map, participant_map, ...
                        fullfile(comparison_results_folder, '01_fixation_duration_condition_group_diff.png'));
end
