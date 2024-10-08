function [norm_diam_around_stim_adhd, norm_diam_around_stim_nonadhd] = getPupilDiamsPerGroup(data, avg, group_labels, conditions)

norm_diam_around_stim_adhd = struct();
norm_diam_around_stim_nonadhd = struct();

% Loop through each condition
for j = 1:length(conditions)
    condition = conditions{j};
    
    % Initialize arrays for storing per-participant averages
    diam_around_stim_adhd_before = [];
    diam_around_stim_adhd_after = [];
    diam_around_stim_nonadhd_before = [];
    diam_around_stim_nonadhd_after = [];
    
    median_baseline_adhd = [];
    mean_baseline_adhd = [];
    median_baseline_nonadhd = [];
    mean_baseline_nonadhd = [];

    % Loop through each participant
    for i = 1:length(data)
        participant = data(i);
        
        % Check if the participant belongs to ADHD or Non-ADHD group
        if strcmp(group_labels{i}, 'ADHD')
            % Extract trials with the current condition
            cond_mask = strcmp(participant.pupilDiam.Condition, condition);
            if sum(cond_mask) > 0
                % Compute the average pupil diameter for this participant across trials with this condition
                if strcmp(avg, 'mean')
                    avg_before = mean(cell2mat(participant.pupilDiam.MeanNormDataPointsBefore(cond_mask, :)), 1, 'omitnan');
                    avg_after = mean(cell2mat(participant.pupilDiam.MeanNormDataPointsAfter(cond_mask, :)), 1, 'omitnan');
                else % median
                    avg_before = median(cell2mat(participant.pupilDiam.MedianNormDataPointsBefore(cond_mask, :)), 1, 'omitnan');
                    avg_after = median(cell2mat(participant.pupilDiam.MedianNormDataPointsAfter(cond_mask, :)), 1, 'omitnan');
                end
                % Append to ADHD group data arrays
                diam_around_stim_adhd_before = [diam_around_stim_adhd_before; avg_before];
                diam_around_stim_adhd_after = [diam_around_stim_adhd_after; avg_after];
                median_baseline_adhd = [median_baseline_adhd; participant.pupilDiam.MedianDiamBlank];
                mean_baseline_adhd = [mean_baseline_adhd; participant.pupilDiam.MeanDiamBlank];
            end
            
        elseif strcmp(group_labels{i}, 'nonADHD')
            % Extract trials with the current condition
            cond_mask = strcmp(participant.pupilDiam.Condition, condition);
            if sum(cond_mask) > 0
                % Compute the average pupil diameter for this participant across trials with this condition
                if strcmp(avg, 'mean')
                    avg_before = mean(cell2mat(participant.pupilDiam.MeanNormDataPointsBefore(cond_mask, :)), 1, 'omitnan');
                    avg_after = mean(cell2mat(participant.pupilDiam.MeanNormDataPointsAfter(cond_mask, :)), 1, 'omitnan');
                else % median
                    avg_before = median(cell2mat(participant.pupilDiam.MedianNormDataPointsBefore(cond_mask, :)), 1, 'omitnan');
                    avg_after = median(cell2mat(participant.pupilDiam.MedianNormDataPointsAfter(cond_mask, :)), 1, 'omitnan');
                end
                % Append to Non-ADHD group data arrays
                diam_around_stim_nonadhd_before = [diam_around_stim_nonadhd_before; avg_before];
                diam_around_stim_nonadhd_after = [diam_around_stim_nonadhd_after; avg_after];
            end
        end
    end
    
    norm_diam_around_stim_adhd.(condition).before = diam_around_stim_adhd_before;
    norm_diam_around_stim_adhd.(condition).after = diam_around_stim_adhd_after;

    norm_diam_around_stim_nonadhd.(condition).before = diam_around_stim_nonadhd_before;
    norm_diam_around_stim_nonadhd.(condition).after = diam_around_stim_nonadhd_after;
end
