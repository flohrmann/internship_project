
function normalized_data = calcProportionalPupilChange(data)
% Initialize a new struct to store the normalized (percentage change) data
normalized_data = data;

% Loop through each participant
for i = 1:length(data)
    participant = data(i);
    
    
    % Loop through each trial and calculate the percentage change relative to the mean baseline
    for trial = 1:height(participant.pupilDiam)
        data_before_mean = participant.pupilDiam.DataPointsBefore(trial,:);
        data_after_mean = participant.pupilDiam.DataPointsAfter(trial,:);
        
        baseline_mean = participant.pupilDiam.MeanDiamBlank(trial);  % Mean baseline pupil diameter
        baseline_median = participant.pupilDiam.MedianDiamBlank(trial);  % Median baseline pupil diameter
        
        % (X - Baseline) / Baseline * 100
        proportional_change_before_mean = ((data_before_mean - baseline_mean) / baseline_mean) * 100;
        proportional_change_after_mean = ((data_after_mean - baseline_mean) / baseline_mean) * 100;
        
        % Store the results in the struct (for mean-based change)
        normalized_data(i).pupilDiam.MeanNormDataPointsBefore{trial} = proportional_change_before_mean;
        normalized_data(i).pupilDiam.MeanNormDataPointsAfter{trial} = proportional_change_after_mean;
        
        % Percentage change relative to the median baseline
        data_before_median = participant.pupilDiam.DataPointsBefore(trial,:);
        data_after_median = participant.pupilDiam.DataPointsAfter(trial,:);
        
        proportional_change_before_median = ((data_before_median - baseline_median) / baseline_median) * 100;
        proportional_change_after_median = ((data_after_median - baseline_median) / baseline_median) * 100;
        
        % Store the results in the struct (for median-based change)
        normalized_data(i).pupilDiam.MedianNormDataPointsBefore{trial} = proportional_change_before_median;
        normalized_data(i).pupilDiam.MedianNormDataPointsAfter{trial} = proportional_change_after_median;
    end
end
% Convert the cell arrays to regular double arrays
normalized_data.MeanNormDataPointsBefore = cell2mat(normalized_data.MeanNormDataPointsBefore);
normalized_data.MeanNormDataPointsAfter = cell2mat(normalized_data.MeanNormDataPointsAfter);
normalized_data.MedianNormDataPointsBefore = cell2mat(normalized_data.MedianNormDataPointsBefore);
normalized_data.MedianNormDataPointsAfter = cell2mat(normalized_data.MedianNormDataPointsAfter);

end