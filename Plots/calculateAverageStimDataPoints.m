function stim_data_points = calculateAverageStimDataPoints(cutData, conditions, analysis_folder)
% count number of datapoints in trials per condition for one participant
% returns means and prints them to cmd

stim_data_points = struct();
disp('Average number of data points in stimulus segment per condition:');

% Loop through each condition
for condIdx = 1:length(conditions)
    condition = conditions{condIdx};
    data_point_counts = [];
    
    % Loop through each trial in cutData and select based on condition
    for trial = 1:size(cutData, 1)
        if strcmp(cutData.Condition{trial}, condition) % Check condition match
            % Extract the current trial data
            current_data = cutData(trial, :);
            
            % Extract the relevant times
            StimulusOnsetTime = current_data.StimulusOnsetTime;
            trialEndTime = current_data.trialEndTime;
            
            % Extract system timestamps
            timestamps = double(current_data.eyeTrial.systemTimeStamp) / 1e6;
            
            % Get the indices corresponding to the key time points
            [~, idx_stimulus] = min(abs(timestamps - StimulusOnsetTime));
            [~, idx_end] = min(abs(timestamps - trialEndTime));
            
            % Count the number of data points in the stimulus segment
            num_stim_points = idx_end - idx_stimulus + 1;
            data_point_counts = [data_point_counts, num_stim_points];
        end
    end
    
    % Calculate the average number of data points for the current condition
    if ~isempty(data_point_counts)
        stim_data_points.condition     = condition;
        stim_data_points.avgDataPoints = mean(data_point_counts);
    else
        stim_data_points.condition     = condition;
        stim_data_points.avgDataPoints = 0;
        
    end
    % Display the average data points per condition
    disp(stim_data_points);
end

    save(fullfile(analysis_folder, 'avg_num_stim_data_points_per_condition.mat'), 'stim_data_points');

end
