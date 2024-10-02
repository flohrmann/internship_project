function avg_stim_30 = calcAvgStimFound30DataPoints(cutData, conditions, num_points)

avg_stim_30 = table();

for trial = 1:size(cutData, 1)
    condition = cutData.Condition{trial}; 
    current_data = cutData(trial, :);
    trialEndTime = current_data.trialEndTime;
    
    % Extract system timestamps and pupil diameters for both eyes
    timestamps = double(current_data.eyeTrial.systemTimeStamp) / 1e6;
    pupil_right = current_data.eyeTrial.right.pupil.diameter;
    pupil_left = current_data.eyeTrial.left.pupil.diameter;
    
    % Average the two eyes
    pupil_avg = nanmean([pupil_right, pupil_left], 1);
    
    % Handle NaN by taking the mean of the previous and next values
    for i = 2:length(pupil_avg)-1
        if isnan(pupil_avg(i))
            pupil_avg(i) = nanmean([pupil_avg(i-1), pupil_avg(i+1)]);
        end
    end
    
    % Get the index of the trialEndTime
    [~, idx_end] = min(abs(timestamps - trialEndTime));
    
    % Extract the 30 points before trialEndTime
    if idx_end > num_points
        % Get the last 30 points before trialEndTime
        points_before_end = pupil_avg(idx_end - num_points + 1:idx_end);
    else
        % If less than 30 points before trialEndTime, pad with NaNs
        points_before_end = [nan(num_points - idx_end, 1); pupil_avg(1:idx_end)];
    end
    
    % Create a new row to append to the table
    new_row = {condition, points_before_end'};  % Store condition and data points as a row
    
    % Append the new row to the table
    avg_stim_30 = [avg_stim_30; new_row];
end

% Set the table's variable names
avg_stim_30.Properties.VariableNames = {'Condition', 'DataPoints'};

% Save the table as a .mat file
save('avg_stim_30_table.mat', 'avg_stim_30');
end