function trialSaccadeTimes = getSaccadeTimes(fixations, participant_data)
numFixations = length(fixations);
trialSaccadeTimes = [];

% Loop over each fixation (starting from the second one to calculate time between fixations)
for fIdx = 2:numFixations
    % Get the end index of the previous fixation and the start index of the current fixation
    prevEndIdx = fixations(fIdx-1).endIdx;
    currStartIdx = fixations(fIdx).startIdx;
    
    % Retrieve the corresponding timestamps from the eye-tracking data
    prevEndTime = double(participant_data.systemTimeStamp(prevEndIdx)) / 1e6;
    currStartTime = double(participant_data.systemTimeStamp(currStartIdx)) / 1e6;
    
    % Calculate the time between fixations (potential saccade duration)
    saccadeTime = currStartTime - prevEndTime;
    
    trialSaccadeTimes = [trialSaccadeTimes, saccadeTime];
end

