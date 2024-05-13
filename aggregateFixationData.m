function return_data = aggregateFixationData(data)
% gets results from experiment and cleans up fixation data

return_data = table();

for n = 1:length(data.sampFix)
    % Initialize a struct to store aggregated timestamps and gaze points
    % in loop to make sure it overwrites for each trial
    clear aggregatedData;
    aggregatedData = struct('deviceTimeStamp', int64([]), 'systemTimeStamp', int64([]), ...
    'left', struct('gazePoint', struct('onDisplayArea', double([]), 'inUserCoords', double([]), 'valid', logical([]), 'available', logical([])), ...
    'pupil', struct('diameter', double([]), 'valid', logical([]), 'available', logical([])), ...
    'gazeOrigin', struct('inUserCoords', double([]), 'inTrackBoxCoords', double([]), 'valid', logical([]), 'available', logical([])), ...
    'eyeOpenness', struct('diameter', double([]), 'valid', logical([]), 'available', logical([]))), ...
    'right', struct('gazePoint', struct('onDisplayArea', double([]), 'inUserCoords', double([]), 'valid', logical([]), 'available', logical([])),...
    'pupil', struct('diameter', double([]), 'valid', logical([]), 'available', logical([])), ...
    'gazeOrigin', struct('inUserCoords', double([]), 'inTrackBoxCoords', double([]), 'valid', logical([]), 'available', logical([])), ...
    'eyeOpenness', struct('diameter', double([]), 'valid', logical([]), 'available', logical([]))));

    current_samp = data.sampFix{n};
    % Loop through each samp in sampFix
    for i = 1:length(current_samp)
        % Append timestamps to the arrays in the struct
        sub_samp = current_samp(i).samp;
        aggregatedData.deviceTimeStamp = [aggregatedData.deviceTimeStamp, sub_samp.deviceTimeStamp];
        aggregatedData.systemTimeStamp = [aggregatedData.systemTimeStamp, sub_samp.systemTimeStamp];

        % Check the number of values in the current deviceTimeStamp
        %numValues = numel(sub_samp.deviceTimeStamp)

        % Append gaze points to the arrays in the struct
        % left eye
        aggregatedData.left.gazePoint.onDisplayArea = [aggregatedData.left.gazePoint.onDisplayArea, sub_samp.left.gazePoint.onDisplayArea];
        aggregatedData.left.gazePoint.inUserCoords = [aggregatedData.left.gazePoint.inUserCoords, sub_samp.left.gazePoint.inUserCoords];
        aggregatedData.left.gazePoint.valid = [aggregatedData.left.gazePoint.valid, sub_samp.left.gazePoint.valid];
        aggregatedData.left.gazePoint.available = [aggregatedData.left.gazePoint.available, sub_samp.left.gazePoint.available];

        aggregatedData.left.gazeOrigin.inUserCoords = [aggregatedData.left.gazeOrigin.inUserCoords, sub_samp.left.gazeOrigin.inUserCoords];
        aggregatedData.left.gazeOrigin.inTrackBoxCoords = [aggregatedData.left.gazeOrigin.inTrackBoxCoords, sub_samp.left.gazeOrigin.inTrackBoxCoords];
        aggregatedData.left.gazeOrigin.valid = [aggregatedData.left.gazeOrigin.valid, sub_samp.left.gazeOrigin.valid];
        aggregatedData.left.gazeOrigin.available = [aggregatedData.left.gazeOrigin.available, sub_samp.left.gazeOrigin.available];

        aggregatedData.left.eyeOpenness.diameter = [aggregatedData.left.eyeOpenness.diameter, sub_samp.left.eyeOpenness.diameter];
        aggregatedData.left.eyeOpenness.valid = [aggregatedData.left.eyeOpenness.valid, sub_samp.left.eyeOpenness.valid];
        aggregatedData.left.eyeOpenness.available = [aggregatedData.left.eyeOpenness.available, sub_samp.left.eyeOpenness.available];

        % right eye
        aggregatedData.right.gazePoint.onDisplayArea = [aggregatedData.right.gazePoint.onDisplayArea, sub_samp.right.gazePoint.onDisplayArea];
        aggregatedData.right.gazePoint.inUserCoords = [aggregatedData.right.gazePoint.inUserCoords, sub_samp.right.gazePoint.inUserCoords];
        aggregatedData.right.gazePoint.valid = [aggregatedData.right.gazePoint.valid, sub_samp.right.gazePoint.valid];
        aggregatedData.right.gazePoint.available = [aggregatedData.right.gazePoint.available, sub_samp.right.gazePoint.available];

        aggregatedData.right.gazeOrigin.inUserCoords = [aggregatedData.right.gazeOrigin.inUserCoords, sub_samp.right.gazeOrigin.inUserCoords];
        aggregatedData.right.gazeOrigin.inTrackBoxCoords = [aggregatedData.right.gazeOrigin.inTrackBoxCoords, sub_samp.right.gazeOrigin.inTrackBoxCoords];
        aggregatedData.right.gazeOrigin.valid = [aggregatedData.right.gazeOrigin.valid, sub_samp.right.gazeOrigin.valid];
        aggregatedData.right.gazeOrigin.available = [aggregatedData.right.gazeOrigin.available, sub_samp.right.gazeOrigin.available];

        aggregatedData.right.eyeOpenness.diameter = [aggregatedData.right.eyeOpenness.diameter, sub_samp.right.eyeOpenness.diameter];
        aggregatedData.right.eyeOpenness.valid = [aggregatedData.right.eyeOpenness.valid, sub_samp.right.eyeOpenness.valid];
        aggregatedData.right.eyeOpenness.available = [aggregatedData.right.eyeOpenness.available, sub_samp.right.eyeOpenness.available];
    end
    return_data(n,1) = {aggregatedData(1,:)};
    % Optionally display the new struct
    %disp(aggregatedData);
end
%return_data.Properties.VariableNames = {'FixSamp'};
