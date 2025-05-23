function [trialwise_gaze, continuous_gaze] = aggregateFixationData(data)
% gets results from experiment and cleans up fixation data

% eyetracking data per trial
trialwise_gaze = table();

% continous eyetracking data
continuous_gaze = struct('deviceTimeStamp', int64([]), 'systemTimeStamp', int64([]), ...
        'left', struct('gazePoint', struct('onDisplayArea', double([]), 'inUserCoords', double([]), 'valid', logical([]), 'available', logical([])), ...
        'pupil', struct('diameter', double([]), 'valid', logical([]), 'available', logical([])), ...
        'gazeOrigin', struct('inUserCoords', double([]), 'inTrackBoxCoords', double([]), 'valid', logical([]), 'available', logical([])), ...
        'eyeOpenness', struct('diameter', double([]), 'valid', logical([]), 'available', logical([]))), ...
        'right', struct('gazePoint', struct('onDisplayArea', double([]), 'inUserCoords', double([]), 'valid', logical([]), 'available', logical([])),...
        'pupil', struct('diameter', double([]), 'valid', logical([]), 'available', logical([])), ...
        'gazeOrigin', struct('inUserCoords', double([]), 'inTrackBoxCoords', double([]), 'valid', logical([]), 'available', logical([])), ...
        'eyeOpenness', struct('diameter', double([]), 'valid', logical([]), 'available', logical([]))));
    
for trial = 1:length(data)
    % Initialize a struct to store aggregated timestamps and gaze points
    % in loop to make sure it overwrites for each trial
    clear aggregatedData;
    aggregatedTrialData = struct('deviceTimeStamp', int64([]), 'systemTimeStamp', int64([]), ...
        'left', struct('gazePoint', struct('onDisplayArea', double([]), 'inUserCoords', double([]), 'valid', logical([]), 'available', logical([])), ...
        'pupil', struct('diameter', double([]), 'valid', logical([]), 'available', logical([])), ...
        'gazeOrigin', struct('inUserCoords', double([]), 'inTrackBoxCoords', double([]), 'valid', logical([]), 'available', logical([])), ...
        'eyeOpenness', struct('diameter', double([]), 'valid', logical([]), 'available', logical([]))), ...
        'right', struct('gazePoint', struct('onDisplayArea', double([]), 'inUserCoords', double([]), 'valid', logical([]), 'available', logical([])),...
        'pupil', struct('diameter', double([]), 'valid', logical([]), 'available', logical([])), ...
        'gazeOrigin', struct('inUserCoords', double([]), 'inTrackBoxCoords', double([]), 'valid', logical([]), 'available', logical([])), ...
        'eyeOpenness', struct('diameter', double([]), 'valid', logical([]), 'available', logical([]))));
    
    % start of trial 
    aggregatedTrialData = addSampToSampStruct(data(trial).sampStart, aggregatedTrialData);
    continuous_gaze = addSampToSampStruct(data(trial).sampStart, continuous_gaze);
    
    % during fixation
    samp_fix = data(trial).sampFix; %{trial};    
    for i = 1:length(samp_fix)
        % Append timestamps to the arrays in the struct
        samp_new = samp_fix(i).samp;
        aggregatedTrialData = addSampToSampStruct(samp_new, aggregatedTrialData);
        continuous_gaze = addSampToSampStruct(samp_new, continuous_gaze);
    end
    
    % stimulation and blank screen 
    aggregatedTrialData = addSampToSampStruct(data(trial).sampBlankStim, aggregatedTrialData);
    continuous_gaze = addSampToSampStruct(data(trial).sampBlankStim, continuous_gaze);
    
    trialwise_gaze(trial,1) = {aggregatedTrialData(1,:)};    
end
%return_data.Properties.VariableNames = {'FixSamp'};
% Optionally display the new struct
%disp(continuous_gaze);
end


function aggregatedData = addSampToSampStruct(samp_new, aggregatedData)

aggregatedData.deviceTimeStamp = [aggregatedData.deviceTimeStamp, samp_new.deviceTimeStamp];
aggregatedData.systemTimeStamp = [aggregatedData.systemTimeStamp, samp_new.systemTimeStamp];

% Check the number of values in the current deviceTimeStamp
%numValues = numel(sub_samp.deviceTimeStamp)

% Append gaze points to the arrays in the struct
% left eye
aggregatedData.left.gazePoint.onDisplayArea = [aggregatedData.left.gazePoint.onDisplayArea, samp_new.left.gazePoint.onDisplayArea];
aggregatedData.left.gazePoint.inUserCoords = [aggregatedData.left.gazePoint.inUserCoords, samp_new.left.gazePoint.inUserCoords];
aggregatedData.left.gazePoint.valid = [aggregatedData.left.gazePoint.valid, samp_new.left.gazePoint.valid];
aggregatedData.left.gazePoint.available = [aggregatedData.left.gazePoint.available, samp_new.left.gazePoint.available];

aggregatedData.left.gazeOrigin.inUserCoords = [aggregatedData.left.gazeOrigin.inUserCoords, samp_new.left.gazeOrigin.inUserCoords];
aggregatedData.left.gazeOrigin.inTrackBoxCoords = [aggregatedData.left.gazeOrigin.inTrackBoxCoords, samp_new.left.gazeOrigin.inTrackBoxCoords];
aggregatedData.left.gazeOrigin.valid = [aggregatedData.left.gazeOrigin.valid, samp_new.left.gazeOrigin.valid];
aggregatedData.left.gazeOrigin.available = [aggregatedData.left.gazeOrigin.available, samp_new.left.gazeOrigin.available];

aggregatedData.left.eyeOpenness.diameter = [aggregatedData.left.eyeOpenness.diameter, samp_new.left.eyeOpenness.diameter];
aggregatedData.left.eyeOpenness.valid = [aggregatedData.left.eyeOpenness.valid, samp_new.left.eyeOpenness.valid];
aggregatedData.left.eyeOpenness.available = [aggregatedData.left.eyeOpenness.available, samp_new.left.eyeOpenness.available];

% right eye
aggregatedData.right.gazePoint.onDisplayArea = [aggregatedData.right.gazePoint.onDisplayArea, samp_new.right.gazePoint.onDisplayArea];
aggregatedData.right.gazePoint.inUserCoords = [aggregatedData.right.gazePoint.inUserCoords, samp_new.right.gazePoint.inUserCoords];
aggregatedData.right.gazePoint.valid = [aggregatedData.right.gazePoint.valid, samp_new.right.gazePoint.valid];
aggregatedData.right.gazePoint.available = [aggregatedData.right.gazePoint.available, samp_new.right.gazePoint.available];

aggregatedData.right.gazeOrigin.inUserCoords = [aggregatedData.right.gazeOrigin.inUserCoords, samp_new.right.gazeOrigin.inUserCoords];
aggregatedData.right.gazeOrigin.inTrackBoxCoords = [aggregatedData.right.gazeOrigin.inTrackBoxCoords, samp_new.right.gazeOrigin.inTrackBoxCoords];
aggregatedData.right.gazeOrigin.valid = [aggregatedData.right.gazeOrigin.valid, samp_new.right.gazeOrigin.valid];
aggregatedData.right.gazeOrigin.available = [aggregatedData.right.gazeOrigin.available, samp_new.right.gazeOrigin.available];

aggregatedData.right.eyeOpenness.diameter = [aggregatedData.right.eyeOpenness.diameter, samp_new.right.eyeOpenness.diameter];
aggregatedData.right.eyeOpenness.valid = [aggregatedData.right.eyeOpenness.valid, samp_new.right.eyeOpenness.valid];
aggregatedData.right.eyeOpenness.available = [aggregatedData.right.eyeOpenness.available, samp_new.right.eyeOpenness.available];
end