function fixations = analyseFixation(id, plot_these, data, dist_threshold, min_duration, screenXpixels, screenYpixels, safe, analysis_folder, length_fixation)

    % Initialize arrays for storing fixations and target centers
    targetCenters = [];
    eyeFixations = [];
    stimulusFixations = [];
    

    % Loop over each trial
    for trialIdx = 1:length(data.eyeTrial)
        
        % Calculate target center coordinates
        target = data.TargetPosition(trialIdx, :);
        target_x_center = data.x_centers{trialIdx}(target(2), target(1));
        target_y_center = screenYpixels - data.y_centers{trialIdx}(target(2), target(1));
        %target_y_center = data.y_centers{trialIdx}(target(2), target(1));

        % Store the target center for this trial
        targetCenters = [targetCenters; [target_x_center, target_y_center]];
        
        % Convert gaze and target positions for eyeTrial
        eyeFixation = detect_fixations(data.eyeTrial(trialIdx), screenXpixels, screenYpixels, dist_threshold, min_duration);
        
        % Convert gaze and target positions for stimulusTrial
        stimulusFixation = detect_fixations(data.stimulusTrial(trialIdx), screenXpixels, screenYpixels, dist_threshold, min_duration);
        
        % timepoint trial/stimulus start
        trial_start    = double(data.eyeTrial(trialIdx).systemTimeStamp(1,1))/ 1e6;
        stimulus_start = double(data.stimulusTrial(trialIdx).systemTimeStamp(1,1))/ 1e6;
        % Store the fixations for this trial
        eyeFixations = [eyeFixations, struct('trial', trialIdx, 'trialStart', trial_start ,'fixations', eyeFixation)];
        stimulusFixations = [stimulusFixations, struct('trial', trialIdx, 'stimStart', stimulus_start, 'fixations', stimulusFixation)];
    end

    % Store the fixation data along with the target information
    fixations.targetCenters = targetCenters'; % (x,y)
    fixations.eyeFixations = eyeFixations;
    fixations.stimulusFixations = stimulusFixations;
    fixations.length = length_fixation;
    fixations.id = id;
    try
        mkdir(analysis_folder);
    catch % folders already exists
    end
    save(strcat(analysis_folder, '\fixations.mat'), 'fixations');
    %save(fullfile(analysis_folder, 'fixations.mat'), 'fixations');

    % Display the results if safe is true
    if safe
        %disp('Fixations detected:');
        %display_fixations(fixations.eyeFixations);
        %display_fixations(fixations.stimulusFixations);
    end
    
    
    % plot some fixations for checking
    if length(plot_these) >= 1
        plotFixations(data, fixations, screenXpixels, screenYpixels, plot_these, analysis_folder, length_fixation, safe);
        close all
    else
    end
end


%% 
function [fixations] = detect_fixations(trialData, screenXpixels, screenYpixels, distThreshold, minDuration)
    % Detects fixations based on proximity of consecutive gaze points and includes target information.
    % trialData: structure containing trial data (either eyeTrial or stimulusTrial)
    % distThreshold: maximum distance between consecutive points to consider as fixation (in pixels).
    % minDuration: minimum number of consecutive points to qualify as a fixation.

    % Preallocate the fixations array
    fixations = struct('startIdx', {}, 'endIdx', {}, 'center', {}, 'duration', {}, 'fixStart', {}, 'fixEnd', {});

    % Convert gaze points to screen coordinates
    xGaze = (trialData.left.gazePoint.onDisplayArea(1, :)) * screenXpixels;
    yGaze = screenYpixels - (trialData.left.gazePoint.onDisplayArea(2, :)) * screenYpixels;

    % Calculate the Euclidean distance between consecutive gaze points
    distances = sqrt(diff(xGaze).^2 + diff(yGaze).^2);

    % Identify the start and end of fixation periods
    isFixation = distances <= distThreshold;
    fixationGroups = split_fixations(isFixation, minDuration);

    % Store the fixations
    for f = 1:length(fixationGroups)
        startIdx = fixationGroups{f}(1);
        endIdx = fixationGroups{f}(end);
        duration = endIdx - startIdx + 1;
        centerX = mean(xGaze(startIdx:endIdx));
        centerY = mean(yGaze(startIdx:endIdx));
        fix_start = double(trialData.systemTimeStamp(1,startIdx))/ 1e6;
        fix_end = double(trialData.systemTimeStamp(1,endIdx))/ 1e6;
        %fixations(end+1) = struct('startIdx', startIdx, 'endIdx', endIdx, 'center', [centerY,centerX], 'duration', duration, 'fixStart', fix_start, 'fixEnd', fix_end);
        fixations(end+1) = struct('startIdx', startIdx, 'endIdx', endIdx, 'center', [centerX, centerY], 'duration', duration, 'fixStart', fix_start, 'fixEnd', fix_end);
    end
end


%% 
function display_fixations(fixations)
    % Display the fixations
    for i = 1:length(fixations)
        for j = 1:length(fixations(i).fixations)
            disp(['Trial: ', num2str(fixations(i).trial), ...
                  ', Start Index: ', num2str(fixations(i).fixations(j).startIdx), ...
                  ', End Index: ', num2str(fixations(i).fixations(j).endIdx), ...
                  ', Duration: ', num2str(fixations(i).fixations(j).duration)]);
        end
    end
end


%%
function [fixationGroups] = split_fixations(isFixation, minDuration)
    % This helper function groups consecutive indices of fixations
    fixationIdx = find(isFixation);
    d = diff(fixationIdx);
    splitPoints = [0 find(d > 1) length(fixationIdx)];
    fixationGroups = {};
    
    for j = 1:length(splitPoints)-1
        group = fixationIdx(splitPoints(j)+1:splitPoints(j+1));
        if length(group) >= minDuration
            fixationGroups{end+1} = group;
        end
    end
end

function [fixationsByTrial] = reorganize_fixations_by_trial(fixations)
    % This function reorganizes the fixations data so that each trial has its own struct.

    % Get the unique trial numbers
    uniqueTrials = unique([fixations.trial]);

    % Initialize the output struct array
    fixationsByTrial = struct('trial', [], 'fixations', []);
    fixationsByTrial = repmat(fixationsByTrial, 1, length(uniqueTrials));

    % Loop through each trial
    for i = 1:length(uniqueTrials)
        trialNumber = uniqueTrials(i);
        trialFixations = fixations([fixations.trial] == trialNumber);
        fixationsByTrial(i).trial = trialNumber;
        fixationsByTrial(i).fixations = trialFixations;
    end
end


