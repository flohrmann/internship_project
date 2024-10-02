function fixations = analyseFixation(data, distThreshold, minDuration, screenXpixels, screenYpixels, safe, analysis_folder)

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
        
        % Store the target center for this trial
        targetCenters = [targetCenters; [target_x_center, target_y_center]];
        
        % Convert gaze and target positions for eyeTrial
        eyeFixation = detect_fixations_with_target(data.eyeTrial(trialIdx), screenXpixels, screenYpixels, distThreshold, minDuration);
        
        % Convert gaze and target positions for stimulusTrial
        stimulusFixation = detect_fixations_with_target(data.stimulusTrial(trialIdx), screenXpixels, screenYpixels, distThreshold, minDuration);
        
        % Store the fixations for this trial
        eyeFixations = [eyeFixations, struct('trial', trialIdx, 'fixations', eyeFixation)];
        stimulusFixations = [stimulusFixations, struct('trial', trialIdx, 'fixations', stimulusFixation)];
    end

    % Store the fixation data along with the target information
    fixations.targetCenters = targetCenters';
    fixations.eyeFixations = eyeFixations;
    fixations.stimulusFixations = stimulusFixations;

    save(fullfile(analysis_folder, 'fixations.mat'), 'fixations');

    % Display the results if safe is true
    if safe
        disp('Fixations detected:');
        %display_fixations(fixations.eyeFixations);
        display_fixations(fixations.stimulusFixations);
    end
end


%% 
function [fixations] = detect_fixations_with_target(trialData, screenXpixels, screenYpixels, distThreshold, minDuration)
    % Detects fixations based on proximity of consecutive gaze points and includes target information.
    % trialData: structure containing trial data (either eyeTrial or stimulusTrial)
    % distThreshold: maximum distance between consecutive points to consider as fixation (in pixels).
    % minDuration: minimum number of consecutive points to qualify as a fixation.

    % Preallocate the fixations array
    fixations = struct('startIdx', {}, 'endIdx', {}, 'center', {}, 'duration', {});

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
        fixations(end+1) = struct('startIdx', startIdx, 'endIdx', endIdx, 'center', [centerX, centerY], 'duration', duration);
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


