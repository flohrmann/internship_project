function saccades = detectSaccadesFromFixations(fixations, stim_start, screenXpixels, screenYpixels)
    % Detects saccades based on the centers of consecutive fixations
    % fixations: detected fixations for a given trial
    % screenXpixels, screenYpixels: screen resolution in pixels

    saccades = struct('startIdx', {}, 'endIdx', {}, 'startCenter', {}, 'endCenter', {}, ...
                      'distance', {}, 'duration', {}, 'velocity', {}, ...
                      'saccStartAfterStimOnset', {},'saccEndAfterStimOnset', {});
    
    for i = 1:length(fixations) - 1
        % Get the centers of the current fixation and the next fixation
        startCenter = fixations(i).center;
        endCenter = fixations(i+1).center;

        % Calculate the saccade distance as the Euclidean distance between fixation centers
        saccadeDistance = sqrt((endCenter(1) - startCenter(1)).^2 + (endCenter(2) - startCenter(2)).^2);
        
        % Calculate the saccade duration as the time between the end of the current fixation
        % and the start of the next fixation (in terms of index points)
        saccadeDuration = fixations(i+1).startIdx - fixations(i).endIdx;
                
        % Calculate the saccade velocity (distance/time)
        if saccadeDuration > 0
            saccadeVelocity = saccadeDistance / (saccadeDuration * (1/60)); % datapoints into ms
        else
            saccadeVelocity = NaN;  % Handle division by zero
        end
        
        sacc_stim_start = fixations(i).fixEnd - stim_start;
        sacc_stim_end   = fixations(i+1).fixStart - stim_start;
        
        % Store saccade details
        saccades(end+1) = struct('startIdx', fixations(i).endIdx, 'endIdx', fixations(i+1).startIdx, ...
                                 'startCenter', startCenter, 'endCenter', endCenter, ...
                                 'distance', saccadeDistance, 'duration', saccadeDuration, ...
                                 'velocity', saccadeVelocity, ...
                                 'saccStartAfterStimOnset', sacc_stim_start ,'saccEndAfterStimOnset',sacc_stim_end);
    end
end
