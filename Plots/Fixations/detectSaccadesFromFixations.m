function saccades = detectSaccadesFromFixations(fixations, screenXpixels, screenYpixels)
    % Detects saccades based on the centers of consecutive fixations
    % fixations: detected fixations for a given trial
    % screenXpixels, screenYpixels: screen resolution in pixels

    saccades = struct('startIdx', {}, 'endIdx', {}, 'startCenter', {}, 'endCenter', {}, ...
                      'distance', {}, 'duration', {}, 'velocity', {});
    
    for i = 1:length(fixations) - 1
        % Get the centers of the current fixation and the next fixation
        startCenter = fixations(i).center;
        endCenter = fixations(i+1).center;

        % Calculate the saccade distance as the Euclidean distance between fixation centers
        saccadeDistance = sqrt((endCenter(1) - startCenter(1)).^2 + (endCenter(2) - startCenter(2)).^2);
        
        % Calculate the saccade duration as the time between the end of the current fixation
        % and the start of the next fixation (in terms of index points)
        saccadeDuration = fixations(i+1).startIdx - fixations(i).endIdx -1;
        
        % Calculate the saccade velocity (distance/time)
        if saccadeDuration > 0
            saccadeVelocity = saccadeDistance / (saccadeDuration * (1/60)); % datapoints into ms
        else
            saccadeVelocity = NaN;  % Handle division by zero
        end
        
        % Store saccade details
        saccades(end+1) = struct('startIdx', fixations(i).endIdx + 1, 'endIdx', fixations(i+1).startIdx - 1, ...
                                 'startCenter', startCenter, 'endCenter', endCenter, ...
                                 'distance', saccadeDistance, 'duration', saccadeDuration, ...
                                 'velocity', saccadeVelocity);
    end
end
