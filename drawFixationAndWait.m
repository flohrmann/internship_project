function samp_fix = drawFixationAndWait(window, xCenter, yCenter, color, dotSizePix, eye_tracker, screenXpixels, screenYpixels, fixationRadius, minFixationTime)
    % Draw the fixation dot
    Screen('DrawDots', window, [xCenter; yCenter], dotSizePix, color, [], 2);
    Screen('Flip', window);
    
    fixStartTime = GetSecs;
    samp_fix = struct();
    count = 1;  % Initialize the counter for indexing the struct array

    while true
        eye_tracker.buffer.start('gaze'); % Restart gaze_tracking
        WaitSecs(0.01); % Track gaze shortly
        samp = eye_tracker.buffer.consumeN('gaze'); % Fetch latest gaze data
        % Average gaze coordinates if both eyes are tracked, ignoring NaNs
        if ~isempty(samp) % Process if there's data
            samp_fix(count).samp = samp;  % Store each samp in a struct array
            count = count + 1;  % Increment the index
            gazeX = mean([samp.left.gazePoint.onDisplayArea(1,:); samp.right.gazePoint.onDisplayArea(1,:)], 'omitnan');
            gazeY = mean([samp.left.gazePoint.onDisplayArea(2,:); samp.right.gazePoint.onDisplayArea(2,:)], 'omitnan');
            gazeX = gazeX * screenXpixels; % Scale to pixel space
            gazeY = gazeY * screenYpixels; % Scale to pixel space

            % Check if gaze is within the fixation radius
            if sqrt((gazeX - xCenter).^2 + (gazeY - yCenter).^2) <= fixationRadius
                if GetSecs - fixStartTime >= minFixationTime
                    break; % Break the loop if fixation is maintained long enough
                end
            else
                fixStartTime = GetSecs; % Reset the timer if gaze wanders off
            end % for debugging:
%             disp(['Gaze Coordinates: (' num2str(gazeX) ', ' num2str(gazeY) ')']);  % Displays the gaze coordinates
%             disp(['Distance from Center: ' num2str(sqrt((gazeX - xCenter).^2 + (gazeY - yCenter).^2))]);  % Distance from center
%             disp(['Time Fixated: ' num2str(GetSecs - fixStartTime)]);  % How long the gaze has been within the radius

        end
        WaitSecs(0.01); % Small delay to prevent high CPU usage
    end
end
