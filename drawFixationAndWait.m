function [samp_start, samp_fix, fixStartTime] = drawFixationAndWait(window, xCenter, yCenter, crossColour, crossLineWidthPix, crossSizePix, eye_tracker, screenXpixels, screenYpixels, fixationRadius, minFixationTime)
    % draw the fixation cross
    fixStartTime = drawFixation(window, xCenter, yCenter, crossColour, crossLineWidthPix, crossSizePix);
    
    %%% 21.5.2025 wait until fixation cross was fixated for 0.5 seconds to continue
    % Define the area around the fixation point where gaze must fall
    fixationRadius = 400; % pixels, defines a circle radius around the fixation point
    minFixationTime = 0.5; % seconds, the minimum time gaze must remain within the radius
    % samp_fix = drawFixationAndWait(window, xCenter, yCenter, color, lineWidthPix, crossSizePix, eye_tracker, screenXpixels, screenYpixels, fixationRadius, minFixationTime)
        % draw the fixation cross
        %fixStartTime = drawFixation(window, xCenter, yCenter, color, lineWidthPix, crossSizePix);
        
        % get gaze until now/from up to last time we consumed it
        samp_start = eye_tracker.buffer.consumeN('gaze');
        %samp_all(trial).samp = samp_trial;

        % struct for all the small samps during fixation in this trial
        samp_fix = struct();
        count = 1;  % Initialize the counter for indexing the struct array
        fixCountTime = fixStartTime;    
        while true
            % 0. restart gaze_tracking
            eye_tracker.buffer.start('gaze');
            
            % 1. track gaze shortly
            WaitSecs(0.01);
            
            % 2. fetch latest gaze data and analyze it
            samp = eye_tracker.buffer.consumeN('gaze'); 
            
            % ==== 🔴 CHECK FOR ESCAPE KEY ====
            [keyIsDown, ~, keyCode] = KbCheck;
            if keyIsDown
                if keyCode(KbName('ESCAPE'))
                    disp('Escape key pressed. Exiting gaze loop.');
                    break;
                end
            end
    
            % average gaze coordinates if both eyes are tracked, ignoring NaNs
            if ~isempty(samp) % Process if there's data
                samp_fix(count).samp = samp;  % store each samp in a struct array
                count = count + 1;  % Increment the index
                gazeX = mean([samp.left.gazePoint.onDisplayArea(1,:); samp.right.gazePoint.onDisplayArea(1,:)], 'omitnan');
                gazeY = mean([samp.left.gazePoint.onDisplayArea(2,:); samp.right.gazePoint.onDisplayArea(2,:)], 'omitnan');
                gazeX = gazeX * screenXpixels; % Scale to pixel space
                gazeY = gazeY * screenYpixels; % Scale to pixel space

                % check if gaze is within the fixation radius
                if sqrt((gazeX - xCenter).^2 + (gazeY - yCenter).^2) <= fixationRadius
                    if GetSecs - fixCountTime >= minFixationTime
                        break; % BREAK the loop if fixation is maintained long enough
                    end
                else % reset the timer if gaze wanders off
                    fixCountTime = GetSecs; 
                end % for debugging:
%                 disp(['Gaze Coordinates: (' num2str(gazeX) ', ' num2str(gazeY) ')']);  % Displays the gaze coordinates
%                 disp(['Distance from Center: ' num2str(sqrt((gazeX - xCenter).^2 + (gazeY - yCenter).^2))]);  % Distance from center
%                 disp(['Time Fixated: ' num2str(GetSecs - fixCountTime)]);  % How long the gaze has been within the radius

            end
            % 3. small delay to prevent high CPU usage
            WaitSecs(0.01); 
        end