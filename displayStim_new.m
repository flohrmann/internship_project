function trial = displayStim_new(window, bar_wh_ratio, jitter_x, jitter_y, current_stim, current_cond, current_target_pos, screenXpixels, screenYpixels, color_bg, color_stim)
% Define the keyboard keys that are listened for. We will be using the left
% and right arrow keys as response keys for the task and the escape key as
% a exit/reset key
escapeKey = KbName('ESCAPE');
leftKey = KbName('LeftArrow');
rightKey = KbName('RightArrow');
% Hide the mouse cursor
%HideCursor; % disabled for testing

num_cols = size(current_stim, 2);
num_rows = size(current_stim, 1);

cellWidth = screenXpixels / num_cols;
cellHeight = screenYpixels / num_rows;
% Define a line length based on the cell size
lineLength = min(cellWidth, cellHeight) * 0.5;
lineWidth = 0.1 * lineLength;

trial = [];

for row = 1:num_rows %y
    for col = 1:num_cols %x
        
        % Calculate center of each cell
        xCenter = (col - 0.5) * cellWidth;
        yCenter = (row - 0.5) * cellHeight;

        % Retrieve angles for this cell
        angles = current_stim{row, col};

        % Loop over each angle and draw the lines
        for angle = angles
            dx = (lineLength / 2) * cosd(angle);
            dy = (lineLength / 2) * sind(angle);
            lineCoords = [-dx, dy; dx, -dy]'; % Column vector for x and y coordinates

            % Draw the line
            Screen('DrawLines', window, lineCoords, lineWidth, color_stim, [xCenter, yCenter], 2);
        end

        % Calculate base position of the bar
%         xPos = (col - 0.5) * (screenXpixels / num_rows);
%         yPos = (row - 0.5) * (screenYpixels / num_cols);
%         % Add random jitter to the position
%         jitterX = randi([-jitter_x, jitter_x]); % Adjusted variable names
%         jitterY = randi([-jitter_y, jitter_y]);
%         xPosJittered = xPos + jitterX; % Use jittered position for drawing
%         yPosJittered = yPos + jitterY;
% 
%         % Fetch the angles from the cell array
%         angles = current_stim{row, col};
%         if ~iscell(angles) % If the stored value is not a cell, make it a cell for uniform processing
%             angles = {angles};
%         end
% 
%         % Loop over each angle in the cell to draw one or more bars
%         for k = 1:length(angles)
%             angleRad = deg2rad(angles{k}); % Convert angle to radians
%             xEnd = xPosJittered + bar_height * cos(angleRad);
%             yEnd = yPosJittered - bar_height * sin(angleRad); % Subtract because Y-axis is inverted
% 
%             % Draw the line to the screen
%             Screen('DrawLine', window, color, xPosJittered, yPosJittered, xEnd, yEnd, bar_width);
%         end
    end
end

% Flip to the screen for each trial
Screen('Flip', window);

% Record the start time of the trial
trial_start_time = GetSecs;

%%% Code from https://peterscarfe.com/poserCuingExperiment.html :
% Now we wait for a keyboard button signaling the observers response.
% The left arrow key signals a "left" response and the right arrow key
% a "right" response. You can also press escape if you want to exit the
% program
respToBeMade = true;
startResp = GetSecs;
while respToBeMade
    [keyIsDown,secs, keyCode] = KbCheck(-1);
    if keyCode(escapeKey)
        ShowCursor;
        sca;
        return
    elseif keyCode(leftKey)
        response = 'L';
        respToBeMade = false;
    elseif keyCode(rightKey)
        response = 'R';
        respToBeMade = false;
    end
end
endResp = GetSecs;
rt = endResp - startResp;

% Work out if the location of the target was identified corrcetly
if current_target_pos == response
    correctness = 1;
elseif current_target_pos ~= response
    correctness = 0;
end
%%%%%%

% Record user response and response time
trial.condition = current_cond;
trial.target = current_target_pos;
trial.response = response;
trial.correct = correctness;
trial.start = startResp;
trial.end = endResp;
trial.rt = rt;

end