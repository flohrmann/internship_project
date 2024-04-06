function trial = displayStim(window, bar_width, bar_height, jitter, current_stim, current_cond, current_target_pos, screenXpixels, screenYpixels, color)
% Define the keyboard keys that are listened for. We will be using the left
% and right arrow keys as response keys for the task and the escape key as
% a exit/reset key
escapeKey = KbName('ESCAPE');
leftKey = KbName('LeftArrow');
rightKey = KbName('RightArrow');

trial = [];

for row = 1:size(current_stim, 1)
    for col = 1:size(current_stim, 2)
        % Calculate position of the bar
        xPos = (col - 0.5) * (screenXpixels / size(current_stim, 2));
        yPos = (row - 0.5) * (screenYpixels / size(current_stim, 1));

        % Add random jitter to the position
        xPos = xPos + randi([-jitter, jitter]);
        yPos = yPos + randi([-jitter, jitter]);

        % Calculate end point of the line based on angle
        angleRad = deg2rad(current_stim(row, col));
        xEnd = xPos + bar_height * cos(angleRad);
        yEnd = yPos - bar_height * sin(angleRad); % Subtract because Y-axis is inverted

        % Draw the line to the screen
        Screen('DrawLine', window, color, xPos, yPos, xEnd, yEnd, 2);
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
trial.target = current_target_pos;
trial.response = response;
trial.correct = correctness;
trial.start = startResp;
trial.end = endResp;
trial.rt = rt;

end
