function fixStartTime = drawFixation(window, xCenter, yCenter, color, lineWidthPix, crossSizePix)
% Set the default color to red if not provided
if nargin < 4
    color = [255 0 0]; % Red color
end

% Set the default size of the cross if not provided
if nargin < 5
    lineWidthPix = 10; % Line width of the cross
end

% Set the default line width of the cross if not provided
if nargin < 6
    crossSizePix = 40; % Size of the cross
end

% Coordinates for the horizontal and vertical rectangles of the cross
horizontalRect = [xCenter - crossSizePix / 2, yCenter - lineWidthPix / 2, xCenter + crossSizePix / 2, yCenter + lineWidthPix / 2];
verticalRect = [xCenter - lineWidthPix / 2, yCenter - crossSizePix / 2, xCenter + lineWidthPix / 2, yCenter + crossSizePix / 2];

% Draw the horizontal part of the cross
Screen('FillRect', window, color, horizontalRect);

% Draw the vertical part of the cross
Screen('FillRect', window, color, verticalRect);

% Flip the screen to show the fixation cross
fixStartTime = Screen('Flip', window);

% Screen('DrawDots', window, [xCenter; yCenter], dotSizePix, color, [], 2);
% Screen('Flip', window);

% end with keystroke
% KbStrokeWait(-1);