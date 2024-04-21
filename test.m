% Setup Psychtoolbox
PsychDefaultSetup(2);
screens = Screen('Screens');
screenNumber = max(screens);

% Open a full screen window
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, 1);

% Read the image from a file
img = imread('C:\Users\idm\Desktop\Semester4\Internship\Matlab\ExperimentRepo\Images\example_1.png');

% Create a texture from the image
texture = Screen('MakeTexture', window, img);

% Get the size of the screen
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Calculate the new dimensions for the image to be 70% of the screen size
newWidth = screenXpixels * 0.7;
newHeight = screenYpixels * 0.7;

% Calculate the position of the scaled image to center it
dstX = (screenXpixels - newWidth) / 2;
dstY = (screenYpixels - newHeight) / 2;
dstRect = [dstX dstY dstX + newWidth dstY + newHeight];

% Draw the scaled texture to the window
Screen('DrawTexture', window, texture, [], dstRect);

% Add text on top of the image
Screen('TextSize', window, 30);
DrawFormattedText(window, 'Check this out!', 'center', dstY * 0.9, 0);

% Update the display
Screen('Flip', window);

% Pause to view the result
WaitSecs(3);

% Clean up
Screen('Close', texture);
sca;
