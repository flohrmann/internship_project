% Clear the workspace and the screen
close all;
clear;
sca;

% test data
data_struct = load('C:\Users\idm\Desktop\Semester4\Internship\Matlab\Data\test_1_20240421_151825\rand_trials.mat');
data = data_struct.rand_trials;

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% Get the screen numbers
screens = Screen('Screens');

% Select the external screen if it is present, else revert to the native screen
screenNumber = max(screens);

% Define black, white and grey
black = BlackIndex(screenNumber);
white = WhiteIndex(screenNumber);
grey = white / 2;

% Open an on screen window and color it grey
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey);

% Set the blend function for the screen
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Get the size of the on screen window in pixels
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Define grid size (assuming a 9x12 grid from your example)
numRows = 9;
numCols = 12;
cellWidth = screenXpixels / numCols;
cellHeight = screenYpixels / numRows;

% Define a line length based on the cell size
lineLength = min(cellWidth, cellHeight) * 0.5;
lineWidth = 0.1 * lineLength; % Line width

anglesMatrix = data.AngleMatrix{1};

% Loop through each cell and draw lines based on angles
for row = 1:numRows
    for col = 1:numCols
        % Calculate center of each cell
        xCenter = (col - 0.5) * cellWidth;
        yCenter = (row - 0.5) * cellHeight;

        % Retrieve angles for this cell
        angles = anglesMatrix{row, col};

        % Loop over each angle and draw the lines
        for angle = angles
            dx = (lineLength / 2) * cosd(angle);
            dy = (lineLength / 2) * sind(angle);
            lineCoords = [-dx, dy; dx, -dy]'; % Column vector for x and y coordinates

            % Draw the line
            Screen('DrawLines', window, lineCoords, lineWidth, white, [xCenter, yCenter], 2);
        end
    end
end

% Flip to the screen
Screen('Flip', window);

% Wait for a key press
KbStrokeWait;

% Clear the screen
sca;
