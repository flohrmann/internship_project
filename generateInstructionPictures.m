% save two pictures of stims
 
%close all;
clear;
sca;

matrix = {
    [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0];
    [45 0], [45 90], [135 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90];
    [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0];
    [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90];
    [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0];
    [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90];
    [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0];
    [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90];
    [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0];
    };

simple_matrix = {
    [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45];
    [45], [45], [135], [45], [45], [45], [45], [45], [45], [45], [45], [45];
    [45], [45], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ];
    [45], [45], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ];
    [45], [45], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ];
    [45], [45], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ];
    [45], [45], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ];
    [45], [45], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ];
    [45], [45], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ], [45 ];
    };

PsychDefaultSetup(2);
screens = Screen('Screens');
[window, windowRect] = PsychImaging('OpenWindow', 0, 1);
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Define grid size
numRows = 9;
numCols = 12;
cellWidth = screenXpixels / numCols;
cellHeight = screenYpixels / numRows;

% Define a line length based on the cell size
lineLength = min(cellWidth, cellHeight) * 0.5;
lineWidth = 0.1 * lineLength; % Line width

anglesMatrix = {simple_matrix, matrix};
for a = 1:2
    % Loop through each cell and draw lines based on angles
    for row = 1:numRows
        for col = 1:numCols
            % Calculate center of each cell
            xCenter = (col - 0.5) * cellWidth;
            yCenter = (row - 0.5) * cellHeight;

            % Retrieve angles for this cell

            angles = anglesMatrix{a}{row, col};

            % Loop over each angle and draw the lines
            for angle = angles
                dx = (lineLength / 2) * cosd(angle);
                dy = (lineLength / 2) * sind(angle);
                lineCoords = [-dx, dy; dx, -dy]'; % Column vector for x and y coordinates

                % Draw the line
                Screen('DrawLines', window, lineCoords, lineWidth, 0, [xCenter, yCenter], 2);
            end
        end
    end
    % Flip to the screen
    Screen('Flip', window);

    % Wait for a key press
    KbStrokeWait;

    current_display = Screen('GetImage', window);
    name = strcat('C:\Users\idm\Desktop\Semester4\Internship\Matlab\example_', int2str(a), '.jpg');
    imwrite(current_display, name);
end
% Clear the screen
sca;



