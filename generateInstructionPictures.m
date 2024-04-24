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
num_rows = 9;
num_cols = 12;
cell_width = screenXpixels / num_cols;
cell_height = screenYpixels / num_rows;

% Define a line length based on the cell size
lineLength = min(cell_width, cell_height) * 0.5;
lineWidth = 0.1 * lineLength; % Line width

% Define jittter
jitter_x = round(0.1*screenXpixels/num_cols);
jitter_y = round(0.1*screenYpixels/num_rows);

anglesMatrix = {simple_matrix, matrix};
for a = 1:2
    % Loop through each cell and draw lines based on angles
    for row = 1:num_rows
        for col = 1:num_cols

            x_jitter = randi([-jitter_x, jitter_x]);
            y_jitter = randi([-jitter_y, jitter_y]);

            % Calculate center of each cell & add jitter
            x_center = ((col - 0.5) * cell_width) + x_jitter;
            y_center = ((row - 0.5) * cell_height)+ y_jitter;

            % Retrieve angles for this cell
            angles = anglesMatrix{a}{row, col};

            % Loop over each angle and draw the lines
            for angle = angles
                dx = (lineLength / 2) * cosd(angle);
                dy = (lineLength / 2) * sind(angle);
                lineCoords = [-dx, dy; dx, -dy]'; % Column vector for x and y coordinates

                % Draw the line
                Screen('DrawLines', window, lineCoords, lineWidth, 0, [x_center, y_center], 2);
            end
        end
    end
    % Flip to the screen
    Screen('Flip', window);

    % Wait for a key press
    KbStrokeWait;

    current_display = Screen('GetImage', window);
    name = strcat('C:\Users\idm\Desktop\Semester4\Internship\Matlab\ExperimentRepo\Images\example_jittered', int2str(a), '.jpg');
    imwrite(current_display, name);
end
% Clear the screen
sca;



