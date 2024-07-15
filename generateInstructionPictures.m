% Define the matrices
a = { 
    [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0];%, [45 90], [45 0];
    [45 0], [45 90], [135 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90];%, [45 0], [45 90];
    [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0];%, [45 90], [45 0];
    [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90];%, [45 0], [45 90];
    [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0];%, [45 90], [45 0];
    [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90];%, [45 0], [45 90];
    [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0];%, [45 90], [45 0];
    [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90];%, [45 0], [45 90];
    [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0];%, [45 90], [45 0];
    %[45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90];
    %[45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0], [45 90], [45 0];
};

a_simple = { 
    [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45];%, [45], [45];
    [45], [45], [135], [45], [45], [45], [45], [45], [45], [45], [45], [45];%, [45], [45];
    [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45];%, [45], [45];
    [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45];%, [45], [45];
    [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45];%, [45], [45];
    [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45];%, [45], [45];
    [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45];%, [45], [45];
    [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45];%, [45], [45];
    [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45];%, [45], [45];
    %[45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45];
    %[45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45], [45];
};

b_simple = { 
    [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135];%, [135], [135];
    [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135];%, [135], [135];
    [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135];%, [135], [135];
    [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135];%, [135], [135];
    [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135];%, [135], [135];
    [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135];%, [135], [135];
    [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135];%, [135], [135];
    [135], [135], [135], [135], [135], [135], [135], [135], [135], [70], [135], [135];%, [135], [135];
    [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135];%, [135], [135];
    %[135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135];
    %[135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135], [135];
};


b = { 
    [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0];%, [135 90], [135 0];
    [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90];%, [135 0], [135 90];
    [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0];%, [135 90], [135 0];
    [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90];%, [135 0], [135 90];
    [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0];%, [135 90], [135 0];
    [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90];%, [135 0], [135 90];
    [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0];%, [135 90], [135 0];
    [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [20 0], [135 0], [135 90];%, [135 0], [135 90];
    [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0];%, [135 90], [135 0];
    %[135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90];
    %[135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0], [135 90], [135 0];
};
       
PsychDefaultSetup(2);
screens = Screen('Screens');
screenNumber = 0; % 0 default screen, 1 for external screen
[window, windowRect] = PsychImaging('OpenWindow', 0, 1);
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Define grid size
num_rows = 9;
num_cols = 11;
cell_width = screenXpixels / num_cols;
cell_height = screenYpixels / num_rows;

% Define a line length based on the cell size
bar_wh_ratio = 0.08; 
lineLength = round(min(cell_width, cell_height) * 0.5);
lineWidth = bar_wh_ratio * lineLength;

% Define jitter
jitter_x = round(0.1 * screenXpixels / num_cols);
jitter_y = round(0.1 * screenYpixels / num_rows);

% Precompute jitter
x_jitter_matrix = randi([-jitter_x, jitter_x], num_rows, num_cols);
y_jitter_matrix = randi([-jitter_y, jitter_y], num_rows, num_cols);

anglesMatrix = {a, a_simple, b, b_simple};


for a = 1:4
    % Loop through each cell and draw lines based on angles
    for row = 1:num_rows
        for col = 1:num_cols
            x_jitter = x_jitter_matrix(row, col);
            y_jitter = y_jitter_matrix(row, col);

            % Calculate center of each cell & add jitter
            x_center = ((col - 0.5) * cell_width) + x_jitter;
            y_center = ((row - 0.5) * cell_height) + y_jitter;

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

    % Save the current display
    current_display = Screen('GetImage', window);
    name = strcat('C:\Users\flohrmann\Documents\MATLAB\internship_project\Images\example_jittered_', int2str(a), '.jpg');
    imwrite(current_display, name);
end

Screen('Flip', window);

crossSizePix = 60;
crossLineWidthPix = 15;
crossColour = [255 0 0]; % Red 
fixStartTime = drawFixation(window, screenXpixels/2, screenYpixels/2, crossColour, crossLineWidthPix, crossSizePix);
KbStrokeWait;
current_display = Screen('GetImage', window);
name = strcat('C:\Users\flohrmann\Documents\MATLAB\internship_project\Images\example_fixation_cross.jpg');
imwrite(current_display, name);



% Screen('Flip', window);
% dotSizePix = 15; % Size of fixation point
% xCenter = screenXpixels / 2;
% yCenter = screenYpixels / 2;
% black = BlackIndex(screenNumber);
% Screen('DrawDots', window, [xCenter, yCenter], dotSizePix, black, [], 2);
% Screen('Flip', window);
% KbStrokeWait;
% current_display = Screen('GetImage', window);
% name = strcat('C:\Users\flohrmann\Documents\MATLAB\internship_project\Images\example_fixation.jpg');
% imwrite(current_display, name);

% Clear the screen
sca;