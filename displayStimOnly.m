function [data, stop] = displayStimOnly(window, data, bar_wh_ratio, jitter_x, jitter_y, screenXpixels, screenYpixels, color_stim, timeoutDuration, blankStartTime, monitorFlipInterval, color_bg)

% Define the keyboard keys that are listened for
% TODO add more keys for slipping
escapeKey = KbName('ESCAPE');
%leftKey = KbName('LeftArrow');
%rightKey = KbName('RightArrow');
leftKey = KbName('LeftShift');
rightKey = KbName('RightShift');
% leftKey = leftShift;
% rightKey = rightShift;

current_stim = data.AngleMatrix{1};
current_target_pos = data.TargetSide{1};

num_cols = size(current_stim, 2);
num_rows = size(current_stim, 1);
cell_width = screenXpixels / num_cols;
cell_height = screenYpixels / num_rows;

% Define a line length based on the cell size
line_length = min(cell_width, cell_height) * 0.5;
line_width = bar_wh_ratio * line_length;

%safe
data.cell_width = cell_width;
data.cell_height = cell_height;
data.line_length = line_length;
data.line_width = line_width;

% Initialize jitter storage
x_jitters = zeros(num_rows, num_cols);
y_jitters = zeros(num_rows, num_cols);
x_centers = zeros(num_rows, num_cols);
y_centers = zeros(num_rows, num_cols);

%trial_coords_list = cell(num_rows, num_cols);


% flip to another blank screen just to make sure eyetracker starts before stimulation starts
% Screen('FillRect', window, color_bg); 
% Screen('Flip', window, blankStartTime+0.2-1.5*monitorFlipInterval); % wait 0.2 sec - 1 screen flip
%  
% Start recording eye data if with eyetracking
% if ~continue_without_eyetracking
%     eye_tracker.buffer.start('gaze');
% end   


for row = 1:num_rows %y
    for col = 1:num_cols %x
        % Get random jitter for this position
        x_jitter = randi([-jitter_x, jitter_x]);
        y_jitter = randi([-jitter_y, jitter_y]);
        % Store jitter values
        x_jitters(row, col) = x_jitter;
        y_jitters(row, col) = y_jitter;
        
        % Calculate center of each cell & add jitter
        x_center = ((col - 0.5) * cell_width) + x_jitter;
        y_center = ((row - 0.5) * cell_height)+ y_jitter;
        x_centers(row, col) = x_center; 
        y_centers(row, col) = y_center; % safe
        angles = current_stim{row, col}; % Retrieve angles for this cell
        %cell_line_coords = [];
        % Loop over each angle and draw the lines
        for angle = angles
            dx = (line_length / 2) * cosd(angle); % cosd/sind can use degrees, so dont
            dy = (line_length / 2) * sind(angle); % need to convert to radians like for cos/sin
            line_coords = [-dx, dy; dx, -dy]'; % Column vector for x and y coordinates
            % Draw the line
            Screen('DrawLines', window, line_coords, line_width, color_stim, [x_center, y_center], 2);
            %cell_line_coords = [cell_line_coords, line_coords];
        end
        %trial_coords_list{row, col} = cell_line_coords;
    end
end
%data.TrialCenterStimCoords = {trial_coords_list};
      
% Start time of the trial
trial_start_time = blankStartTime+0.2-0.5*monitorFlipInterval;
% Flip to the screen for each trial
StimulusOnsetTime = Screen('Flip', window, trial_start_time); % to make sure it doesnt flip a frame too late
    
response = 'none'; % Default response if no key is pressed

%%% Code from https://peterscarfe.com/poserCuingExperiment.html :
% Now we wait for a keyboard button signaling the observers response.
% The left arrow key signals a "left" response and the right arrow key
% a "right" response. You can also press escape if you want to exit the
% program
respToBeMade = true;
stop = 0;
while respToBeMade && (GetSecs - trial_start_time < timeoutDuration)
%while respToBeMade
    [keyIsDown,secs, keyCode] = KbCheck(-1);
    if keyCode(escapeKey)
        ShowCursor;
        stop = 1;
        sca;
        return
    elseif keyCode(leftKey)
        response = 'L';
        respToBeMade = false;
        stop = 0;
    elseif keyCode(rightKey)
        response = 'R';
        respToBeMade = false;
        stop = 0;
    end
end
% Get eye data, since start/last getting eye data & safe
% if ~continue_without_eyetracking
%     % Get eye data, since start/last getting eye data & save
%     samp = eye_tracker.buffer.consumeN('gaze');
% else
%     samp = struct('deviceTimeStamp', [], 'systemTimeStamp', [], 'left', struct(), 'right', struct()); % No eye-tracking data
% end
%     
%trial_resp_time = GetSecs; % number of seconds since system start up

% flip to blank screen to get system time instead
Screen('FillRect', window, color_bg); 
trial_resp_time = Screen('Flip', window);
%  
% noise
fs = 5000; t = 0:0.00002:0.02;
LowToneSoundwave =  sin(2*pi*fs/2*t);
sound(LowToneSoundwave, fs); pause(0.2);
%

% Work out if the location of the target was identified corrcetly
if current_target_pos == response
    correctness = 1;
elseif current_target_pos ~= response
    correctness = 0;
end

%%%

% Store the jitters in the data struct
data.x_jitters = {x_jitters};
data.y_jitters = {y_jitters};
data.x_centers = {x_centers};
data.y_centers = {y_centers};

% Record user response, response time
data.target = current_target_pos;
data.response = {response};
data.correct = correctness;
data.blankStartTime = blankStartTime;
data.trialStartTime = trial_start_time;
data.StimulusOnsetTime = StimulusOnsetTime;
data.trialEndTime = trial_resp_time;
data.rt = trial_resp_time - StimulusOnsetTime;
%data.eyeTrial = samp; % eyetracker 



end
