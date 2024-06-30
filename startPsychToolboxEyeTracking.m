
function results = startPsychToolboxEyeTracking(data, folder_name, NbX, NbY, timeout, language)
global exp_folder
global ptb_drawformattedtext_oversize

% make subfolder for results
subfolder_name = [folder_name '\results'];
disp(subfolder_name)
mkdir(subfolder_name);
% Define key names
KbName('UnifyKeyNames');
escapeKey = KbName('ESCAPE');
rightKey = KbName('RightArrow');
leftKey = KbName('LeftArrow');
leftShift = KbName('LeftShift');
rightShift = KbName('RightShift');
leftKey = leftShift;
rightKey = rightShift;

% Hide the mouse cursor
%HideCursor; % TODO disable, for testing

bar_wh_ratio = 0.08; % bar width to height ratio
screenNumber = 0; % 0 default screen, 1 for external screen
color_bg = WhiteIndex(screenNumber);% colors for background and stim
color_stim = BlackIndex(screenNumber);
dotSizePix = 15; % Size of fixation point
text_size = 30;

% Open Psychtoolbox window
%PsychDebugWindowConfiguration(0,0.5); % TODO disable, makes window transparent
%Screen('Preference', 'SkipSyncTests', 1); % TODO disable, skips synchronization tests
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, color_bg);
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
xCenter = screenXpixels / 2;
yCenter = screenYpixels / 2;

% Maximum amount of jitter for each bar
jitter_x = round(0.1*screenXpixels/NbX);
jitter_y = round(0.1*screenYpixels/NbY);

% save infos
infos = {bar_wh_ratio, color_bg, color_stim, text_size, dotSizePix, jitter_x, jitter_y, xCenter, yCenter};
infos_table = cell2table(infos, 'VariableNames', {'bar_wh_ratio', 'color_bg', 'color_stim', ...
     'text_size', 'dotSizePix', 'jitter_x', 'jitter_y', 'xCenter', 'yCenter'});
saveData(infos_table, folder_name, 'screen_info.csv');

% Set text properties
Screen('TextSize', window, text_size);
Screen('TextColor', window, color_stim);
% Enable alpha blending for anti-aliasing
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


%% Instructions & calibration
if language == 1 % choose instructions language
    text_path = strcat(exp_folder, '\Images\instructions_en');
else
    text_path = strcat(exp_folder, '\Images\instructions_ger');
end


% Show the first 2 pictures/instructions
t = 1;
while t <= 2
    % Display first two images based on the current index
    img = imread(strcat(text_path, sprintf('\\instructions-0%0d.jpg', t)));
    displayImage(window, img, screenXpixels, screenYpixels);
    while 1 % Check for key presses to navigate through texts
        [~, ~, keyCode] = KbCheck;
        if keyCode(rightKey) % Check for right arrow press
            if t == 2 % Check if it's the second image
                WaitSecs(0.2); % Prevent rapid key presses
                t = t + 1; % Move index beyond the second to exit the loop
                break; % Exit the inner while loop
            else
                t = t + 1; % Otherwise, just move to the next image
                WaitSecs(0.2);
                break;
            end
        elseif keyCode(leftKey) && t > 1
            t = t - 1; % Navigate to the previous image
            WaitSecs(0.2);
            break;
        elseif keyCode(escapeKey)
            t = 3; % Move index beyond the second to exit the loop
            break; % Exit the loop if the escape key is pressed
        end
    end
    if t > 2
        break; % Break the outer loop if past the second image
    end
end


% Initialize Titta and set up the eye tracker 
try
    settings = Titta.getDefaults('Tobii Pro Nano');
    eye_tracker = Titta(settings);
    eye_tracker.init();
catch
    error('No eye trackers found. Make sure the eye tracker is connected and try again.');
end
% calibrate & safe results
try
    calibration = eye_tracker.calibrate(window);
    calibration_file = [subfolder_name, '\calibration_results_' , datestr(now, 'yyyymmdd_HHMM'), '.mat'];
    save(calibration_file, 'calibration');
catch
    error('Eye tracker couldnt be calibrated. Try again or restart experiment without eye tracking.');
end

% show rest of instructions
image_number = 3; % Start from the third image
while image_number <= 11 % until picture 11    
    try
        img = imread(strcat(text_path, sprintf('\\instructions-0%0d.jpg', image_number)));
    catch
        break; % Exit the loop if no more images are found
    end
    displayImage(window, img, screenXpixels, screenYpixels);
    while 1 % Check for key presses to navigate through images
        [~, ~, keyCode] = KbCheck;
        if keyCode(rightKey) % Check for right arrow press
            image_number = image_number + 1; % Move to the next image
            WaitSecs(0.2); % Prevent rapid key presses
            break;
        elseif keyCode(leftKey) && image_number > 3
            image_number = image_number - 1; % Navigate to the previous image
            WaitSecs(0.2);
            break;
        elseif keyCode(escapeKey)
            image_number = inf; % Move index beyond the last to exit the loop
            break; % Exit the loop if the escape key is pressed
        end
    end
    if image_number == inf
        break; % Break the outer loop if escape key was pressed
    end
end



%% Experiment loop starts here
disp('Start Experiment');
results_file_name = [subfolder_name, '\trial_results_' , datestr(now, 'yyyymmdd_HHMM'), '.mat'];
trial_results = table();

% TODO disable/change back
% for i = 1:size(data, 1)
for i = 1:4
%     % Display fixation dot & empty screen after
%     fixStartTime = GetSecs;
%     % Define the area around the fixation point where gaze must fall
%     fixationRadius = 50; % pixels, defines a circle radius around the fixation point
%     minFixationTime = 1.0; % seconds, the minimum time gaze must remain within the radius
%     samp_fix = drawFixationAndWait(window, xCenter, yCenter, color_stim, dotSizePix, eye_tracker, screenXpixels, screenYpixels, fixationRadius, minFixationTime);
%         %drawFixation(window, xCenter, yCenter, color_stim, dotSizePix); %without eyetrackgin
%         %WaitSecs(0.5); 
%     fixEndTime = GetSecs;
%     drawEmptyScreen(window, color_bg, 1)
%     WaitSecs(0.2); % wait before drawing stim
% %     end

    % Display fixation dot & empty screen after
    fixStartTime = GetSecs;
    eye_tracker.buffer.start('gaze');
    Screen('TextSize', window, text_size);
    displayText(window, 'Press any button to continue with the next trial.');
    while true
    [~, ~, keyCode] = KbCheck;
        if keyCode(leftKey) || keyCode(rightKey)
            break;
        end
    end
    Screen('Flip', window);
    drawFixation(window, xCenter, yCenter, color_stim, dotSizePix); %without eyetrackgin
    WaitSecs(0.5); 
    fixEndTime = GetSecs;
    drawEmptyScreen(window, color_bg, 1)
    WaitSecs(0.2); % wait before drawing stim
    samp_fix = eye_tracker.buffer.consumeN('gaze');

    current_data = data(i,:);
    % Start recording eye data
    eye_tracker.buffer.start('gaze');
    % Display stimulus
    trial_result = displayStim_new(window, current_data, bar_wh_ratio, jitter_x, jitter_y, ...
        screenXpixels, screenYpixels, color_stim, timeout);
        % Get eye data, since start/last getting eye data & safe
    samp = eye_tracker.buffer.consumeN('gaze');

    % Add data from fixation
    trial_result.fixStartTime = fixStartTime;
    trial_result.fixEndTime = fixEndTime;
    trial_result.fixTime = fixEndTime - fixStartTime;
    trial_result.samp_fix = samp_fix; 

    trial_result.samp = samp; % add eye data

    % save/overwrite each loop in case experiment crashes/gets aborted
    trial_results(i,:) = trial_result(1,:);
    save(results_file_name, 'trial_results');
    
    % Check for Escape key press to close the experiment
    if keyCode('ESCAPE')
        disp('Experiment Aborted: Escape Key was pressed')
        break;
    end
end
% Clean up: stop eye tracker, save results, and close
eye_tracker.buffer.stop('gaze');
eye_tracker.deInit();
results = trial_results;
%% End
% End screen
% TODO CHANGE BACK FOR EXPERIMENT
%if i == size(data, 1)
if i == 4
    Screen('TextSize', window, text_size);
    if language == '2' % choose instructions language     
        displayText(window, 'Experiment Abgeschlossen\n\n Vielen Dank fürs Mitmachen! <3');
    else
        displayText(window, 'Experiment Complete\n\n Thank you for your participation! <3');
    end
    WaitSecs(2);
    KbStrokeWait(-1);  % Wait to end with key press
end

disp('Experiment Finished')
sca; % Close the window
end
