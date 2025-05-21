function [trial_results, samp] = startPsychToolboxEyeTrackingFixationTracked(data, folder_name, NbX, NbY, timeout, language, continue_without_eyetracking)
global exp_folder
global ptb_drawformattedtext_oversize

%% make subfolder for results
subfolder_name = [folder_name '\results'];
disp(subfolder_name)
mkdir(subfolder_name);
trial_results = table(); % to collect all the data/results (also eye tracker snippets)
samp = struct(); % to collect all the eye tracker snippets

%% Define key names
KbName('UnifyKeyNames');
escapeKey = KbName('ESCAPE');
rightKey = KbName('RightArrow');
leftKey = KbName('LeftArrow');
leftShift = KbName('LeftShift');
rightShift = KbName('RightShift');
leftKey = leftShift;
rightKey = rightShift;
spaceKey = KbName('space');
deleteKey = KbName('Delete'); % ENTF key to stop eyetracking

%continue_without_eyetracking = false; % flag to check if the experiment should continue without eye-tracking

%% Hide the mouse cursor
HideCursor; % TODO disable, for testing


%% What to show:
bar_wh_ratio = 0.08; % bar width to height ratio
screenNumber = 0; % 0 default screen, 1 for external screen
color_bg = WhiteIndex(screenNumber);% colors for background and stim
color_stim = BlackIndex(screenNumber);

crossSizePix = 60; % fixation cross size
crossLineWidthPix = 15;
crossColour = [255 0 0]; % Red

text_size = 40;

minFixationTime = 0.5; % how many seconds to look at fixation cross for experiment to continue
fixationRadius = 200; % radius how many pixels away from center counts as fixating the cross
wait_time_s_after_button_press = 0.5; % how many seconds to stay in stimulation screen after button was pressed

%% Open Psychtoolbox window
%PsychDebugWindowConfiguration(0,0.5); % TODO disable, makes window transparent
Screen('Preference', 'SkipSyncTests', 1); % TODO disable, skips synchronization tests
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, color_bg);
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
xCenter = screenXpixels / 2;
yCenter = screenYpixels / 2;

% Maximum amount of jitter for each bar
jitter_x = round(0.1*screenXpixels/NbX);
jitter_y = round(0.1*screenYpixels/NbY);

% save infos % TODO save cross
infos = {bar_wh_ratio, color_bg, color_stim, text_size, jitter_x, jitter_y, xCenter, yCenter, wait_time_s_after_button_press, ...
        crossSizePix , crossLineWidthPix, crossColour};
infos_table = cell2table(infos, 'VariableNames', {'bar_wh_ratio', 'color_bg', 'color_stim', ...
                    'text_size', 'jitter_x', 'jitter_y', 'xCenter', 'yCenter', 'wait_time_s_after_button_press', ...
                    'crossSizePix', 'crossLineWidthPix', 'crossColour'});
saveData(infos_table, folder_name, 'screen_info.csv');

monitorFlipInterval = Screen('GetFlipInterval', window); % to get inter flip interval for break

% Set text properties
Screen('TextSize', window, text_size);
Screen('TextColor', window, color_stim);
% Enable alpha blending for anti-aliasing
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%% Instructions Part 1
if language == 1 % choose instructions language
    text_path = strcat(exp_folder, '\Images\instructions_en');
else
    text_path = strcat(exp_folder, '\Images\instructions_ger');
end


% Show the first 2 pictures/instructions
t = 1;

if ~continue_without_eyetracking
    t_end = 2;
else
    t_end = 1;
end

while t <= t_end
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
            disp('Experiment Aborted During First Instructions');
            sca; % Close the window
            return; % Close function
        end
    end
    if t > 2
        break; % Break the outer loop if past the second image
    end
end


%% Initialize Titta and set up the eye tracker
if ~continue_without_eyetracking
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
        %error('Eye tracker couldnt be calibrated. Try again or restart experiment without eye tracking.');
        % Calibration failed, allow user to continue without eye-tracking
        Screen('TextSize', window, text_size);
        displayText(window, 'Calibration failed. Press ENTF to continue without eye-tracking or the Spacebar to try again.');
        %Screen('Flip', window);
        WaitSecs(1);
        while 1
            [~, ~, keyCode] = KbCheck;
            if keyCode(deleteKey)
                continue_without_eyetracking = true;
                break;
            elseif keyCode(spaceKey)
                try
                    settings = Titta.getDefaults('Tobii Pro Nano');
                    eye_tracker = Titta(settings);
                    eye_tracker.init();
                    calibration = eye_tracker.calibrate(window);
                    calibration_file = [subfolder_name, '\calibration_results_' , datestr(now, 'yyyymmdd_HHMM'), '.mat'];
                    save(calibration_file, 'calibration');
                catch
                    Screen('TextSize', window, text_size);
                    displayText(window, 'Calibration failed. Press ENTF to continue without eye-tracking.');
                end
                break;
            end
        end
    end
else
    % do nothing
end

% %% DEBUG: instead of calibration do this while debugging:
% settings = Titta.getDefaults('Tobii Pro Nano');
% eye_tracker = Titta(settings);
% eye_tracker.init();
% %%%


%% start eyetracker 
if ~continue_without_eyetracking
    eye_tracker.buffer.start('gaze');
end


% Instructions Part 2
image_number = 3; % Start from the third image
while image_number <= 16 % until picture 16
    try % einstellige nr
        img = imread(strcat(text_path, sprintf('\\instructions-0%0d.jpg', image_number)));
    catch
        try % zweistellige nr
            img = imread(strcat(text_path, sprintf('\\instructions-%02d.jpg', image_number)));
        catch
            break; % Exit the loop if no more images are found
        end
        %break; % Exit the loop if no more images are found
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
            disp('Experiment Aborted During Instructions');
            sca; % Close the window
            return; % Close function
        end
    end
end

%% Experiment Loop
disp('Start Experiment');
results_file_name = [subfolder_name, '\trial_results_' , datestr(now, 'yyyymmdd_HHMM'), '.mat'];
eye_results_file_name = [subfolder_name, '\eyetracking_results_' , datestr(now, 'yyyymmdd_HHMM'), '.mat'];

for trial = 1:4% size(data, 1)
    current_data = data(trial,:);
    
    %% Display Trial Start Instructions until button press
    Screen('FillRect', window, color_bg);
    Screen('TextSize', window, text_size);
    DrawFormattedText(window, 'Put your two index fingers on the shift buttons to get ready \n\n\n and press either of them to continue with the next trial', 'center', 'center');
    trialStartTime = Screen('Flip', window);
    %displayText(window, 'Put your two index fingers on the shift buttons to get ready \n\n\n and press either of them to continue with the next trial');
    
    %% Option to disable eye tracking during this menu
    while true
        [~, ~, keyCode] = KbCheck;
        if keyCode(leftKey) || keyCode(rightKey)
            break;
        elseif keyCode(deleteKey)
            disp('Experiment switched to no eyetracking: Entf/Delete Key was pressed')
            continue_without_eyetracking = true;
            % safe continous eyetracker data
            samp = eye_tracker.buffer.consumeN('gaze');
            save(eye_results_file_name, 'samp');
            eye_tracker.buffer.stop('gaze');
            eye_tracker.deInit();
            break;
        elseif keyCode(escapeKey)
            disp('Experiment Aborted: Escape Key was pressed')
            % safe continous eyetracker data
            %continue_without_eyetracking = true;
            if ~continue_without_eyetracking
                samp = eye_tracker.buffer.consumeN('gaze');
                save(eye_results_file_name, 'samp');
                eye_tracker.buffer.stop('gaze');
                eye_tracker.deInit();
                sca; % close window
                return; % close function
            else
                sca; % close window
                return; % close function
            end
        end
    end
    
    %% Display fixation cross & check for gaze/ if focusing on cross
    %fixStartTime = drawFixation(window, xCenter, yCenter, crossColour, crossLineWidthPix, crossSizePix);
    [samp_start, samp_fix, fixStartTime] = drawFixationAndWait(window, xCenter, yCenter, crossColour, crossLineWidthPix, crossSizePix, eye_tracker, screenXpixels, screenYpixels, fixationRadius, minFixationTime);

    % time of end of fixation > should be approximaltely the same as blankStartTime    
    %fixEndTime = GetSecs;
    
    % restart eyetracker
    eye_tracker.buffer.start('gaze');
    
    %% Display empty screen (for 0.2 seconds)
    Screen('FillRect', window, color_bg); % blank screen
    blankStartTime = Screen('Flip', window); % use the next possible time flip  , fixStartTime+0.5-0.5*monitorFlipInterval); % wait 0.5 s
    
    %% Display stimulus (after waiting 0.2 s for screen flip)
    % after pressing button: makes beep noise and waits for 0.2 s for beep to play/end
    [trial_result, stop] = displayStimOnly(window, current_data, bar_wh_ratio, jitter_x, jitter_y, screenXpixels, screenYpixels,...
        color_stim, timeout, blankStartTime, monitorFlipInterval, color_bg);
    
    % Check for Escape key press to close the experiment
    if stop == 1
        disp('Experiment Aborted: Escape Key was pressed')
        break;
    end
    
    %% wait for 0.5 seconds on this screen just so that it doesnt change in luminance for gaze analyis
    WaitSecs(wait_time_s_after_button_press);

    %% safe trial results    
    
    % get rest of trial eye tracking data
    samp_blank_stim = eye_tracker.buffer.consumeN('gaze');
        
    % Add all the timestamps
    trial_result.trialStartTime = trialStartTime; % timestamp trial start
    trial_result.fixStartTime = fixStartTime; % timestamp fixation start
    trial_result.UserInPressAnyButtonScreenDuration = fixStartTime - trialStartTime; % duration until fixation started this trial (pause of user basically)
    
    trial_result.eyeStart = samp_start; % eye start of trial/before fixation
    trial_result.eyeFix = {samp_fix}; % eye fixation
    trial_result.eyeBlankAndStim = samp_blank_stim; % eye blank screen and stimulation
    
    % doppelt gemoppelt: same for eyetracking data
    samp(trial).sampStart = samp_start;
    samp(trial).sampFix = samp_fix; 
    samp(trial).sampBlankStim = samp_blank_stim;
    
    % save/overwrite each loop in case experiment crashes/gets aborted
    trial_results(trial,:) = trial_result(1,:);
    save(results_file_name, 'trial_results');
    save(eye_results_file_name, 'samp'); 
    
    clear trial_result
end
% catch
%     disp('Experiment Crashed')
% end



%% End screen
% clean up: stop eye tracker and close
eye_tracker.deInit();

% TODO CHANGE BACK FOR EXPERIMENT
if trial == size(data, 1)
    %if i == 4
    Screen('FillRect', window, color_bg);
    Screen('TextSize', window, text_size);
    if language == 2 % choose instructions language
        displayText(window, 'Experiment Abgeschlossen\n\n Vielen Dank fürs Mitmachen!');
    else
        img = imread(strcat(text_path, '\instructions-18.jpg'));
        displayImage(window, img, screenXpixels, screenYpixels);
    end
    WaitSecs(2);
    KbStrokeWait(-1);  % Wait to end with key press
end

disp('Experiment Finished')
sca; % Close the window
end
