function startPsychToolbox(data, folder_name, NbX, NbY, timeout, language)
% doesnt use eyetracking
% make subfolder for results
subfolder_name = [folder_name '\results'];
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
dotSizePix = 8; % Size of fixation point
text_size = 30;

% Open Psychtoolbox window
PsychDebugWindowConfiguration(0,0.5); % TODO disable, makes window transparent
Screen('Preference', 'SkipSyncTests', 1); % TODO disable, skips synchronization tests
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

%% Instructions
texts_en = {'Welcome :)', ...
    'In this experiment you will be presented a grid of bars \n\n Press the Left Arrow Key if you can see a uniquely oriented bar on the Left side of the screen \n and the Right Arrow Key if you can see a uniquely oriented bar on the Right side of the screen', ...
    'Lets look at a possible images you will see later', ...
    'The trials can look like this for example. Can you find the uniquely oriented bar?',...
    'The trials can also look like this. Can you still see the uniquely oriented bar?', ...
    'Before each trial you will see a fixation dot in the middle of the screen, please remember to look at it for as long as it is there \n\n The next trial will start automatically after half a second/if you look at it for around 2 seconds.', ...
    'During the experiment once you have found the uniquely oriented bar remember to press the corresponding button \n\n Please try to be as fast and as accurate as possible!', ...
    'Do you have any questions? \n \n If you feel ready: Press the Right Arrow Key To Begin The Experiment'};
texts_ger = {'Willkommen :)', ...
    'In diesem Experiment wird Ihnen ein Gitter aus Balken präsentiert. \n\n Drücken Sie die linke Shift Taste, wenn Sie links auf dem Bildschirm einen einzigartig orientierten Balken sehen können, \n und die rechte Shift-Taste, wenn Sie rechts auf dem Bildschirm einen einzigartig orientierten Balken sehen können', ...
    'Lassen Sie uns einige mögliche Bilder anschauen, die Sie später sehen werden', ...
    'Ein Durchlauf könnte könnte zum Beispiel so aussehen. Können Sie den einzigartig orientierten Balken finden?', ...
    'Ein Durchlauf könnte auch so aussehen. Können Sie den einzigartig orientierten Balken immer noch finden?', ...
    'Vor jedem Versuch sehen Sie einen Fixierungspunkt in der Mitte des Bildschirms, bitte schauen Sie so lange darauf, wie er dort ist \n\n Der nächste Versuch beginnt automatisch nach einer halben Sekunde.', ...
    'Während des Experiments, sobald Sie den einzigartig orientierten Balken gefunden haben, denken Sie daran, die entsprechende Taste zu drücken \n (linke/rechte Shift-Taste) \n Bitte versuchen Sie, so schnell und so genau wie möglich zu sein!', ...
    'Haben Sie irgendwelche Fragen? \n \n Wenn Sie bereit sind: Drücken Sie die rechte Shift-Taste, um mit dem Experiment zu beginnen'};

if language == 2 % choose instructions language
    texts = texts_ger;
else
    texts = texts_en;
end

t = 1;
while t <= numel(texts)
    % Display text or image based on the current index
    if t == 4 % Add images to certain 'slides'
        img = imread('C:\Users\idm\Desktop\Semester4\Internship\Matlab\ExperimentRepo\Images\example_jittered1.jpg'); 
        displayImage(window, texts{t}, img, screenXpixels, screenYpixels);
    elseif t == 5
        img = imread('C:\Users\idm\Desktop\Semester4\Internship\Matlab\ExperimentRepo\Images\example_jittered2.jpg'); 
        displayImage(window, texts{t}, img, screenXpixels, screenYpixels);
    else
        displayText(window, texts{t}); % Display text only
    end
    while 1 % Check for key presses to navigate through texts
        [~, ~, keyCode] = KbCheck;
        if keyCode(rightKey) % Check for right arrow press
            if t == numel(texts) % Check if it's the last text
                WaitSecs(0.2); % Prevent rapid key presses
                t = t + 1; % Move index beyond the last to exit the loop
                break; % Exit the inner while loop
            else
                t = t + 1; % Otherwise, just move to the next text
                WaitSecs(0.2); 
                break; 
            end
        elseif keyCode(leftKey) && t > 1
            t = t - 1; % Navigate to the previous text
            WaitSecs(0.2); 
            break; 
        elseif keyCode(escapeKey)
            t = numel(texts) + 1; % Move index beyond the last to exit the loop
            break; % Exit the loop if the escape key is pressed
        end
    end
    if t > numel(texts) 
        break; % Break the outer loop if past the last text
    end
end

%% Start experiment
disp('Start Experiment');
results_file_name = [subfolder_name, '\trial_results_' , datestr(now, 'yyyymmdd_HHMM'), '.mat'];
trial_results = table();
% TODO disable/change back
% for i = 1:size(data, 1)
for i = 1:4
    % Display fixation dot & empty screen after
    drawFixation(window, xCenter, yCenter, color_stim, dotSizePix);
    WaitSecs(0.5); % show for 0.5 seconds
    drawEmptyScreen(window, color_bg, 1)
    WaitSecs(0.2); % wait before drawing stim

    current_data = data(i,:);
    % Display stimulus
    trial_result = displayStim_new(window, current_data, bar_wh_ratio, jitter_x, jitter_y, ...
        screenXpixels, screenYpixels, color_stim, timeout);
         
    % save/overwrite each loop in case experiment crashes/gets aborted
    trial_results(i,:) = trial_result(1,:);
    save(results_file_name, 'trial_results');

    % Check for Escape key press to close the experiment
    if keyCode('ESCAPE')
        disp('Experiment Aborted: Escape Key was pressed')
        break;
    end
end

%% End screen
% TODO CHANGE BACK FOR EXPERIMENT
%if i == size(data, 1)
if i == 4
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