function startPsychToolbox(trial_data)
    % Clear the workspace and the screen
    sca;
    close all;

    % Open Psychtoolbox window
    Screen('Preference', 'SkipSyncTests', 1); % Skip synchronization tests
    PsychDefaultSetup(2);

    % Get the screen numbers
    screens = Screen('Screens');
    
    % Draw to the external screen if avaliable
    screenNumber = max(screens);
    
    % Define black and white
    white = WhiteIndex(screenNumber);
    black = BlackIndex(screenNumber);

    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, white); % Open window with a white background

    % Set text properties
    Screen('TextSize', window, 24);
    Screen('TextColor', window, [0 0 0]); % Black text color

    % Get the size of the on screen window
    [screenXpixels, screenYpixels] = Screen('WindowSize', window);

    % Define texts
    texts = {'Text 1: This is the first text.', ...
             'Text 2: This is the second text.', ...
             'Text 3: This is the third text.'};

    % Display each text
    for i = 1:numel(texts)
        % Display text
        displayText(window, texts{i});

        % Wait for a key press to continue
        while KbCheck; end % Wait for all keys to be released
        while ~KbCheck; end % Wait for a key press
        [~, ~, keyCode] = KbCheck; % Check for key press
        
        % Check for Escape key press to close the experiment
        if keyCode(KbName('Escape'))
            break;
        end
    end

    % Close the window
    sca;
end

