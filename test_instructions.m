% Define screen parameters
exp_folder = 'C:\Users\flohrmann\Documents\MATLAB\internship_project';
Screen('Preference', 'SkipSyncTests', 1);
%PsychDebugWindowConfiguration(0,0.5);
screenNumber = 0;
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, 0);
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
Screen('TextSize', window, 24);

% Define text content
texts = {
    'Welcome to the experiment. Press the right arrow key to start.',...
    'This is slide 2. Press the right arrow key to continue.', ...
    'This is slide 3. Press the right arrow key to continue.',...
    'Slide 4 with an image.',...
    'Slide 5 with another image.',...
    'This is slide 6. Press the right arrow key to continue.',...
    'Thank you for participating. Press the right arrow key to finish.'...
    };

% Define keys
KbName('UnifyKeyNames');
rightKey = KbName('RightArrow');
leftKey = KbName('LeftArrow');
escapeKey = KbName('ESCAPE');

t = 1;
while t <= numel(texts)
    % Display text or image based on the current index
    if t == 4 % Add images to certain 'slides'
        img = imread(strcat(exp_folder,'\Images\example_jittered1.jpg'));
        displayImage(window, texts{t}, img, screenXpixels, screenYpixels);
    elseif t == 5
        img = imread(strcat(exp_folder,'\Images\example_jittered2.jpg'));
        displayImage(window, texts{t}, img, screenXpixels, screenYpixels);
    else
        DrawFormattedText(window, texts{t}, 'center', 'center');
        Screen('Flip', window);
        
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
end
