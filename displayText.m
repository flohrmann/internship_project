function displayText(window, text)
    % Display text on the Psychtoolbox window
    Screen('TextSize', window, 24);
    DrawFormattedText(window, text, 'center', 'center');
    Screen('Flip', window);
end