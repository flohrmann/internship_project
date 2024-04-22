function displayText(window, text)
    % Display text on the Psychtoolbox window
    %Screen('TextSize', window, text_size);
    DrawFormattedText(window, text, 'center', 'center');
    Screen('Flip', window);
end