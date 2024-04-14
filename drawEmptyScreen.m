function drawEmptyScreen(window, color_bg, duration)

Screen('FillRect', window, color_bg);
Screen('Flip', window);
WaitSecs(duration);% Wait 0.3 seconds