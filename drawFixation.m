function drawFixation(window, xCenter, yCenter, color, dotSizePix)

Screen('DrawDots', window, [xCenter; yCenter], dotSizePix, color, [], 2);
Screen('Flip', window);
KbStrokeWait(-1);

