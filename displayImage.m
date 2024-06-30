function displayImage(window, img, screenXpixels, screenYpixels)
        texture = Screen('MakeTexture', window, img);
        % Calculate the position of the scaled image to center it
        dstX = (screenXpixels - screenXpixels) / 2;
        dstY = (screenYpixels - screenYpixels) / 2;
        dstRect = [dstX dstY dstX + screenXpixels dstY + screenYpixels];
        % Add text on top of the image
        %Screen('TextSize', window, 30);
        %DrawFormattedText(window, text, 'center', dstY * 0.8, 0);
        % Draw the scaled texture to the window
        Screen('DrawTexture', window, texture, [], dstRect);
        %Screen('DrawTexture', window, texture, []);
        Screen('Flip', window);
end