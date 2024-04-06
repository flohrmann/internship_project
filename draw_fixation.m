function draw_fixation(window, xCenter, yCenter, color, dotSizePix)

Screen('DrawDots', window, [xCenter; yCenter], dotSizePix, color, [], 2);
Screen('Flip', window);



% % Set the coordinates (these are all relative to zero we will let
% % the drawing routine center the cross in the center of our monitor for us)
% xCoords = [-size size 0 0];
% yCoords = [0 0 -size size];
% allCoords = [xCoords; yCoords];
% 
% % Draw the fixation cross in white, set it to the center of our screen and
% % set good quality antialiasing
% Screen('DrawLines', window, allCoords,...
%     width, color, [xCenter yCenter], 2);
% 
% % Flip to the screen
% Screen('Flip', window);
% 
% % Wait for a key press
KbStrokeWait;
