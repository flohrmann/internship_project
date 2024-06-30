% ONLY WORKS WITH SKIPSYNCTESTS DISABLED
screenNumber = 0; % 0 default screen, 1 for external screen
color_bg = WhiteIndex(screenNumber);% colors for background and stim
subfolder_name = 'C:\Users\flohrmann\Documents\MATLAB\internship_project\test';
Screen('Preference', 'SkipSyncTests', 1); % TODO disable, skips synchronization tests
[window, ~] = PsychImaging('OpenWindow', screenNumber, color_bg);
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
%screenXpixels = screenXpixels\2;
%screenYpixels = screenYpixels\2;
% Initialize Titta and set up the eye tracker 
%try
    settings = Titta.getDefaults('Tobii Pro Nano');
    eye_tracker = Titta(settings);
    eye_tracker.init();
%catch
%    error('No eye trackers found. Make sure the eye tracker is connected and try again.');
%end
% calibrate & safe results
try
    calibration = eye_tracker.calibrate(window);
    calibration_file = [subfolder_name, '\calibration_results.mat'];
    save(calibration_file, 'calibration');
catch
    error('Eye tracker couldnt be calibrated. Try again or restart experiment without eye tracking.');
end
% for testing only: fanis calibration
% retrieve_calibration_data('C:\Users\idm\Desktop\Semester4\Internship\Matlab\ExperimentRepo\results\calibration_results_20240513_1144')	;
% apply_calibration_data(calibration_file);	
sca; % Close the window



%% test 2
Screen('Preference', 'SkipSyncTests', 1); % TODO disable, skips synchronization tests
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, color_bg);

% Add Titta to MATLAB path if not already added
settings = Titta.getDefaults('Tobii Pro Nano');
eye_tracker = Titta(settings);
eye_tracker.init();

eye_tracker.calibrate(window);

sca; % Close the window

