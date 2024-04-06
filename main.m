function main()
    % Clear the workspace and the screen
    sca;
    close all;
    clear;
    %% Create a folder with the current date and time 
    folder_name = ['C:\Users\idm\Desktop\Semester4\Internship\Matlab\Data\test_' datestr(now, 'yyyymmdd_HHMMSS')];
    mkdir(folder_name);
    
    %% Parameters for this run
    n_rows = 12; % number stimulus x axis
    n_columns = 20; % number stimulus y axis
    grid_visual_angle = [34, 46]; % in degrees
    stim_size = [0.12, 1.1]; % in degrees
    ec_circle = 15; % circle of stim pos in degrees
    ec_min = 12; % minimum horizontal eccentricity
    fix_stim_dia = 0.3; % in degrees
    %stimulus_brightness = 48; % in cd/m² -> dont need, white bg & black stim?
    %conditions = [a_simple, b_simple, a, b]; % conditions you wanna use

    % save
    param = {n_rows, n_columns, grid_visual_angle, stim_size, ec_circle, ec_min, fix_stim_dia};
    param_table = cell2table(param, 'VariableNames', {'n_rows', 'n_columns', 'grid_visual_angle', 'stim_size', 'eccentricity_circle', ...
                                                      'eccentricity_min', 'fixation_stimulus_diameter'});
    saveData(param_table, folder_name, 'parameters.csv');
 
    % set the seed for reproducibility
    rng(42);

    %% Function Calls
    % Ask user for input parameters
    data = startGUI();
    
    % Save input parameters in a table
    param = {data.date, data.id, data.name, data.n_trials, data.conditions, data.notes};
    param_table = cell2table(param, 'VariableNames', {'Date', 'ID', 'Name', 'Trials', 'Conditions', 'Notes'});
    saveData(param_table, folder_name, 'info.csv');
 
    %% Generate trials
    trials = generateTrials(data.n_trials, data.conditions, n_rows, n_columns, grid_visual_angle, ec_circle, ec_min);
    trial_data_file_name = fullfile(folder_name, 'trials.mat');
    save(trial_data_file_name, 'trials');
    
    %% Fill trials with angles
    trial_data = createTrialsByCondition(data.n_trials, trials, data.conditions);
    trial_data_file_name = fullfile(folder_name, 'trials_filled.mat');
    save(trial_data_file_name, 'trial_data');

    %% randomize order of trials
    rand_trials = randomize_trials(trial_data, folder_name); % struct to table
    rand_trials_file_name = fullfile(folder_name, 'rand_trials.mat');
    save(rand_trials_file_name, 'rand_trials');

    %% Start Psychtoolbox and display instructions
    %startPsychToolbox(rand_trials, folder_name); 
    
    %%%% Load custom path for testing 
    data_struct = load('C:\Users\idm\Desktop\Semester4\Internship\Matlab\Data\test_20240406_145054\rand_trials.mat');
    data = data_struct.rand_trials;
    folder_name = 'C:\Users\idm\Desktop\Semester4\Internship\Matlab\Data\test_20240406_145054\';
    startPsychToolbox(data, folder_name); 
    %%%
    % Step 6: End experiment and cleanup
    %endExperiment();
end
%%
function startPsychToolbox(data, folder_name)
    subfolder_name = [folder_name 'results' datestr(now, 'yyyymmdd_HHMMSS')];
    mkdir(subfolder_name);
    % Clear the workspace and the screen
    sca;
    close all;
    % Open Psychtoolbox window
    Screen('Preference', 'SkipSyncTests', 1); % Skip synchronization tests, disable for experiment!!
    PsychDefaultSetup(2); % Some default variables

    % Get the screen numbers
    screens = Screen('Screens');
    % Draw to the external screen if avaliable
    screenNumber = max(screens);

    % Define colors for background and stim
    color_bg = WhiteIndex(screenNumber);
    color_stim = BlackIndex(screenNumber);
    
    
    % Calculate center position and screen size
%     screenSize = get(0, 'ScreenSize');
%     screenXpixels = screenSize(3); 
%     screenYpixels = screenSize(4);

    % Open window
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, color_bg); 
    [screenXpixels, screenYpixels] = Screen('WindowSize', window);
    xCenter = screenXpixels / 2;
    yCenter = screenYpixels / 2;
    % Set text properties
    Screen('TextSize', window, 26);
    Screen('TextColor', window, color_stim);
    % Enable alpha blending for anti-aliasing
    Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    %%%%%%% Define texts
    texts = {'Text 1: This is the first text.', ...
             'Text 2: This is the second text.', ...
             'Text 3: This is the third text.'};
    % Display each instruction/text
    for t = 1:numel(texts)
        % Display text
        displayText(window, texts{t});
        % Wait for a key press to continue
        while KbCheck; end % Wait for all keys to be released
        while ~KbCheck; end % Wait for a key press
        [~, ~, keyCode] = KbCheck; % Check for key press
        % Check for Escape key press to close the experiment
        if keyCode(KbName('Escape'))
            break;
        end
    end    

    trial_results = [];

    %%%%%% start experiment
    for i = 1:size(data, 1)
        % Display fixation cross
        %draw_fixation(xCenter, yCenter, color_stim, 40, 4); % colour, size, width

        % Wait 0.3 seconds
        WaitSecs(1);
        %Screen('FillRect', window, color_bg);
        %Screen('Flip', window);

        % Display stimulus
        current_stim = data.TrialMatrix{i};
        current_cond = data.Condition{i};
        current_target_pos = data.TargetSide{i};

        disp(['Current Stim ' data.TargetPosition(i)]);
        disp(['Current Target ' current_target_pos]);

        % Define parameters for the bars
        bar_width = 5; % Width of each bar
        bar_height = 20; % Height of each bar
        jitter = 3; % Maximum amount of jitter for each bar
        % start stim
        trial_result = displayStim(window, bar_width, bar_height, jitter, ...
            current_stim, current_cond, current_target_pos, screenXpixels, screenYpixels, color_stim);
    
        trial_results = [trial_results, trial_result];
        % save/overwrite each loop in case experiment crashes/gets aborted
        results_file_name = fullfile(subfolder_name, 'trial_results.mat');
        save(results_file_name, 'trial_results');

        
        % Check for Escape key press to close the experiment
        if keyCode(KbName('Escape'))
            break;
        end
    end

    disp('Experiment Finished')
    % Close the window
    sca;
end

% function endExperiment()
%     % Cleanup Psychtoolbox and any other resources
% end
