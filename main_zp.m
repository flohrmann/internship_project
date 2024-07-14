% start this function to start the experiment

% make sure the DrawFormattedText in this Project is in your path and NOT
% the one from psychtoolbox

% strg f for "TODO disable" for real experiment

% TODOs 
% add tactical stickers to buttons
% extra file with inputs like jitter (everyhting that DOESNT change), make
% main that loads this as input


function main_zp()
% Clear the workspace and the screen
sca; close all; clear;

%% Standard parameters for this experiment
%PsychDefaultSetup(2); % standard setup 
%rand('state',42); % seed for reproducibility
rand('state', sum(clock)); % seed AGAINST reproducibility

global exp_folder backup_folder ptb_drawformattedtext_oversize
exp_folder = 'C:\Users\flohrmann\Documents\MATLAB\internship_project';
results_folder = 'C:\Users\flohrmann\Documents\Results';
backup_folder = 'C:\Users\flohrmann\Documents\Backup';

ptb_drawformattedtext_oversize = 2;

% unused
grid_visual_angle = [34, 46]; % in degrees
stim_size = [0.12, 1.1]; % in degrees
ec_circle = 15; % circle of stim pos in degrees
ec_min = 12; % minimum horizontal eccentricity
fix_stim_dia = 0.3; % in degrees

%% Stimuli
% All angles are degrees counter clockwise from vertical
% Bsp. - = 0° this is the standard drawBar
%      / = 70°   | = 90°   \ = 110°
% tuples of angles are crossed bars
% Bsp. + = [0, 90]     x = [45, 135]
% In the struct each condition has subconditions that look like:
% [stimulus, distractor 1, distractor 2]
conditions = struct();
% simple conditions have two equal looking distractors
conditions.a_simple = [135, 45, 45;
    % kinda looks like: \   /   /
    45, 135, 135];
conditions.b_simple = [110, 45, 45;
    70, 135, 135;
    160, 45, 45;
    20, 135, 135];
% normal conditions have two different looking distractors
conditions.a = { % flip distractors to make sure they appear at each position
    [135, 90],[45, 0],[45, 90];
    [135, 90],[45, 90],[45, 0]; % distractors flipped
    [135, 0],[45, 90],[45, 0];
    [135, 0],[45, 0],[45, 90];  % distractors flipped
    [45, 90],[135, 0],[135, 90];
    [45, 90],[135, 90],[135, 0];% distractors flipped
    [45, 0],[135, 90],[135, 0];
    [45, 0],[135, 0],[135, 90]; % distractors flipped
    };
conditions.b = {
    [110, 90],[45, 0],[45, 90];
    [110, 90],[45, 90],[45, 0]; % distractors flipped
    [160, 0],[45, 90],[45, 0];
    [160, 0],[45, 0],[45, 90];  % distractors flipped
    [70, 90],[135, 0],[135, 90];
    [70, 90],[135, 90],[135, 0];% distractors flipped
    [20, 0],[135, 90],[135, 0];
    [20, 0],[135, 0],[135, 90]; % distractors flipped
    };

NConditions = length(fieldnames(conditions));

%% User Input 
SetSize_prompt = {'Numbers of columns (minimum 4) and rows (mininum 4) in the search array'};
SetSize_dialog_title='Row_column_Num';
num_lines=1;
SetSize_default_answer={'12 9'}; % this 
%etSize_default_answer={'14 11'};

SetSize_info = inputdlg(SetSize_prompt,SetSize_dialog_title,num_lines,SetSize_default_answer);

OneSetSize = str2num(SetSize_info{1});
n_rows = OneSetSize(2); % number stimulus x axis
n_columns = OneSetSize(1); % number stimulus y axis

%---- check resolution, brightness, scale of the display ok or not
Check_Display_prompt = 'Is display set at brightness 80 percent, size scale 200 percent, and display resolution 3240 x 2160?';
Resolution_Scale_Brightness_Info = inputdlg(Check_Display_prompt, 'Check Display Setting', 1, {'ok'});

%---- Test audio volumns for feedback
fs = 5000; t = 0:0.00002:0.02;
LowToneSoundwave =  sin(2*pi*fs/2*t);
HighToneSoundwave = sin(2*pi*fs*2*t);
VolumnOK =0;
while VolumnOK ==0
    sound(LowToneSoundwave, fs); pause(0.5);   sound(HighToneSoundwave, fs); pause(0.5);
    SoundVolumn_prompt = {'Enter 1 or 0 if the volumn of the low and high tune for feedback are ok or otherwise', 'Check audio volumn for feedback use'};
    SoundVolumn_prompt_title='Check audio speaker';
    num_lines =1;
    SoundVolumn_default_answer={'1', 'ok'};
    SoundVolumn_Info = inputdlg(SoundVolumn_prompt, SoundVolumn_prompt_title, num_lines, SoundVolumn_default_answer);
    VolumnOK = str2num(SoundVolumn_Info{1});
    AudioCheckedOK  = SoundVolumn_Info{2};
end

%% subject information
Subject_prompt={'Subject name (no space)', ...
    'Subject ID', ...
    'Session Number (for this subject)', ...
    'Gender[lowercase one word]', ...
    'Age', ...
    'Left eye sight', ...
    'Right eye sight', ...
    'Other vision info. (e.g., Depth vision)' ...
    'ADHD (yes/no)',...
    'If ADHD, diagnosis when?',...
    'Autism (yes/no)', ...
    'Gaming per week',...
    };
dialog_title='Give_Subject_Information';
num_lines=1;
Subject_default_answer={'Fani','1', '1', '','27','normal','normal', 'experimenter', '', '', '', ''};
subject_info=inputdlg(Subject_prompt,dialog_title,num_lines,Subject_default_answer);
SubjectName  = subject_info{1};
SubjectID = str2num(subject_info{2});
SessionNumber = str2num(subject_info{3});
Subject_Gender = subject_info{4};
Subject_Age = str2num(subject_info{5});
Subject_LeftEyeSight = subject_info{6};
Subject_RightEyeSight = subject_info{7};
Subject_OtherVisionInfo = subject_info{8};
Subject_ADHD = subject_info{9};
Subject_ADHD_when = subject_info{10};
Subject_Autism = subject_info{11};
Subject_Gaming = subject_info{12};


%% experiment infos
Exp_prompt={'Lights on?', 'other info.', ...
    'Language of Instructions, enter English or German', ...
    'Number of blocks', ...
    'Number of trials for each condition in each block (4 numbers for the 4 conditions, min. 8 each)', ...
    'Numbers of Trials for each condition for training (4 numbers for the 4 conditions [UNUSED])', ...
    'Duration (in seconds) for time out in a trial without response', ...
    'Number of breaks in this session = ', 'Viewing distance', 'Display screen width', ...
    'Display screen height', 'Display screen description', 'With eyetracking? (1 yes, 0 no)'};
Exp_dialog_title='Give_Exp_Information';
num_lines=1;
%Exp_default_answer={'Indoor dim light', 'none', 'English', '1', '8 8 8 8', '0 0 0 0', '60', '0', '50 cm', '32 cm', '21 cm',  '', '1'};
Exp_default_answer={'Indoor dim light', 'none', 'English', '1', '40 40 40 40', '0 0 0 0', '60', '0', '50 cm', '32 cm', '21 cm',  '', '1'};

%Exp_default_answer={'Done', 'Indoor dim light', 'yes/irrelevant', 'none', '30 30 30 30', '1 1 1 1', '0', '4', '+', '1', '40 cm', '40.9 cm', '25.6 cm', 'Attwood lab display screen'};
Exp_info=inputdlg(Exp_prompt,Exp_dialog_title,num_lines,Exp_default_answer);
Exp_RoomLights = Exp_info{1};
Exp_OtherInfo = Exp_info{2};
Exp_InstructionLanguage = Exp_info{3};
if strcmp(Exp_info{3}, 'English') ==1
    EnglishOrGerman = 1;
elseif strcmp(Exp_info{3}, 'German') ==1
    EnglishOrGerman = 2;
else
    R=input('Input language is neither English or German, enter to continue in default German');
    EnglishOrGerman = 2;
end
NumberOfBlocks= str2num(Exp_info{4});
NTrialsEachCondition = str2num(Exp_info{5});
NTrialsEachConditionTraining = str2num(Exp_info{6});
%TouchToleranceFraction = str2num(Exp_info{7});   % a number <1 for the tolerance of touch location as a fraction of the window dimension.
TimeOut_DurationInSeconds = str2num(Exp_info{7});   % duration in second for time out of no response in a trial.
use_eyetracking = str2num(Exp_info{13});   % fraction of screen size in x and y dimension from top-left to touch to terminate a trial

if length(NTrialsEachCondition) ~= NConditions |  ...
        length(NTrialsEachConditionTraining) ~= NConditions
    'Need as many numbers as the number of experimental conditions';
    return;
end
NumBreaks = str2num(Exp_info{8}); % changed numbers to fit w/o touch
ViewingDistance = Exp_info{9};
DisplayScreenWidth = Exp_info{10};
DisplayScreenHeight = Exp_info{11};
DisplayScreenDescription = Exp_info{12};

%% Create a folder with User ID and the current date and time
folder_name = strcat(results_folder, ['\fani_' num2str(SubjectID) '_' datestr(now, 'yyyymmdd_HHMMSS')]);
disp(folder_name)
%folder_name = ['C:\Users\flohrmann\Documents\MATLAB\ExperimentRepo\test_' num2str(SubjectID) '_' datestr(now, 'yyyymmdd_HHMMSS')];
mkdir(folder_name);

%% Save all params into autogenerated folder
param = {n_rows, n_columns, grid_visual_angle, stim_size, ec_circle, ec_min, fix_stim_dia};
param_table = cell2table(param, 'VariableNames', {'n_rows', 'n_columns', 'grid_visual_angle', 'stim_size', 'eccentricity_circle', ...
    'eccentricity_min', 'fixation_stimulus_diameter'});
saveData(param_table, folder_name, 'parameters.csv');

% Save subject info parameters in a table
infos = {SubjectName, SubjectID, SessionNumber, Subject_Gender, Subject_Age, Subject_LeftEyeSight, Subject_RightEyeSight, Subject_OtherVisionInfo, Subject_ADHD, Subject_ADHD_when, Subject_Autism, Subject_Gaming};
infos_table = cell2table(infos, 'VariableNames', {'SubjectName', 'SubjectID', 'SessionNumber', 'Subject_Gender', ...
    'Subject_Age', 'Subject_LeftEyeSight', 'Subject_RightEyeSight', 'Subject_OtherVisionInfo', 'Subject_ADHD', 'Subject_ADHD_when', 'Subject_Autism', 'Subject_Gaming'});
saveData(infos_table, folder_name, 'info.csv');

% Save exp infos in table
exps = {Exp_RoomLights, Exp_OtherInfo, Exp_InstructionLanguage, NumberOfBlocks, NTrialsEachCondition, NTrialsEachConditionTraining,  ...
    TimeOut_DurationInSeconds, NumBreaks, ViewingDistance, DisplayScreenWidth, DisplayScreenHeight, DisplayScreenDescription};
exps_table = cell2table(exps, 'VariableNames', {'Exp_RoomLights', 'Exp_OtherInfo', 'Exp_InstructionLanguage', 'NumberOfBlocks', 'NTrialsEachCondition', 'NTrialsEachConditionTraining', ...
    'TimeOut_DurationInSeconds', 'NumBreaks', 'ViewingDistance', 'DisplayScreenWidth', 'DisplayScreenHeight', 'DisplayScreenDescription'});
saveData(exps_table, folder_name, 'exp_info.csv');

%% Generate trials
n_trials = NumberOfBlocks * sum(NTrialsEachCondition) + sum(NTrialsEachConditionTraining);
trials = generateTrials_new(n_trials, n_rows, n_columns, grid_visual_angle, ec_circle, ec_min);
trial_data_file_name = fullfile(folder_name, 'trials.mat');
save(trial_data_file_name, 'trials');

%% Fill trials with angles
trial_data = createTrialsByCondition_new(NumberOfBlocks, NTrialsEachCondition, trials, conditions);
trial_data_file_name = fullfile(folder_name, 'trials_filled.mat');
save(trial_data_file_name, 'trial_data');

%% randomize order of trials
rand_trials = randomize_trials(trial_data, folder_name); % struct to table
rand_trials_file_name = fullfile(folder_name, 'rand_trials.mat');
save(rand_trials_file_name, 'rand_trials');

%%  Start Psychtoolbox and display instructions
if use_eyetracking == 1 % with eyetracking
    results = startPsychToolboxEyeTracking(rand_trials, folder_name, n_columns, n_rows, TimeOut_DurationInSeconds, EnglishOrGerman);
else % without eyetracking
    startPsychToolbox(rand_trials, folder_name, n_columns, n_rows, TimeOut_DurationInSeconds, EnglishOrGerman);
end

%% Clean up Fixation Data and Safe
%clean_data = aggregateFixationData(results);
%results.sampFixAll = table2array(clean_data);
%results_file_name = [folder_name, 'results\trial_results_fixed.mat'];
%save(results_file_name, 'results');

%% Ask for user feedback
Ending_prompt={'Type in subject report', ...
    'write down experimenter comment/observations'};
dialog_title='Ending_Information';
num_lines=1;
ending_default_answer={'none','none'};
Ending_info=inputdlg(Ending_prompt, dialog_title,num_lines,ending_default_answer);
Subject_EndingReport  = Ending_info{1};
Experimenter_Comments = Ending_info{2};
end_table = cell2table({Subject_EndingReport, Experimenter_Comments}, 'VariableNames', {'Subject_EndingReport', 'Experimenter_Comments'});
saveData(end_table, folder_name, 'end_comments.csv');

%% User questionaire of (in)attention
createQuestionnaire(folder_name, num2str(SubjectID));

%% Step 6: Analyse Data
plotConditionSpreadAndStimPosition(rand_trials, n_rows, n_columns);



end

function trial_results = startPsychToolboxEyeTracking(data, folder_name, NbX, NbY, timeout, language)
global exp_folder
global ptb_drawformattedtext_oversize

% make subfolder for results
subfolder_name = [folder_name '\results'];
disp(subfolder_name)
mkdir(subfolder_name);
trial_results = table();

% Define key names
KbName('UnifyKeyNames');
escapeKey = KbName('ESCAPE');
rightKey = KbName('RightArrow');
leftKey = KbName('LeftArrow');
leftShift = KbName('LeftShift');
rightShift = KbName('RightShift');
leftKey = leftShift;
rightKey = rightShift;
spaceKey = KbName('space');
deleteKey = KbName('Delete'); % ENTF key to stop eyetracking

continue_without_eyetracking = false; % flag to check if the experiment should continue without eye-tracking

% Hide the mouse cursor
%HideCursor; % TODO disable, for testing

bar_wh_ratio = 0.08; % bar width to height ratio
screenNumber = 0; % 0 default screen, 1 for external screen
color_bg = WhiteIndex(screenNumber);% colors for background and stim
color_stim = BlackIndex(screenNumber);
dotSizePix = 15; % Size of fixation point
text_size = 40;

% Open Psychtoolbox window
%PsychDebugWindowConfiguration(0,0.5); % TODO disable, makes window transparent
Screen('Preference', 'SkipSyncTests', 1); % TODO disable, skips synchronization tests
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, color_bg);
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
xCenter = screenXpixels / 2;
yCenter = screenYpixels / 2;

% Maximum amount of jitter for each bar
jitter_x = round(0.1*screenXpixels/NbX);
jitter_y = round(0.1*screenYpixels/NbY);

% save infos
infos = {bar_wh_ratio, color_bg, color_stim, text_size, dotSizePix, jitter_x, jitter_y, xCenter, yCenter};
infos_table = cell2table(infos, 'VariableNames', {'bar_wh_ratio', 'color_bg', 'color_stim', ...
     'text_size', 'dotSizePix', 'jitter_x', 'jitter_y', 'xCenter', 'yCenter'});
saveData(infos_table, folder_name, 'screen_info.csv');
% TODO safe
monitorFlipInterval = Screen('GetFlipInterval', window); % to get inter flip interval for break

    
% Set text properties
Screen('TextSize', window, text_size);
Screen('TextColor', window, color_stim);
% Enable alpha blending for anti-aliasing
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


%% Instructions Part 1 
if language == 1 % choose instructions language
    text_path = strcat(exp_folder, '\Images\instructions_en');
else
    text_path = strcat(exp_folder, '\Images\instructions_ger');
end


% Show the first 2 pictures/instructions
t = 1;
while t <= 2
    % Display first two images based on the current index
    img = imread(strcat(text_path, sprintf('\\instructions-0%0d.jpg', t)));
    displayImage(window, img, screenXpixels, screenYpixels);
    while 1 % Check for key presses to navigate through texts
        [~, ~, keyCode] = KbCheck;
        if keyCode(rightKey) % Check for right arrow press
            if t == 2 % Check if it's the second image
                WaitSecs(0.2); % Prevent rapid key presses
                t = t + 1; % Move index beyond the second to exit the loop
                break; % Exit the inner while loop
            else
                t = t + 1; % Otherwise, just move to the next image
                WaitSecs(0.2);
                break;
            end
        elseif keyCode(leftKey) && t > 1
            t = t - 1; % Navigate to the previous image
            WaitSecs(0.2);
            break;
        elseif keyCode(escapeKey)
            disp('Experiment Aborted During First Instructions');
            sca; % Close the window
            return; % Close function   
        end
    end
    if t > 2
        break; % Break the outer loop if past the second image
    end
end


%% Initialize Titta and set up the eye tracker 
try
    settings = Titta.getDefaults('Tobii Pro Nano');
    eye_tracker = Titta(settings);
    eye_tracker.init();
catch
    error('No eye trackers found. Make sure the eye tracker is connected and try again.');
end

% calibrate & safe results
try
    calibration = eye_tracker.calibrate(window);
    calibration_file = [subfolder_name, '\calibration_results_' , datestr(now, 'yyyymmdd_HHMM'), '.mat'];
    save(calibration_file, 'calibration');
catch
    %error('Eye tracker couldnt be calibrated. Try again or restart experiment without eye tracking.');
    % Calibration failed, allow user to continue without eye-tracking  
    Screen('TextSize', window, text_size);
    displayText(window, 'Calibration failed. Press ENTF to continue without eye-tracking or the Spacebar to try again.');
    %Screen('Flip', window);
    WaitSecs(1);
    while 1
        [~, ~, keyCode] = KbCheck;
        if keyCode(deleteKey)
            continue_without_eyetracking = true;
            break;
        elseif keyCode(spaceKey)
            try
                settings = Titta.getDefaults('Tobii Pro Nano');
                eye_tracker = Titta(settings);
                eye_tracker.init();
                calibration = eye_tracker.calibrate(window);
                calibration_file = [subfolder_name, '\calibration_results_' , datestr(now, 'yyyymmdd_HHMM'), '.mat'];
                save(calibration_file, 'calibration');
            catch
                Screen('TextSize', window, text_size);
                displayText(window, 'Calibration failed. Press ENTF to continue without eye-tracking.');
            end
            break;
        end
    end
end

%% Instructions Part 2
image_number = 3; % Start from the third image
while image_number <= 16 % until picture 16
    try % einstellige nr
        img = imread(strcat(text_path, sprintf('\\instructions-0%0d.jpg', image_number)));
    catch
        try % zweistellige nr
            img = imread(strcat(text_path, sprintf('\\instructions-%02d.jpg', image_number)));
        catch
            break; % Exit the loop if no more images are found
        end
        %break; % Exit the loop if no more images are found
    end
    displayImage(window, img, screenXpixels, screenYpixels);
    while 1 % Check for key presses to navigate through images
        [~, ~, keyCode] = KbCheck;
        if keyCode(rightKey) % Check for right arrow press
            image_number = image_number + 1; % Move to the next image
            WaitSecs(0.2); % Prevent rapid key presses
            break;
        elseif keyCode(leftKey) && image_number > 3
            image_number = image_number - 1; % Navigate to the previous image
            WaitSecs(0.2);
            break;
        elseif keyCode(escapeKey)
            disp('Experiment Aborted During Instructions');
            sca; % Close the window
            return; % Close function         
        end
    end
end


%% Experiment Loop
disp('Start Experiment');
results_file_name = [subfolder_name, '\trial_results_' , datestr(now, 'yyyymmdd_HHMM'), '.mat'];

%try
    for i = 1:size(data, 1)
        current_data = data(i,:);
        %for i = 1:4
        %     % Display fixation dot & empty screen after
        %     fixStartTime = GetSecs;
        %     % Define the area around the fixation point where gaze must fall
        %     fixationRadius = 50; % pixels, defines a circle radius around the fixation point
        %     minFixationTime = 1.0; % seconds, the minimum time gaze must remain within the radius
        %     samp_fix = drawFixationAndWait(window, xCenter, yCenter, color_stim, dotSizePix, eye_tracker, screenXpixels, screenYpixels, fixationRadius, minFixationTime);
        %         %drawFixation(window, xCenter, yCenter, color_stim, dotSizePix); %without eyetrackgin
        %         %WaitSecs(0.5);
        %     fixEndTime = GetSecs;
        %     drawEmptyScreen(window, color_bg, 1)
        %     WaitSecs(0.2); % wait before drawing stim
        % %     end
        if ~continue_without_eyetracking
            eye_tracker.buffer.start('gaze');
        end
        % Display Trial Start Instructions until button press
        Screen('FillRect', window, color_bg);
        Screen('TextSize', window, text_size);
        DrawFormattedText(window, 'Put your two index fingers on the shift buttons to get ready \n\n\n and press either of them to continue with the next trial', 'center', 'center');
        trialStartTime = Screen('Flip', window);
        %displayText(window, 'Put your two index fingers on the shift buttons to get ready \n\n\n and press either of them to continue with the next trial');
        while true
            [~, ~, keyCode] = KbCheck;
            if keyCode(leftKey) || keyCode(rightKey)
                break;
            elseif keyCode(deleteKey)
                disp('Experiment switched to no eyetracking: Entf/Delete Key was pressed')
                continue_without_eyetracking = true;
                eye_tracker.buffer.stop('gaze');
                eye_tracker.deInit();
                break;
            elseif keyCode(escapeKey)
                disp('Experiment Aborted: Escape Key was pressed')
                eye_tracker.buffer.stop('gaze');
                eye_tracker.deInit(); % stop eye tracker
                sca; % close window
                return; % close function
            end
        end
        % Display fixation dot for 0.5 seconds
        Screen('DrawDots', window, [xCenter; yCenter], dotSizePix, color_stim, [], 2);
        fixStartTime = Screen('Flip', window); % next possible screen flip
        %drawFixation(window, xCenter, yCenter, color_stim, dotSizePix); % fixationcross without eyetrackgin
        
        % Display empty screen for 0.2 seconds (wait 0.5 for screen flip)
        Screen('FillRect', window, color_bg); % blank screen
        blankStartTime = Screen('Flip', window, fixStartTime+0.5-0.5*monitorFlipInterval); % wait 0.5 s
        
        if ~continue_without_eyetracking % eyetracking data during instr/fixation
            samp_fix = eye_tracker.buffer.consumeN('gaze'); % struct with eyetracking data
        else
            samp_fix = struct('deviceTimeStamp', [], 'systemTimeStamp', [], 'left', struct(), 'right', struct()); % No eye-tracking data/empty scruct
        end

        % Display stimulus (wait 0.2 s for screen flip)
        [trial_result, stop] = displayStim_new(window, current_data, bar_wh_ratio, jitter_x, jitter_y, screenXpixels, screenYpixels,...
            color_stim, timeout, blankStartTime, monitorFlipInterval, continue_without_eyetracking, eye_tracker, color_bg);
        
        % Check for Escape key press to close the experiment
        if stop == 1
            disp('Experiment Aborted: Escape Key was pressed')
            break;
        end
        %      % Get eye data, since start/last getting eye data & safe
        %     if ~continue_without_eyetracking
        %         % Get eye data, since start/last getting eye data & save
        %         samp = eye_tracker.buffer.consumeN('gaze');
        %     else
        %         samp = struct('deviceTimeStamp', [], 'systemTimeStamp', [], 'left', struct(), 'right', struct()); % No eye-tracking data
        %     end
        
        % Add data from fixation
        trial_result.trialStartTime = trialStartTime;
        trial_result.fixStartTime = fixStartTime;
        trial_result.fixTime = fixStartTime - trialStartTime;
        trial_result.eyeFix = samp_fix; % eye fixation
        
        % save/overwrite each loop in case experiment crashes/gets aborted
        trial_results(i,:) = trial_result(1,:);
        save(results_file_name, 'trial_results');
        
        clear trial_result
    end
% catch
%     disp('Experiment Crashed')
% end
% Clean up: stop eye tracker, save results, and close
if ~continue_without_eyetracking
    eye_tracker.buffer.stop('gaze');
    eye_tracker.deInit();
end

%% End
% End screen
% TODO CHANGE BACK FOR EXPERIMENT
if i == size(data, 1)
%if i == 4
    Screen('FillRect', window, color_bg);
    Screen('TextSize', window, text_size);
    if language == 2 % choose instructions language     
        displayText(window, 'Experiment Abgeschlossen\n\n Vielen Dank fürs Mitmachen!');
    else
        img = imread(strcat(text_path, '\instructions-18.jpg'));
        displayImage(window, img, screenXpixels, screenYpixels);    
    end
    WaitSecs(2);
    KbStrokeWait(-1);  % Wait to end with key press
end

disp('Experiment Finished')
sca; % Close the window
end
