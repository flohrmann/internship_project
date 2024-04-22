% start this function to start the experiment

% TODOs 
% add cirlce of buttons around button i want pressed
% OR add tactical stickers to buttons?
% add german texts and change depending on input
% add breaks in case of serveral blocks/breaks wanted
% add eye tracking
% make eye tracking optional in case it doesnt work
% extra file with inputs like jitter (everyhting that DOESNT change)
% add jitter for each trial in corresponding row/all infos per trial

% strg f for "TODO disable" for real experiment

function main_zp()
% Clear the workspace and the screen
sca; close all; clear;

%% Standard parameters for this experiment
PsychDefaultSetup(2); % standard setup 
rand('state',42); % seed for reproducibility

% unused
grid_visual_angle = [34, 46]; % in degrees
stim_size = [0.12, 1.1]; % in degrees
ec_circle = 15; % circle of stim pos in degrees
ec_min = 12; % minimum horizontal eccentricity
fix_stim_dia = 0.3; % in degrees
%stimulus_brightness = 48; % in cd/m² -> dont need, white bg & black stim?
% screens = Screen('Screens'); % Get the screen numbers
% % Draw to the external screen if avaliable
% ScreenNumber = max(screens);
% [WindowDimensions] = Screen('rect',  ScreenNumber);
%WhiteLuminance = WhiteIndex(ScreenNumber);
%BlackLuminance= BlackIndex(ScreenNumber);
%FrameRate=Screen('NominalFrameRate', ScreenNumber);
% Window_Width = WindowDimensions(3);
% Window_Height = WindowDimensions(4);

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
    45, 45, 135];
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
SetSize_default_answer={'12 9'};
SetSize_info = inputdlg(SetSize_prompt,SetSize_dialog_title,num_lines,SetSize_default_answer);

OneSetSize = str2num(SetSize_info{1});
n_rows = OneSetSize(2); % number stimulus x axis
n_columns = OneSetSize(1); % number stimulus y axis

% old: data = startGUI();
% Prompts from MainCode_Saved_Feb14_2023
%---- check resolution, brightness, scale of the display ok or not
% TODO change numbers depending on display
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
    'gender[lowercase one word]', ...
    'age', ...
    'Left eye sight', ...
    'Right eye sight', ...
    'Other vision info. (e.g., Depth vision)' ...
    };
dialog_title='Give_Subject_Information';
num_lines=1;
Subject_default_answer={'Fani','1', '1', '','27','normal','normal', 'experimenter'};
subject_info=inputdlg(Subject_prompt,dialog_title,num_lines,Subject_default_answer);
SubjectName  = subject_info{1};
SubjectID = str2num(subject_info{2});
SessionNumber = str2num(subject_info{3});
Subject_Gender = subject_info{4};
Subject_Age = str2num(subject_info{5});
Subject_LeftEyeSight = subject_info{6};
Subject_RightEyeSight = subject_info{7};
Subject_OtherVisionInfo = subject_info{8};

%% experiment infos
Exp_prompt={'Lights on?', ...
    'other info.', ...
    'Language of Instructions, enter English or German', ...
    'Number of blocks', ...
    'Number of trials for each condition in each block (4 numbers for the 4 conditions, min. 8 each)', ...
    'Numbers of Trials for each condition for training (4 numbers for the 4 conditions [UNUSED])', ...
    'Duration (in seconds) for time out in a trial without response', ...
    'Number of breaks in this session = ', ...
    'Viewing distance', 'Display screen width', 'Display screen height', 'Display screen descrption'};
Exp_dialog_title='Give_Exp_Information';
num_lines=1;
Exp_default_answer={'Indoor dim light',  'none', 'English', '1', '8 8 8 8', '0 0 0 0', '60', '0', '50 cm', '32 cm', '21 cm',  'TODO'};
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
%TouchOut_XYboundsFraction = str2num(Exp_info{9});   % fraction of screen size in x and y dimension from top-left to touch to terminate a trial

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
folder_name = ['C:\Users\idm\Desktop\Semester4\Internship\Matlab\Data\test_' num2str(SubjectID) '_' datestr(now, 'yyyymmdd_HHMMSS')];
mkdir(folder_name);

%% Save all params into autogenerated folder
param = {n_rows, n_columns, grid_visual_angle, stim_size, ec_circle, ec_min, fix_stim_dia};
param_table = cell2table(param, 'VariableNames', {'n_rows', 'n_columns', 'grid_visual_angle', 'stim_size', 'eccentricity_circle', ...
    'eccentricity_min', 'fixation_stimulus_diameter'});
saveData(param_table, folder_name, 'parameters.csv');

% Save subject info parameters in a table
infos = {SubjectName, SubjectID, SessionNumber, Subject_Gender, Subject_Age, Subject_LeftEyeSight, Subject_RightEyeSight, Subject_OtherVisionInfo};
infos_table = cell2table(infos, 'VariableNames', {'SubjectName', 'SubjectID', 'SessionNumber', 'Subject_Gender', ...
    'Subject_Age', 'Subject_LeftEyeSight', 'Subject_RightEyeSight', 'Subject_OtherVisionInfo'});
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

%% Start Psychtoolbox and display instructions
startPsychToolbox(rand_trials, folder_name, n_columns, n_rows, TimeOut_DurationInSeconds);

%% Load custom path for testing
% data_struct = load('C:\Users\idm\Desktop\Semester4\Internship\Matlab\Data\test_1_20240421_151825\rand_trials.mat');
% data = data_struct.rand_trials;
% folder_name = 'C:\Users\idm\Desktop\Semester4\Internship\Matlab\Data\test_1_20240421_151825';
% startPsychToolbox(data, folder_name, 9, 12, 10);

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

%% Step 6: Analyse Data
end

%% 
function startPsychToolbox(data, folder_name, NbX, NbY, timeout)
% make subfolder for results
subfolder_name = [folder_name '\results'];
mkdir(subfolder_name);
% Define key names
KbName('UnifyKeyNames');
escapeKey = KbName('ESCAPE');
rightKey = KbName('RightArrow');
leftKey = KbName('LeftArrow');

leftShift = KbName('LeftShift');
rightShift = KbName('RightShift');

leftKey = leftShift;
rightKey = rightShift;

% Hide the mouse cursor
HideCursor; % TODO disabled for testing

bar_wh_ratio = 0.08; % bar width to height ratio
screenNumber = 0; % 0 default screen, 1 for external screen
color_bg = WhiteIndex(screenNumber);% colors for background and stim
color_stim = BlackIndex(screenNumber);
dotSizePix = 8; % Size of fixation point
text_size = 30;

% Open Psychtoolbox window
PsychDebugWindowConfiguration(0,0.5); % TODO disable, makes window transparent
Screen('Preference', 'SkipSyncTests', 1); % TODO disable, skips synchronization tests
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, color_bg);
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
xCenter = screenXpixels / 2;
yCenter = screenYpixels / 2;

% Maximum amount of jitter for each bar
jitter_x = round(0.1*screenXpixels/NbX);
jitter_y = round(0.1*screenYpixels/NbY);

% save infos
infos = {color_bg, color_stim, text_size, dotSizePix, jitter_x, jitter_y, xCenter, yCenter};
infos_table = cell2table(infos, 'VariableNames', {'color_bg', 'color_stim', ...
     'text_size', 'dotSizePix', 'jitter_x', 'jitter_y', 'xCenter', 'yCenter'});
saveData(infos_table, folder_name, 'screen_info.csv');

% Set text properties
Screen('TextSize', window, text_size);
Screen('TextColor', window, color_stim);
% Enable alpha blending for anti-aliasing
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%% Instructions
texts = {'Welcome :)', ...
    'In this experiment you will be presented a grid of bars \n\n Press the Left Arrow Key if you can see a uniquely oriented bar on the Left side of the screen \n and the Right Arrow Key if you can see a uniquely oriented bar on the Right side of the screen', ...
    'Lets look at a possible images you will see later', ...
    'The trials can look like this for example. Can you find the uniquely oriented bar?',...
    'The trials can also look like this. Can you still see the uniquely oriented bar?', ...
    'Before each trial you will see a fixation dot in the middle of the screen, please remember to look at it for as long as it is there \n\n The next trial will start automatically after half a second/if you look at it for around 2 seconds.', ...
    'During the experiment once you have found the uniquely oriented bar remember to press the corresponding button \n\n Please try to be as fast and as accurate as possible!', ...
    'Do you have any questions? \n \n If you feel ready: Press the Right Arrow Key To Begin The Experiment'};
t = 1;
while t <= numel(texts)
    % Display text or image based on the current index
    if t == 4 % Add images to certain 'slides'
        img = imread('C:\Users\idm\Desktop\Semester4\Internship\Matlab\ExperimentRepo\Images\example_1.png'); 
        displayImage(window, texts{t}, img, screenXpixels, screenYpixels);
    elseif t == 5
        img = imread('C:\Users\idm\Desktop\Semester4\Internship\Matlab\ExperimentRepo\Images\example_2.png'); 
        displayImage(window, texts{t}, img, screenXpixels, screenYpixels);
    else
        displayText(window, texts{t}); % Display text only
    end
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

%% Start experiment
disp('Start Experiment');
results_file_name = [subfolder_name, '\trial_results_' , datestr(now, 'yyyymmdd_HHMM'), '.mat'];
trial_results = [];
% TODO disable/change back
% for i = 1:size(data, 1)
for i = 1:10
    % Display fixation dot & empty screen after
    drawFixation(window, xCenter, yCenter, color_stim, dotSizePix);
    WaitSecs(0.5); % show for 0.5 seconds
    drawEmptyScreen(window, color_bg, 1)
    WaitSecs(0.2); % wait before drawing stim

    current_stim = data.AngleMatrix{i};
    current_cond = data.Condition{i};
    current_target_pos = data.TargetSide{i};

    % Display stimulus
    trial_result = displayStim_new(window, bar_wh_ratio, jitter_x, jitter_y, current_stim,...
        current_cond, current_target_pos, screenXpixels, screenYpixels, color_stim, timeout);

    trial_results = [trial_results, trial_result];
    % save/overwrite each loop in case experiment crashes/gets aborted
    save(results_file_name, 'trial_results');

    % Check for Escape key press to close the experiment
    if keyCode('ESCAPE')
        disp('Experiment Aborted: Escape Key was pressed')
        break;
    end
end

%% End screen
% TODO CHANGE BACK FOR EXPERIMENT
%if i == size(data, 1)
if i == 10
    displayText(window, 'Experiment Complete\n\n Thank you for your participation! <3');
    WaitSecs(2);
    KbStrokeWait(-1);  % Wait to end with key press
end
disp('Experiment Finished')
sca; % Close the window
end