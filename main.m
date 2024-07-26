% start this function to start the experiment

% make sure the DrawFormattedText in this Project is in your path and NOT
% the one from psychtoolbox

% TODOs 
% add tactical stickers to buttons
% extra file with inputs like jitter (everyhting that DOESNT change), make
% main that loads this as input


function main()
% Clear the workspace and the screen
sca; close all; clear;

%% Standard parameters for this experiment
%PsychDefaultSetup(2); % standard setup 
%rand('state',42); % seed for reproducibility
seed = sum(clock); 
rand('state', seed); % seed AGAINST reproducibility

global exp_folder backup_folder ptb_drawformattedtext_oversize
exp_folder = 'C:\Users\flohrmann\Documents\MATLAB\internship_project';
results_folder = 'C:\Users\flohrmann\Documents\Results\';
backup_folder = 'C:\Users\flohrmann\Documents\Backup';

ptb_drawformattedtext_oversize = 2;

% unused
% grid_visual_angle = [34, 46]; % in degrees
% stim_size = [0.12, 1.1]; % in degrees
% ec_circle = 15; % circle of stim pos in degrees
% ec_min = 12; % minimum horizontal eccentricity
% fix_stim_dia = 0.3; % in degrees

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
    'Which Medication used?',...
    'Medication used how many days per week?',...
    'Medication taken today?',...
    'Autism (yes/no)', ...
    'Gaming per week',...
    };
dialog_title='Give_Subject_Information';
num_lines=1;
Subject_default_answer={'','', '', '','','','', '', '', '', '', '', '', '', ''};
%Subject_default_answer={'Fani','1', '1', '','27','normal','normal', 'experimenter', '', '', '', '', '', ''};
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
Subject_ADHD_meds = subject_info{11};
Subject_ADHD_meds_taken = subject_info{12};
Subject_Autism = subject_info{13};
Subject_Gaming = subject_info{14};


%% experiment infos
Exp_prompt={'Lights on?', ...
            'other info', ...
            'Language of Instructions, enter English or German', ...
            'Number of trials for each condition in each block (4 numbers for the 4 conditions, min. 8 each)', ...
            'Duration (in seconds) for time out in a trial without response', ...
            'Viewing distance', ...
            'Display screen width', ...
            'Display screen height', ...
            'Display screen description', ...
            'With eyetracking? (1 yes, 0 no)'};
Exp_dialog_title='Give_Exp_Information';
num_lines=1;
%Exp_default_answer={'Indoor dim light', 'none', 'English', '1', '8 8 8 8', '0 0 0 0', '60', '0', '50 cm', '32 cm', '21 cm',  '', '1'};
Exp_default_answer={'Indoor dim light', '', 'English', '50 50 50 50', '60', '50 cm', '32 cm', '21 cm',  '', '1'};
Exp_info=inputdlg(Exp_prompt,Exp_dialog_title,num_lines,Exp_default_answer);
Exp_RoomLights = Exp_info{1};
Exp_OtherInfo = Exp_info{2};
Exp_InstructionLanguage = Exp_info{3};
if strcmp(Exp_info{3}, 'English') ==1
    EnglishOrGerman = 1;
elseif strcmp(Exp_info{3}, 'German') ==1
    EnglishOrGerman = 2;
else
    R=input('Input language is neither English or German, enter to continue in default English');
    EnglishOrGerman = 1;
end
NTrialsEachCondition = str2num(Exp_info{4});
TimeOut_DurationInSeconds = str2num(Exp_info{5}); % duration in second for time out of no response in a trial.
ViewingDistance = Exp_info{6};
DisplayScreenWidth = Exp_info{7};
DisplayScreenHeight = Exp_info{8};
DisplayScreenDescription = Exp_info{9};
use_eyetracking = str2num(Exp_info{10});  

% if length(NTrialsEachCondition) ~= NConditions |  ...
%         length(NTrialsEachConditionTraining) ~= NConditions
%     'Need as many numbers as the number of experimental conditions';
%     return;
% end


%% Create a folder with User ID and the current date and time
folder_name = strcat(results_folder, [num2str(SubjectID) '_' datestr(now, 'yyyymmdd_HHMMSS')]);
disp(folder_name)
%folder_name = ['C:\Users\flohrmann\Documents\MATLAB\ExperimentRepo\test_' num2str(SubjectID) '_' datestr(now, 'yyyymmdd_HHMMSS')];
mkdir(folder_name);

%% Save all params into autogenerated folder
param = {seed, n_rows, n_columns};%, grid_visual_angle, stim_size, ec_circle, ec_min, fix_stim_dia};
param_table = cell2table(param, 'VariableNames', {'seed', 'n_rows', 'n_columns'});%, 'grid_visual_angle', 'stim_size', 'eccentricity_circle', ...
    %'eccentricity_min', 'fixation_stimulus_diameter'});
saveData(param_table, folder_name, 'parameters.csv');

% Save subject info parameters in a table
infos = {SubjectName, SubjectID, SessionNumber, Subject_Gender, Subject_Age, Subject_LeftEyeSight, Subject_RightEyeSight,...
    Subject_OtherVisionInfo, Subject_ADHD, Subject_ADHD_when,Subject_ADHD_meds,Subject_ADHD_meds_taken, Subject_Autism, Subject_Gaming};
infos_table = cell2table(infos, 'VariableNames', {'SubjectName', 'SubjectID', 'SessionNumber', 'SubjectGender', ...
    'Subject_Age', 'SubjectLeftEyeSight', 'SubjectRightEyeSight', 'SubjectOtherVisionInfo', 'SubjectADHD', 'SubjectADHD_when', 'SubjecADHD_meds', 'SubjectADHD_meds_taken', 'SubjectAutism', 'SubjectGaming'});
saveData(infos_table, folder_name, 'info.csv');

% Save exp infos in table
exps = {Exp_RoomLights, Exp_OtherInfo, Exp_InstructionLanguage, NTrialsEachCondition,  ...
    TimeOut_DurationInSeconds, ViewingDistance, DisplayScreenWidth, DisplayScreenHeight, DisplayScreenDescription};
exps_table = cell2table(exps, 'VariableNames', {'Exp_RoomLights', 'Exp_OtherInfo', 'Exp_InstructionLanguage', 'NTrialsEachCondition', ...
    'TimeOut_DurationInSeconds', 'ViewingDistance', 'DisplayScreenWidth', 'DisplayScreenHeight', 'DisplayScreenDescription'});
saveData(exps_table, folder_name, 'exp_info.csv');

%% Generate trials
n_trials = sum(NTrialsEachCondition);
trials = generateTrials_new(n_trials, n_rows, n_columns);%, grid_visual_angle, ec_circle, ec_min);
trial_data_file_name = fullfile(folder_name, 'trials.mat');
save(trial_data_file_name, 'trials');

%% Fill trials with angles
trial_data = createTrialsByCondition_new(NTrialsEachCondition, trials, conditions);
trial_data_file_name = fullfile(folder_name, 'trials_filled.mat');
save(trial_data_file_name, 'trial_data');

%% randomize order of trials
rand_trials = randomize_trials(trial_data, folder_name); % struct to table
rand_trials_file_name = fullfile(folder_name, 'rand_trials.mat');
save(rand_trials_file_name, 'rand_trials');

%%  Start Psychtoolbox and display instructions
if use_eyetracking == 1 % with eyetracking
    continue_without_eyetracking = false; % flag to check if the experiment should continue without eye-tracking
    [trial_results, samp] = startPsychToolboxEyeTracking(rand_trials, folder_name, n_columns, n_rows, TimeOut_DurationInSeconds, EnglishOrGerman, continue_without_eyetracking);
else % without eyetracking
    continue_without_eyetracking = true; % flag to check if the experiment should continue without eye-tracking
    [trial_results, samp] = startPsychToolboxEyeTracking(rand_trials, folder_name, n_columns, n_rows, TimeOut_DurationInSeconds, EnglishOrGerman, continue_without_eyetracking);
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


end
