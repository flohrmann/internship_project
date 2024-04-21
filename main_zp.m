% start this function to start the epxeriment
% add cirlce of buttons around button i want pressed

function main()
    %% Clear the workspace and the screen
    sca;
    close all;
    clear;
    PsychDefaultSetup(2);
    % set the seed for reproducibility
    %rng('default') 
    

    global folder_name

    % fani windows desktopauflösung: 1920 x 1080, frame rate 60,033 Hz
    % Get the screen numbers
    screens = Screen('Screens');
    % Draw to the external screen if avaliable
    ScreenNumber = max(screens);

    %set(0,'units','pixels') 
    %Pix_SS = get(0,'screensize')
    [WindowDimensions] = Screen('rect',  ScreenNumber);
    WhiteLuminance = WhiteIndex(ScreenNumber);
    BlackLuminance= BlackIndex(ScreenNumber);
    FrameRate=Screen('NominalFrameRate', ScreenNumber);

    Window_Width = WindowDimensions(3);
    Window_Height = WindowDimensions(4);

    %% Standard parameters for this experiment

    %ZhaopingGuyader_Exp_SetUp()
    
% %         Conditions_A = [29:32];  %condition A in Zhaoping and Guyader 2007
% %         Conditions_B = [33:36];  %condition B in Zhaoping and Guyader 2007
% %         Conditions_Asimple = [49, 50];  %condition Asimple, in Zhaoping and Guyader 2007, repeated twice to match the numbers in A, B etc.
% %         Conditions_Bsimple = [51, 52, 53, 54];  %condition Bsimple, in Zhaoping and Guyader 2007, repeated twice to match the numbers in A, B etc.
% %        
% %         Conditions_ForExp{1} = Conditions_Asimple;
% %         Conditions_ForExp{2} = Conditions_Bsimple;
% %         Conditions_ForExp{3} = Conditions_A;
% %         Conditions_ForExp{4} = Conditions_B;
% %         NConditions = length(Conditions_ForExp);
% %         
        %% Stimuli 
        % All angles are degrees counter clockwise from vertical 
        % Bsp. - = 0° this is the standard drawBar
        %      / = 70°
        %      | = 90°
        %      \ = 110°
        % tuples of angles are crossed bars
        % Bsp. + = [0, 90]
        %      x = [45, 135]
        % In the struct each condition has subconditions that look like:
        %     stimulus, distractor 1, distractor 2
        conditions = struct(); 
        % simple conditions have two equal looking distractors
        conditions.a_simple = [135, 45, 45; 
           % kinda looks like: \    /    / 
                              45, 45, 135];
        conditions.b_simple = [110, 45, 45;
                               70, 135, 135;
                               160, 45, 45;
                               20, 135, 135];
        conditions.a = { % flip distractors to make sure they appear at different positions
            % Bsp 23232323   where 1 is filled with stimulus, and 2, 3 with
            %     32321232   1st and 2nd distractor
            %     23232323
            [135, 90],[45, 0],[45, 90];
            [135, 90],[45, 90],[45, 0]; % distractors flipped
            [135, 0],[45, 90],[45, 0];
            [135, 0],[45, 0],[45, 90]; % distractors flipped
            [45, 90],[135, 0],[135, 90];
            [45, 90],[135, 90],[135, 0]; % distractors flipped
            [45, 0],[135, 90],[135, 0];
            [45, 0],[135, 0],[135, 90]; % distractors flipped
        };
        conditions.b = {
            [110, 90],[70, 0],[70, 90];
            [110, 90],[70, 90],[70, 0]; % distractors flipped
            [160, 0],[70, 90],[70, 0];
            [160, 0],[70, 0],[70, 90]; % distractors flipped
            [70, 90],[110, 0],[120, 90];
            [70, 90],[120, 90],[110, 0]; % distractors flipped
            [20, 0],[110, 90],[110, 0];
            [20, 0],[110, 0],[110, 90]; % distractors flipped
        };
        NConditions = length(fieldnames(conditions));


        %--- set sizes of these NConditions conditions
        %OneSetSize = [12, 9];  %[15, 11];   %--- =[NbX, NbY], the number of columns and rows of the search array
        SetSize_prompt = {'Numbers of columns (minimum 4) and rows (mininum 4) in the search array'}
        SetSize_dialog_title='Row_column_Num';
        num_lines=1;
        SetSize_default_answer={'12 9'};
        SetSize_info = inputdlg(SetSize_prompt,SetSize_dialog_title,num_lines,SetSize_default_answer);
        OneSetSize = str2num(SetSize_info{1}); %
%         % SetSizes_ForExp = repmat(OneSetSize, NConditions, 1);  %== could also give different set sizes for different conditions.
%         %== could also give different set sizes for different conditions.
%         %SetSizes_ForExp = [15, 11; 4, 3; 20, 15; 10, 8];
        n_rows = OneSetSize(2); % number stimulus x axis
        n_columns = OneSetSize(1); % number stimulus y axis

    grid_visual_angle = [34, 46]; % in degrees
    stim_size = [0.12, 1.1]; % in degrees
    ec_circle = 15; % circle of stim pos in degrees
    ec_min = 12; % minimum horizontal eccentricity
    fix_stim_dia = 0.3; % in degrees
    %stimulus_brightness = 48; % in cd/m² -> dont need, white bg & black stim?

    %% Ask user for input parameters
    % Ask user for input/change params
    %data = startGUI();
            % from MainCode_Saved_Feb14_2023
            %---- check resolution, brightness, scale of the display ok or not
            Check_Display_prompt = 'Is display set at brightness 80 percent, size scale 200 percent, and display resolution 3240 x 2160?';
            Resolution_Scale_Brightness_Info = inputdlg(Check_Display_prompt, 'Check Display Setting', 1, {'ok'});
            
            %---- Test audio volumns for feedback
%             fs = 5000; t = 0:0.00002:0.02;
%             LowToneSoundwave =  sin(2*pi*fs/2*t);
%             HighToneSoundwave = sin(2*pi*fs*2*t);
%             VolumnOK =0;
%             while VolumnOK ==0
%                 sound(LowToneSoundwave, fs); pause(0.2);   sound(HighToneSoundwave, fs);
%                 SoundVolumn_prompt = {'Enter 1 or 0 if the volumn of the low and high tune for feedback are ok or otherwise', 'Check audio volumn for feedback use'};
%                 SoundVolumn_prompt_title='Check audio speaker';
%                 num_lines =1;
%                 SoundVolumn_default_answer={'1', 'ok'};
%                 SoundVolumn_Info = inputdlg(SoundVolumn_prompt, SoundVolumn_prompt_title, num_lines, SoundVolumn_default_answer);
%                 VolumnOK = str2num(SoundVolumn_Info{1});
%                 AudioCheckedOK  = SoundVolumn_Info{2};
%             end
            %---- get subject information
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
            
            
            Exp_prompt={'Lights on?', ...
                'other info.', ...
                'Language of Instructions, enter English or German', ...
                'Number of blocks', ...
                'Number of trials for each condition in each block (4 numbers for the 4 conditions)', ...
                'Numbers of Trials for each condition for training (4 numbers for the 4 conditions)', ...
                'Duration (in seconds) for time out in a trial without response', ...
                'Number of breaks in this session = ', ...
                'Viewing distance', 'Display screen width', 'Display screen height', 'Display screen descrption'};
            Exp_dialog_title='Give_Exp_Information';
            num_lines=1;
            Exp_default_answer={'Indoor dim light',  'none', 'German', '1', '8 8 8 8', '0 0 0 0', '60', '0', '50 cm', '32 cm', '21 cm',  'TODO'};
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
                'Need as many numbers as the number of experimental conditions'
                return;
            end
            NumBreaks = str2num(Exp_info{8}); % changed numbers to fit w/o touch
            ViewingDistance = Exp_info{9};
            DisplayScreenWidth = Exp_info{10};
            DisplayScreenHeight = Exp_info{11};
            DisplayScreenDescription = Exp_info{12};
            
            %
            GiveFeedbackOrNot = 0;

            %---- prepare
            KbName('UnifyKeyNames');
            leftarrow=KbName('LeftArrow');
            rightarrow=KbName('RightArrow');
            stopkey=KbName('ESCAPE');
            %HideCursor;
            RandomNumber = sum(100*clock);
            Screen('Preference', 'SkipSyncTests', 1);
            rand('state',RandomNumber);
            
%             day = date; clocktime = clock;
%             fn_out =[num2str(SubjectID), '_', num2str(SessionNumber), '-', day, '-', num2str(clocktime(4)), '-', num2str(clocktime(5)), '.mat'];


    %% Create a folder with User ID and the current date and time 
    folder_name = ['C:\Users\idm\Desktop\Semester4\Internship\Matlab\Data\test_' num2str(SubjectID) '_' datestr(now, 'yyyymmdd_HHMMSS')];
    mkdir(folder_name);
    
    % TODO save new infos
    % Save all params into folder
    %param = {n_rows, n_columns, grid_visual_angle, stim_size, ec_circle, ec_min, fix_stim_dia};
    %param_table = cell2table(param, 'VariableNames', {'n_rows', 'n_columns', 'grid_visual_angle', 'stim_size', 'eccentricity_circle', ...
    %                                                  'eccentricity_min', 'fixation_stimulus_diameter'});
    %saveData(param_table, folder_name, 'parameters.csv');
 
    % Save input parameters in a table
    %param = {data.date, data.id, data.name, data.n_trials, data.conditions, data.notes};
    %param_table = cell2table(param, 'VariableNames', {'Date', 'ID', 'Name', 'Trials', 'Conditions', 'Notes'});
    %saveData(param_table, folder_name, 'info.csv');
 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function Calls
        
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
    %startPsychToolbox(rand_trials, folder_name); 

    %%%% Load custom path for testing 
    data_struct = load('C:\Users\idm\Desktop\Semester4\Internship\Matlab\Data\test_20240406_155217\rand_trials.mat');
    data = data_struct.rand_trials;
    folder_name = 'C:\Users\idm\Desktop\Semester4\Internship\Matlab\Data\test_20240406_155217\';
    startPsychToolbox(data, folder_name); 
    %%%

    % Step 6: Analyse Data
    %endExperiment();

    % optional eye tracking in case it doesnt work
end
%%
function startPsychToolbox(data, folder_name)
    subfolder_name = [folder_name 'results' datestr(now, 'yyyymmdd_HHMMSS')];
    mkdir(subfolder_name);
    
    jitter_x = 0.8*Window_Width/NbX*1/3; % Maximum amount of jitter for each bar
    jitter_y = 0.8*Window_Height/NbY*1/3;


    % Open Psychtoolbox window
    Screen('Preference', 'SkipSyncTests', 1); % Skip synchronization tests, disable for experiment!!

%%%

    % Define colors for background and stim
    color_bg = WhiteIndex(screenNumber);
    color_stim = BlackIndex(screenNumber);
    % Size of the fixation point
    dotSizePix = 8;
    
    % Open window
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, color_bg); 
    [screenXpixels, screenYpixels] = Screen('WindowSize', window);
    xCenter = screenXpixels / 2;
    yCenter = screenYpixels / 2;
    % Set text properties
    Screen('TextSize', window, 30);
    Screen('TextColor', window, color_stim);
    % Enable alpha blending for anti-aliasing
    Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    %%%%%%% Define texts
    texts = {'Welcome :)', ...
             'In this experiment you will be presented a grid of bars \n Press the Left Arrow Key if you can see a uniquely oriented bar on the Left side of the screen \n and the Right Arrow Key if you can see a uniquely oriented bar on the Right side of the screen', ...
             'Before each trial you will see a dot in the middle of the screen \n The next trial will only start if you look at it for around 2 seconds/PRESS A KEY TODO. \n (If you need a break inbetween not looking at it will do the trick ;)', ...
             'Please try to be as fast and as accurate as possible!', ...
             'TODO Insert Exercise Trials', ...
             'Press Any Key To Begin'};
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
    % TODO CHANGE BACK FOR EXPERIMENT 
    % for i = 1:size(data, 1)
    for i = 1:4
        % Display fixation dot & empty screen after
        % TODO still works with keypress, need eyetracking?
        drawFixation(window, xCenter, yCenter, color_stim, dotSizePix);
        drawEmptyScreen(window, color_bg, 1)
       
        % Display stimulus
        current_stim = data.TrialMatrix{i};
        current_cond = data.Condition{i};
        current_target_pos = data.TargetSide{i};

        % Define parameters for the bars
        bar_width = 5; % Width of each bar
        bar_height = 20; % Height of each bar
        


        % start stim
        trial_result = displayStim(window, bar_width, bar_height, jitter_x, jitter_y, ...
            current_stim, current_cond, current_target_pos, screenXpixels, screenYpixels, color_stim);
    
        trial_results = [trial_results, trial_result];
        % save/overwrite each loop in case experiment crashes/gets aborted
        results_file_name = fullfile(subfolder_name, 'trial_results.mat');
        save(results_file_name, 'trial_results');
        
        % Check for Escape key press to close the experiment
        if keyCode(KbName('Escape'))
            disp('Experiment Aborted: Escape Key was pressed')
            break;
        end
    end

    Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    % TODO CHANGE BACK FOR EXPERIMENT 
    %if i == size(data, 1)
    if i == 4 
        % Draw the instructions: in reality the person could press any of
        % the listened to keys to exist. But they do not know that.
        displayText(window, 'Experiment Complete\n Thank you for your participation! <3');
        % Flip to the screen
        %Screen('Flip', window);
        WaitSecs(4);
        % Wait for a key press
        KbStrokeWait(-1);

    end
    disp('Experiment Finished')
    % Close the window
    sca;
end

% function endExperiment()
%     % Cleanup Psychtoolbox and any other resources
% end
