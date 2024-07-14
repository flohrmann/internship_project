%% verschub eyetracker and screen how much?

%load('C:\Users\flohrmann\Documents\Results\fani_1_20240701_135200\results\trial_results_20240701_1352.mat')
load('C:\Users\flohrmann\Documents\Results\zhaoping_test\results\trial_results_20240701_1517.mat')

%% based on sampling rate/bildschirmwiederholrate
fs_screen = 59.995; %hz
fs_screen_time = 1/fs_screen; % ms
fs_eye = 60; %hz

period_diff = 1/fs_screen - 1/fs_eye; % time it takes for one whole screen drift

frames4drift = fs_screen_time/ period_diff; % frames it takes to get 1 full verschub

time4drift = frames4drift/fs_screen;

disp(['takes ', num2str(time4drift), ' seconds, ', num2str(time4drift/60), ' minutes or ', num2str(frames4drift),' screen flips, to be one frame behind tracker']);


%% based on eyetracker
% "Total system latency 17 ms/ 1 frame"
% https://www.srlabs.it/wp-content/uploads/2019/10/Tobii_Pro_Nano.pdf
% The translation to computer/system time compensates for drift between the
% two clocks and the communication latency in USB [...]."
% -> use system time

%% from software: screen flip vs. start eyetracker
%trial = 10;
diff_start_list = [];
diff_end_list   = [];
firstSystemTime_list = [];
lastSystemTime_list = [];
firstDeviceTime_list = [];
lastDeviceTime_list = [];
blank_start_time_list = [];
trial_start_time_list = [];
stim_onset_time_list = [];
stim_end_time_list = [];

for trial=1:height(trial_results)
    % eyetracker: first datapoint: flipped one screen earlier -> should be 16 ms diff
    % first and last timestamps for system time and device time
    firstSystemTime = trial_results.eyeTrial(trial).systemTimeStamp(1);
    lastSystemTime  = trial_results.eyeTrial(trial).systemTimeStamp(end);
    firstSystemTime_list = [firstSystemTime_list, firstSystemTime];
    lastSystemTime_list = [lastSystemTime_list, lastSystemTime];
    
    firstDeviceTime = trial_results.eyeTrial(trial).deviceTimeStamp(1);
    lastDeviceTime  = trial_results.eyeTrial(trial).deviceTimeStamp(end);
    firstDeviceTime_list = [firstDeviceTime_list, firstDeviceTime];
    lastDeviceTime_list = [lastDeviceTime_list, lastDeviceTime];
    
    % duration for system time and device time
    durationSystemTime_us = double(lastSystemTime - firstSystemTime); % microseconds
    durationDeviceTime_us = double(lastDeviceTime - firstDeviceTime); % microseconds
    
    % convert durations to ms/s
    durationSystemTime_ms = durationSystemTime_us / 1e3;
    durationDeviceTime_ms = durationDeviceTime_us / 1e3;
    durationSystemTime_s  = durationSystemTime_us / 1e6;
    durationDeviceTime_s  = durationDeviceTime_us / 1e6;
    
    if trial == 1 || trial == 20
        disp(['Duration System Time: ', num2str(durationSystemTime_ms), ' milliseconds or ', num2str(durationSystemTime_s), ' seconds']);
        disp(['Duration Device Time: ', num2str(durationDeviceTime_ms), ' milliseconds or ', num2str(durationSystemTime_s), ' seconds']);
    else
    end
    
    % from screen flip to stim until button press time
    blank_start_time = trial_results.blankStartTime(trial); % seconds
    trial_start_time = trial_results.trialStartTime(trial); % seconds
    stim_onset_time  = trial_results.StimulusOnsetTime(trial); % seconds
    stim_end_time    = trial_results.trialEndTime(trial); % seconds
    blank_start_time_list = [blank_start_time_list, blank_start_time];
    trial_start_time_list = [trial_start_time_list, trial_start_time];
    stim_onset_time_list = [stim_onset_time_list, stim_onset_time];
    stim_end_time_list = [stim_end_time_list, stim_end_time];
    
    duration_stim = stim_end_time - stim_onset_time;
    if trial == 1 || trial == 20
        disp(['Duration Screen Time: ', num2str(duration_stim), ' seconds']);
    else
    end
    % Calculate the time difference between system time and device time
    diff_start = double(firstSystemTime)/1e6 - stim_onset_time;
    diff_end   = double(firstSystemTime)/1e6  - stim_end_time;
    
    % safe
    diff_start_list = [diff_start_list, diff_start];
    diff_end_list = [diff_end_list, diff_end];
    % Display the results
    if trial==1 || trial==20
        disp(['Eyetracker Starts: ', num2str(diff_start), ' s sooner/later than the Stimulation Screen']);
        disp(['Eyetracker Ends: ', num2str(diff_end), ' s sooner/later than the Stimulation Screen']);
    else
    end
end

mean_diff_start = mean(diff_start_list);
min_diff_start = min(diff_start_list);
max_diff_start = max(diff_start_list);

mean_diff_end = mean(diff_end_list);
min_diff_end = min(diff_end_list);
max_diff_end = max(diff_end_list);

% Display the mean differences
disp(['Start Time Differences: Mean:', num2str(mean_diff_start), ' s, Min: ', num2str(min_diff_start), ' s, Max: ', num2str(max_diff_start), ' s']);
disp(['End Time Difference: Mean',     num2str(mean_diff_end), ' s, Min:',    num2str(min_diff_end), ' s, Max: ',   num2str(max_diff_end), ' s']);


%%


% seconds for plotting
firstSystemTime_list = double(firstSystemTime_list') / 1e6;
lastSystemTime_list = double(lastSystemTime_list') / 1e6;
firstDeviceTime_list = double(firstDeviceTime_list') / 1e6;
lastDeviceTime_list = double(lastDeviceTime_list') / 1e6;
blank_start_time_list = double(blank_start_time_list');
trial_start_time_list = double(trial_start_time_list');
stim_onset_time_list = double(stim_onset_time_list');
stim_end_time_list = double(stim_end_time_list');



% Create subplots
figure;
subplot(3, 3, 1);
plot(1:length(firstSystemTime_list), firstSystemTime_list, '-o');
xlabel('Trial Number');
ylabel('Time (s)');
title('First System Time');

subplot(3, 3, 2);
plot(1:length(lastSystemTime_list), lastSystemTime_list, '-o');
xlabel('Trial Number');
ylabel('Time (s)');
title('Last System Time');

subplot(3, 3, 3);
plot(1:length(firstDeviceTime_list), firstDeviceTime_list, '-o');
xlabel('Trial Number');
ylabel('Time (s)');
title('First Device Time');

subplot(3, 3, 4);
plot(1:length(lastDeviceTime_list), lastDeviceTime_list, '-o');
xlabel('Trial Number');
ylabel('Time (s)');
title('Last Device Time');

subplot(3, 3, 5);
plot(1:length(blank_start_time_list), blank_start_time_list, '-o');
xlabel('Trial Number');
ylabel('Time (s)');
title('Blank Start Time');

subplot(3, 3, 6);
plot(1:length(trial_start_time_list), trial_start_time_list, '-o');
xlabel('Trial Number');
ylabel('Time (s)');
title('Trial Start Time');

subplot(3, 3, 7);
plot(1:length(stim_onset_time_list), stim_onset_time_list, '-o');
xlabel('Trial Number');
ylabel('Time (s)');
title('Stim Onset Time');

subplot(3, 3, 8);
plot(1:length(stim_end_time_list), stim_end_time_list, '-o');
xlabel('Trial Number');
ylabel('Time (s)');
title('Stim End Time');

% Adjust layout
sgtitle('Comparison of Times per Trial');
saveas(gcf,'C:\Users\flohrmann\Documents\MATLAB\internship_project\Plots\img\onsets.png')





%% per trial

figure;
hold on;
plot(firstSystemTime_list, '-o', 'DisplayName', 'First System Time');
plot(lastSystemTime_list, '-o', 'DisplayName', 'Last System Time');
%plot(firstDeviceTime_list, '-o', 'DisplayName', 'First Device Time');
%plot(lastDeviceTime_list, '-o', 'DisplayName', 'Last Device Time');
plot(blank_start_time_list, '-o', 'DisplayName', 'Blank Start Time');
plot(trial_start_time_list, '-o', 'DisplayName', 'Trial Start Time');
plot(stim_onset_time_list, '-o', 'DisplayName', 'Stim Onset Time');
plot(stim_end_time_list, '-o', 'DisplayName', 'Stim End Time');
hold off;

% Add labels and legend
xlabel('Trial Number');
ylabel('Time (s)');
title('Comparison of Times per Trial');
legend('show');

%% diff
diff_stim_onset_firstSystem = stim_onset_time_list - firstSystemTime_list;
diff_stim_onset_system = trial_start_time_list - firstSystemTime_list;
diff_end = stim_end_time_list - lastSystemTime_list;
figure;
hold on;
plot(diff_stim_onset_firstSystem, '-o', 'DisplayName', 'StimOnset - FirstSystemTime');
%plot(diff_stim_onset_system, '-o', 'DisplayName', 'trial_start_time_list - First System Time');
plot(diff_end, '-o', 'DisplayName', 'StimEndTime - LastSystemTime');

xlabel('Trial Number');
ylabel('Time Difference (s)');
title('Time Differences Between Stim Onset and System/Trial Start Times');
legend('show');
hold off;

saveas(gcf,'C:\Users\flohrmann\Documents\MATLAB\internship_project\Plots\img\diffs_onset.png')


%%
% genauer als getsecs
% move eytracker to beginning
% splice eyetracker data per trial
% plot trials gaze with stim in background
% plot times/differences
% change questionaire/put invalid option offscreen
% look for doku maybe on yetracker latencies