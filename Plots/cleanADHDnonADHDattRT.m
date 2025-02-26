function ttestResults = cleanADHDnonADHDattRT(merged_data, conditions)
% Initialize empty arrays to hold mean RTs and accuracies
mean_rt_adhd = [];
mean_accuracy_adhd = [];
mean_rt_nonadhd = [];
mean_rt_a_adhd = [];
mean_rt_b_adhd = [];
mean_rt_asimple_adhd = [];
mean_rt_bsimple_adhd = [];
mean_rt_a_nonadhd = [];
mean_rt_b_nonadhd = [];
mean_rt_asimple_nonadhd = [];
mean_rt_bsimple_nonadhd = [];
mean_accuracy_nonadhd = [];
mean_6inatt_adhd = [];
count_6inatt_adhd = [];
mean_allinatt_adhd = [];
mean_inatt_nonadhd = [];
count_6inatt_nonadhd = [];
mean_allinatt_nonadhd = [];
mean_con_adhd = [];
mean_con_nonadhd = [];
rtv_a_adhd = [];
rtv_b_adhd = [];
rtv_asimple_adhd = [];
rtv_bsimple_adhd = [];
rtv_a_nonadhd = [];
rtv_b_nonadhd = [];
rtv_asimple_nonadhd = [];
rtv_bsimple_nonadhd = [];

inattCategorical = [];

% Loop through each participant and calculate mean RT and accuracy
for i = 1:size(merged_data,1)
    % Check if the participant belongs to the ADHD group
    if strcmp(merged_data.group_merged_data{i}, 'ADHD')
        % the mean RT
        mean_rt_a_adhd = [mean_rt_a_adhd; merged_data.nRTa(i)];
        mean_rt_b_adhd = [mean_rt_b_adhd; merged_data.nRTb(i)];
        mean_rt_asimple_adhd = [mean_rt_asimple_adhd; merged_data.nRTasimple(i)];
        mean_rt_bsimple_adhd = [mean_rt_bsimple_adhd; merged_data.nRTbsimple(i)];
        % rtv
        rtv_a_adhd = [rtv_a_adhd; merged_data.a(i).RT_ButtonPress_var];
        rtv_b_adhd = [rtv_b_adhd; merged_data.b(i).RT_ButtonPress_var];
        rtv_asimple_adhd = [rtv_asimple_adhd; merged_data.a_simple(i).RT_ButtonPress_var];
        rtv_bsimple_adhd = [rtv_bsimple_adhd; merged_data.b_simple(i).RT_ButtonPress_var];
        % mean accuracy
        mean_accuracy_adhd = [mean_accuracy_adhd; mean(merged_data.accuracy{i}, 'omitnan')];
        % confusion
        mean_con_adhd  = [mean_con_adhd; merged_data.nRTa(i)/merged_data.nRTb(i)];
        % inattention level
        mean_6inatt_adhd  = [mean_6inatt_adhd; merged_data.PartAMeanSymptoms(i)];
        count_6inatt_adhd = [count_6inatt_adhd; merged_data.PartASymptomsCount(i)];
        mean_allinatt_adhd = [mean_allinatt_adhd; merged_data.MeanSymptoms(i)];      
        
        inattCategorical  = [inattCategorical; merged_data.InattentionLevel(i)];      
        
    % Check if the participant belongs to the non-ADHD group
    elseif strcmp(merged_data.group_merged_data{i}, 'nonADHD')
        % mean RT
        %rt_values = [merged_data.nRTa(i), merged_data.nRTb(i), merged_data.nRTasimple(i), merged_data.nRTbsimple(i)];
        %mean_rt_nonadhd = [mean_rt_nonadhd; rt_values]; % mean(rt_values, 'omitnan')];
        mean_rt_a_nonadhd = [mean_rt_a_nonadhd; merged_data.nRTa(i)];
        mean_rt_b_nonadhd = [mean_rt_b_nonadhd; merged_data.nRTb(i)];
        mean_rt_asimple_nonadhd = [mean_rt_asimple_nonadhd; merged_data.nRTasimple(i)];
        mean_rt_bsimple_nonadhd = [mean_rt_bsimple_nonadhd; merged_data.nRTbsimple(i)];
        % rtv
        rtv_a_nonadhd = [rtv_a_nonadhd; merged_data.a(i).RT_ButtonPress_var];
        rtv_b_nonadhd = [rtv_b_nonadhd; merged_data.b(i).RT_ButtonPress_var];
        rtv_asimple_nonadhd = [rtv_asimple_nonadhd; merged_data.a_simple(i).RT_ButtonPress_var];
        rtv_bsimple_nonadhd = [rtv_bsimple_nonadhd; merged_data.b_simple(i).RT_ButtonPress_var];
        % confusion
        mean_con_nonadhd  = [mean_con_nonadhd; merged_data.nRTa(i)/merged_data.nRTb(i)];
        % mean accuracy
        mean_accuracy_nonadhd = [mean_accuracy_nonadhd; mean(merged_data.accuracy{i}, 'omitnan')];
        % inattention level
        mean_inatt_nonadhd = [mean_inatt_nonadhd; merged_data.PartAMeanSymptoms(i)];
        count_6inatt_nonadhd = [count_6inatt_nonadhd; merged_data.PartASymptomsCount(i)];
        mean_allinatt_nonadhd = [mean_allinatt_nonadhd; merged_data.MeanSymptoms(i)];
    end
end




% %% ttest 
ttestResults = runTTests_ADHD_vs_nonADHD(...
    mean_rt_a_adhd, mean_rt_b_adhd, ...
    mean_rt_asimple_adhd, mean_rt_bsimple_adhd, ...
    mean_rt_a_nonadhd, mean_rt_b_nonadhd, ...
    mean_rt_asimple_nonadhd, mean_rt_bsimple_nonadhd, ...
    rtv_a_adhd, rtv_a_nonadhd,...
    rtv_b_adhd,rtv_b_nonadhd, ... 
    rtv_asimple_adhd,rtv_asimple_nonadhd,... 
    rtv_bsimple_adhd,rtv_bsimple_nonadhd, ...
    mean_accuracy_adhd, mean_accuracy_nonadhd, ...
    mean_con_adhd, mean_con_nonadhd, ...
    mean_6inatt_adhd, mean_inatt_nonadhd, ...
    count_6inatt_adhd, count_6inatt_nonadhd, ...
    mean_allinatt_adhd, mean_allinatt_nonadhd);



%% Correlation Analysis
% Plot both correlations

% Correlation Analysis for BOTH group
[r_adhd_rt, p_adhd_rt]   = corr(mean_6inatt_adhd, [mean_rt_a_adhd, mean_rt_b_adhd, mean_rt_asimple_adhd, mean_rt_bsimple_adhd]);
[r_adhd_rtv, p_adhd_rtv] = corr(mean_6inatt_adhd, [rtv_a_adhd, rtv_b_adhd, rtv_asimple_adhd, rtv_bsimple_adhd]);
[r_adhd_con, p_adhd_con] = corr(mean_6inatt_adhd, mean_con_adhd);
[r_adhd_acc, p_adhd_acc] = corr(mean_6inatt_adhd, mean_accuracy_adhd);
 
% % Display results for BOTH group
for i = 1:length(r_adhd_rt)
    fprintf('%s Correlation between inattention levels and RT: r = %.2f, p = %.3f\n', conditions{i}, r_adhd_rt(1, i), p_adhd_rt(1, i));
    fprintf('%s Correlation between inattention levels and RTV: r = %.2f, p = %.3f\n', conditions{i}, r_adhd_rtv(1, i), p_adhd_rtv(1, i));
end
fprintf('Correlation between inattention levels and Accuracy: r = %.2f, p = %.3f\n', r_adhd_acc, p_adhd_acc);
fprintf('Correlation between inattention levels and Confusion: r = %.2f, p = %.3f\n', r_adhd_con, p_adhd_con);



%% both groups in one, correlation for categorical inatt

% Initialize empty arrays to hold mean RTs and accuracies
mean_rt_adhd = [];
mean_accuracy_adhd = [];
mean_rt_a_adhd = [];
mean_rt_b_adhd = [];
mean_rt_asimple_adhd = [];
mean_rt_bsimple_adhd = [];
mean_6inatt_adhd = [];
count_6inatt_adhd = [];
mean_allinatt_adhd = [];
mean_con_adhd = [];
rtv_a_adhd = [];
rtv_b_adhd = [];
rtv_asimple_adhd = [];
rtv_bsimple_adhd = [];
inattCategorical = [];

% Loop through each participant and calculate mean RT and accuracy
for i = 1:size(merged_data,1)
        % the mean RT
        mean_rt_a_adhd = [mean_rt_a_adhd; merged_data.nRTa(i)];
        mean_rt_b_adhd = [mean_rt_b_adhd; merged_data.nRTb(i)];
        mean_rt_asimple_adhd = [mean_rt_asimple_adhd; merged_data.nRTasimple(i)];
        mean_rt_bsimple_adhd = [mean_rt_bsimple_adhd; merged_data.nRTbsimple(i)];
        % rtv
        rtv_a_adhd = [rtv_a_adhd; merged_data.a(i).RT_ButtonPress_var];
        rtv_b_adhd = [rtv_b_adhd; merged_data.b(i).RT_ButtonPress_var];
        rtv_asimple_adhd = [rtv_asimple_adhd; merged_data.a_simple(i).RT_ButtonPress_var];
        rtv_bsimple_adhd = [rtv_bsimple_adhd; merged_data.b_simple(i).RT_ButtonPress_var];
        % mean accuracy
        mean_accuracy_adhd = [mean_accuracy_adhd; mean(merged_data.accuracy{i}, 'omitnan')];
        % confusion
        mean_con_adhd  = [mean_con_adhd; merged_data.nRTa(i)/merged_data.nRTb(i)];
        % inattention level
        mean_6inatt_adhd  = [mean_6inatt_adhd; merged_data.PartAMeanSymptoms(i)];
        count_6inatt_adhd = [count_6inatt_adhd; merged_data.PartASymptomsCount(i)];
        mean_allinatt_adhd = [mean_allinatt_adhd; merged_data.MeanSymptoms(i)];      
        inattCategorical  = [inattCategorical; merged_data.InattentionLevel(i)];      
end


% same but for categorical inattention values
% Correlation Analysis for BOTH group
[r_adhd_rt, p_adhd_rt]   = corr(inattCategorical, [mean_rt_a_adhd, mean_rt_b_adhd, mean_rt_asimple_adhd, mean_rt_bsimple_adhd]);
[r_adhd_rtv, p_adhd_rtv] = corr(inattCategorical, [rtv_a_adhd, rtv_b_adhd, rtv_asimple_adhd, rtv_bsimple_adhd]);
[r_adhd_con, p_adhd_con] = corr(inattCategorical, mean_con_adhd);
[r_adhd_acc, p_adhd_acc] = corr(inattCategorical, mean_accuracy_adhd);
 
% % Display results for BOTH group
for i = 1:length(r_adhd_rt)
    fprintf('%s Correlation between categorical inattention levels and RT: r = %.2f, p = %.3f\n', conditions{i}, r_adhd_rt(1, i), p_adhd_rt(1, i));
    fprintf('%s Correlation between categorical inattention levels and RTV: r = %.2f, p = %.3f\n', conditions{i}, r_adhd_rtv(1, i), p_adhd_rtv(1, i));
end
fprintf('Correlation between categorical inattention levels and Accuracy: r = %.2f, p = %.3f\n', r_adhd_acc, p_adhd_acc);
fprintf('Correlation between categorical inattention levels and Confusion: r = %.2f, p = %.3f\n', r_adhd_con, p_adhd_con);



