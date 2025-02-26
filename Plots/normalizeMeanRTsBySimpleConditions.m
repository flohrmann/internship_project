function data = normalizeMeanRTsBySimpleConditions(data, average, unique_conditions, comparison_results_folder)
% Function to normalize mean/median reaction times (RTs) and SEM per condition by the
% mean/median RT of two specified simple conditions and store the normalized values back into the data struct.
% Inputs:
%   data: struct array containing the data for each observer
%   average: string for if/else parts, decides whether mean/median gets calculated
%   comparison_results_folder: string, folder path to save the normalized data
% Output:
%   data: struct array with added fields 'nRTa', 'nRTb', 'nRTasimple', 'nRTbsimple', 'norm',
%         'nSEMa', 'nSEMb', 'nSEMasimple', 'nSEMbsimple',
%         'nRTa_eye', 'nRTb_eye', 'nRTasimple_eye', 'nRTbsimple_eye',
%         'nSEMa_eye', 'nSEMb_eye', 'nSEMasimple_eye', 'nSEMbsimple_eye'

for i = 1:length(data)
    rt        = data(i).rt; % button press RT
    rt_eye    = data(i).rt_eye; % RT of first eye reaching stim
    condition = data(i).Condition;
    
    %% Calculate mean RTs and SEMs for each condition (button press)
    if strcmp(average, 'mean')
        avg_rt_a       = nanmean(rt(strcmp(condition, unique_conditions{1})));
        avg_rt_b       = nanmean(rt(strcmp(condition, unique_conditions{2})));
        avg_rt_asimple = nanmean(rt(strcmp(condition, unique_conditions{3})));
        avg_rt_bsimple = nanmean(rt(strcmp(condition, unique_conditions{4})));
    elseif strcmp(average, 'median')
        avg_rt_a       = nanmedian(rt(strcmp(condition, unique_conditions{1})));
        avg_rt_b       = nanmedian(rt(strcmp(condition, unique_conditions{2})));
        avg_rt_asimple = nanmedian(rt(strcmp(condition, unique_conditions{3})));
        avg_rt_bsimple = nanmedian(rt(strcmp(condition, unique_conditions{4})));
    else
        disp('use mean or median')
    end    
    sem_rt_a       = nanstd(rt(strcmp(condition, unique_conditions{1}))) / sqrt(sum(~isnan(rt(strcmp(condition, unique_conditions{1})))));
    sem_rt_b       = nanstd(rt(strcmp(condition, unique_conditions{2}))) / sqrt(sum(~isnan(rt(strcmp(condition, unique_conditions{2})))));
    sem_rt_asimple = nanstd(rt(strcmp(condition, unique_conditions{3}))) / sqrt(sum(~isnan(rt(strcmp(condition, unique_conditions{3})))));
    sem_rt_bsimple = nanstd(rt(strcmp(condition, unique_conditions{4}))) / sqrt(sum(~isnan(rt(strcmp(condition, unique_conditions{4})))));
    
    data(i).RTa       = avg_rt_a;
    data(i).RTb       = avg_rt_b;
    data(i).RTasimple = avg_rt_asimple;
    data(i).RTbsimple = avg_rt_bsimple;
    
    data(i).SEMa = sem_rt_a;
    data(i).SEMb = sem_rt_b;
    data(i).SEMasimple = sem_rt_asimple;
    data(i).SEMbsimple = sem_rt_bsimple;
    
    % Calculate the normalization factor based on simple conditions (button press)
    norm_factor = (avg_rt_asimple + avg_rt_bsimple) / 2;
    
    % Normalize the mean RTs and SEMs (button press)
    data(i).nRTa       = avg_rt_a / norm_factor;
    data(i).nRTb       = avg_rt_b / norm_factor;
    data(i).nRTasimple = avg_rt_asimple / norm_factor;
    data(i).nRTbsimple = avg_rt_bsimple / norm_factor;
    
    data(i).nSEMa       = sem_rt_a / norm_factor;
    data(i).nSEMb       = sem_rt_b / norm_factor;
    data(i).nSEMasimple = sem_rt_asimple / norm_factor;
    data(i).nSEMbsimple = sem_rt_bsimple / norm_factor;
    
    %% Calculate mean RTs and SEMs for each condition (eye movement)
    if strcmp(average, 'mean')
        avg_rt_a_eye       = nanmean(rt_eye(strcmp(condition, unique_conditions{1})));
        avg_rt_b_eye       = nanmean(rt_eye(strcmp(condition, unique_conditions{2})));
        avg_rt_asimple_eye = nanmean(rt_eye(strcmp(condition, unique_conditions{3})));
        avg_rt_bsimple_eye = nanmean(rt_eye(strcmp(condition, unique_conditions{4})));
    elseif strcmp(average, 'median')
        % yes for id 2 the medians for a and b are the same, no this is not
        % a mistake, their means are different and this routine works you
        % have already checked this twice
        avg_rt_a_eye       = nanmedian(rt_eye(strcmp(condition, unique_conditions{1})));
        avg_rt_b_eye       = nanmedian(rt_eye(strcmp(condition, unique_conditions{2})));
        avg_rt_asimple_eye = nanmedian(rt_eye(strcmp(condition, unique_conditions{3})));
        avg_rt_bsimple_eye = nanmedian(rt_eye(strcmp(condition, unique_conditions{4})));
    else
        disp('use mean or median')
    end
    
    sem_rt_a_eye        = nanstd(rt_eye(strcmp(condition, unique_conditions{1}))) / sqrt(sum(~isnan(rt_eye(strcmp(condition, unique_conditions{1})))));
    sem_rt_b_eye        = nanstd(rt_eye(strcmp(condition, unique_conditions{2}))) / sqrt(sum(~isnan(rt_eye(strcmp(condition, unique_conditions{2})))));
    sem_rt_asimple_eye  = nanstd(rt_eye(strcmp(condition, unique_conditions{3}))) / sqrt(sum(~isnan(rt_eye(strcmp(condition, unique_conditions{3})))));
    sem_rt_bsimple_eye  = nanstd(rt_eye(strcmp(condition, unique_conditions{4}))) / sqrt(sum(~isnan(rt_eye(strcmp(condition, unique_conditions{4})))));
    
    data(i).RTa_eye       = avg_rt_a_eye;
    data(i).RTb_eye       = avg_rt_b_eye;
    data(i).RTasimple_eye = avg_rt_asimple_eye;
    data(i).RTbsimple_eye = avg_rt_bsimple_eye;
    
    data(i).SEMa_eye       = sem_rt_a_eye;
    data(i).SEMb_eye       = sem_rt_b_eye;
    data(i).SEMasimple_eye = sem_rt_asimple_eye;
    data(i).SEMbsimple_eye = sem_rt_bsimple_eye;
    
    % Normalize the mean RTs and SEMs (eye movement)
    data(i).nRTa_eye = avg_rt_a_eye / norm_factor;
    data(i).nRTb_eye = avg_rt_b_eye / norm_factor;
    data(i).nRTasimple_eye = avg_rt_asimple_eye / norm_factor;
    data(i).nRTbsimple_eye = avg_rt_bsimple_eye / norm_factor;
    
    data(i).nSEMa_eye = sem_rt_a_eye / norm_factor;
    data(i).nSEMb_eye = sem_rt_b_eye / norm_factor;
    data(i).nSEMasimple_eye = sem_rt_asimple_eye / norm_factor;
    data(i).nSEMbsimple_eye = sem_rt_bsimple_eye / norm_factor;
    
    % Store the normalization factor in the data struct
    data(i).norm = norm_factor;
end

% Save the normalized data
save(fullfile(comparison_results_folder, strcat('data_struct_norm_', average, '_and_SEM.mat')), 'data');
disp(strcat('Normalization of ', average, ' RTs and SEM by simple conditions and saving as struct is complete.'));
end
