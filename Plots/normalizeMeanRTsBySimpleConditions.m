function data = normalizeMeanRTsBySimpleConditions(data, comparison_results_folder)
% Function to normalize mean reaction times (RTs) and SEM per condition by the
% mean RT of two specified simple conditions and store the normalized values back into the data struct.
% Inputs:
%   data: struct array containing the data for each observer
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
        mean_rt_a = nanmean(rt(strcmp(condition, 'a')));
        mean_rt_b = nanmean(rt(strcmp(condition, 'b')));
        mean_rt_asimple = nanmean(rt(strcmp(condition, 'a_simple')));
        mean_rt_bsimple = nanmean(rt(strcmp(condition, 'b_simple')));

        sem_rt_a = nanstd(rt(strcmp(condition, 'a'))) / sqrt(sum(~isnan(rt(strcmp(condition, 'a')))));
        sem_rt_b = nanstd(rt(strcmp(condition, 'b'))) / sqrt(sum(~isnan(rt(strcmp(condition, 'b')))));
        sem_rt_asimple = nanstd(rt(strcmp(condition, 'a_simple'))) / sqrt(sum(~isnan(rt(strcmp(condition, 'a_simple')))));
        sem_rt_bsimple = nanstd(rt(strcmp(condition, 'b_simple'))) / sqrt(sum(~isnan(rt(strcmp(condition, 'b_simple')))));

        data(i).RTa = mean_rt_a;
        data(i).RTb = mean_rt_b;
        data(i).RTasimple = mean_rt_asimple;
        data(i).RTbsimple = mean_rt_bsimple;

        data(i).SEMa = sem_rt_a;
        data(i).SEMb = sem_rt_b;
        data(i).SEMasimple = sem_rt_asimple;
        data(i).SEMbsimple = sem_rt_bsimple;
        
        % Calculate the normalization factor based on simple conditions (button press)
        norm_factor = (mean_rt_asimple + mean_rt_bsimple) / 2;

        % Normalize the mean RTs and SEMs (button press)        
        data(i).nRTa = mean_rt_a / norm_factor;
        data(i).nRTb = mean_rt_b / norm_factor;
        data(i).nRTasimple = mean_rt_asimple / norm_factor;
        data(i).nRTbsimple = mean_rt_bsimple / norm_factor;

        data(i).nSEMa = sem_rt_a / norm_factor;
        data(i).nSEMb = sem_rt_b / norm_factor;
        data(i).nSEMasimple = sem_rt_asimple / norm_factor;
        data(i).nSEMbsimple = sem_rt_bsimple / norm_factor;

        %% Calculate mean RTs and SEMs for each condition (eye movement)
        mean_rt_a_eye = nanmean(rt_eye(strcmp(condition, 'a')));
        mean_rt_b_eye = nanmean(rt_eye(strcmp(condition, 'b')));
        mean_rt_asimple_eye = nanmean(rt_eye(strcmp(condition, 'a_simple')));
        mean_rt_bsimple_eye = nanmean(rt_eye(strcmp(condition, 'b_simple')));

        sem_rt_a_eye = nanstd(rt_eye(strcmp(condition, 'a'))) / sqrt(sum(~isnan(rt_eye(strcmp(condition, 'a')))));
        sem_rt_b_eye = nanstd(rt_eye(strcmp(condition, 'b'))) / sqrt(sum(~isnan(rt_eye(strcmp(condition, 'b')))));
        sem_rt_asimple_eye = nanstd(rt_eye(strcmp(condition, 'a_simple'))) / sqrt(sum(~isnan(rt_eye(strcmp(condition, 'a_simple')))));
        sem_rt_bsimple_eye = nanstd(rt_eye(strcmp(condition, 'b_simple'))) / sqrt(sum(~isnan(rt_eye(strcmp(condition, 'b_simple')))));
        
        data(i).RTa_eye = mean_rt_a_eye;
        data(i).RTb_eye = mean_rt_b_eye;
        data(i).RTasimple_eye = mean_rt_asimple_eye;
        data(i).RTbsimple_eye = mean_rt_bsimple_eye;

        data(i).SEMa_eye = sem_rt_a_eye;
        data(i).SEMb_eye = sem_rt_b_eye;
        data(i).SEMasimple_eye = sem_rt_asimple_eye;
        data(i).SEMbsimple_eye = sem_rt_bsimple_eye;

        % Normalize the mean RTs and SEMs (eye movement)
        data(i).nRTa_eye = mean_rt_a_eye / norm_factor;
        data(i).nRTb_eye = mean_rt_b_eye / norm_factor;
        data(i).nRTasimple_eye = mean_rt_asimple_eye / norm_factor;
        data(i).nRTbsimple_eye = mean_rt_bsimple_eye / norm_factor;

        data(i).nSEMa_eye = sem_rt_a_eye / norm_factor;
        data(i).nSEMb_eye = sem_rt_b_eye / norm_factor;
        data(i).nSEMasimple_eye = sem_rt_asimple_eye / norm_factor;
        data(i).nSEMbsimple_eye = sem_rt_bsimple_eye / norm_factor;

        % Store the normalization factor in the data struct
        data(i).norm = norm_factor;
    end

    % Save the normalized data
    save(fullfile(comparison_results_folder, 'data_struct_norm_mean_and_SEM.mat'), 'data');
    disp('Normalization of mean RTs and SEM by simple conditions and saving as struct is complete.');
end
