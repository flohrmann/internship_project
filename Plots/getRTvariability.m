function [rt_median_eye,rt_median_button]  = getRTvariability(data, group_labels, ids, unique_conditions, condition_labels, measure, color_map, color_map_individual, comp_results_fix)
% This function computes the variability (std dev or variance) of reaction times (button press and eye arrival time)
% for each participant and condition, then the groups, and then plots the results using line plots with error bars.
% measure: 'std' for standard deviation, 'var' for variance.

numParticipants = length(data); % Number of participants

% [participant x condition]
rtv_eye          = zeros(numParticipants, length(unique_conditions));
rtv_button       = zeros(numParticipants, length(unique_conditions));

rt_median_eye       = zeros(numParticipants, length(unique_conditions));
rt_median_button    = zeros(numParticipants, length(unique_conditions));
rt_mean_eye         = zeros(numParticipants, length(unique_conditions));
rt_mean_button      = zeros(numParticipants, length(unique_conditions));
rt_diff_cond_median = zeros(numParticipants, length(unique_conditions));
rt_diff_cond_mean   = zeros(numParticipants, length(unique_conditions));

% todo make this one work
rt_mu_eye       = zeros(numParticipants, length(unique_conditions));
rt_sigma_eye    = zeros(numParticipants, length(unique_conditions));
rt_tau_eye      = zeros(numParticipants, length(unique_conditions));
rt_mu_button    = zeros(numParticipants, length(unique_conditions));
rt_sigma_button = zeros(numParticipants, length(unique_conditions));
rt_tau_button   = zeros(numParticipants, length(unique_conditions));

% Loop through each participant
for p = 1:numParticipants
    % Loop through each condition
    for c = 1:length(unique_conditions)
        condition = unique_conditions{c};
        condition_trials = strcmp(data(p).Condition, condition);

        rt_eye_cond = data(p).rt_eye(condition_trials);
        rt_button_cond = data(p).rt(condition_trials);
            
        rt_eye_median    = nanmedian(rt_eye_cond);
        rt_button_median = nanmedian(rt_button_cond);
        rt_eye_mean      = nanmean(rt_eye_cond);
        rt_button_mean   = nanmedian(rt_eye_cond);
        
        rt_diff_cond_median(p, c) = nanmedian(rt_button_cond-rt_eye_cond);
        rt_diff_cond_mean(p, c)   = nanmean(rt_button_cond-rt_eye_cond);
      
        % Variability measures
        if strcmp(measure, 'Standard Deviation')
            % Standard Deviation (SD)
            rtv_eye(p, c)          = std(rt_button_cond, 'omitnan');
            rtv_button(p, c)       = std(rt_eye_cond, 'omitnan');
        elseif strcmp(measure, 'Variance')
            % Variance
            rtv_eye(p, c) = var(rt_button_cond, 'omitnan');
            rtv_button(p, c)    = var(rt_eye_cond, 'omitnan');
            
        elseif strcmp(measure, 'Ex-Gaussian Analysis')
            % Ex-Gaussian Analysis (τ, Mu, Sigma)
            % todo read up on what exactly this means/ find other way to do
            % this/check if tau of 0 mathematically makes sense if no tail exists?
            % https://www.sciencedirect.com/science/article/abs/pii/S0891422213003260?via%3Dihub
            if length(rt_button_cond(~isnan(rt_button_cond))) > 3
                pd_button = fitdist(rt_button_cond(~isnan(rt_button_cond)), 'normal');
                mu_button = pd_button.mu;
                sigma_button = pd_button.sigma;
                tau_button = nanmean(rt_button_cond) - mu_button; % Approximate Tau
            else
                mu_button = NaN; sigma_button = NaN; tau_button = NaN;
            end
            
            if length(rt_eye_cond(~isnan(rt_eye_cond))) > 3
                pd_eye = fitdist(rt_eye_cond(~isnan(rt_eye_cond)), 'normal');
                mu_eye = pd_eye.mu;
                sigma_eye = pd_eye.sigma;
                tau_eye = nanmean(rt_eye_cond) - mu_eye; % Approximate Tau
            else
                mu_eye = NaN; sigma_eye = NaN; tau_eye = NaN;
            end
            
            % Store Ex-Gaussian parameters
            rt_mu_eye(p, c)    = mu_eye;
            rt_sigma_eye(p, c) = sigma_eye;
            rt_tau_eye(p, c)   = tau_eye;
            
            rt_mu_button(p, c)    = mu_button;
            rt_sigma_button(p, c) = sigma_button;
            rt_tau_button(p, c)   = tau_button;

        elseif strcmp(measure, 'Mean Squared Successive Difference')
            % Mean Squared Successive Difference (MSSD)
            if length(rt_button_cond(~isnan(rt_button_cond))) > 1
                rt_button_diff = diff(rt_button_cond(~isnan(rt_button_cond)));
                rt_button_var = nanmean(rt_button_diff.^2);
            else
                rt_button_var = NaN;
            end
            
            if length(rt_eye_cond(~isnan(rt_eye_cond))) > 1
                rt_eye_diff = diff(rt_eye_cond(~isnan(rt_eye_cond)));
                rt_eye_var = nanmean(rt_eye_diff.^2);
            else
                rt_eye_var = NaN;
            end
            rtv_eye(p, c)          = rt_eye_var;
            rtv_button(p, c)       = rt_button_var;
        else
            error('Invalid measure. Use "Standard Deviation", "Ex-Gaussian Analysis", or "Mean Squared Successive Difference".');
        end
        rt_median_eye(p, c)    = rt_eye_median;
        rt_median_button(p, c) = rt_button_median;
        rt_mean_eye(p, c)      = rt_eye_mean;
        rt_mean_button(p, c)   = rt_button_mean;
    end
end

       

%% rt median button
adhd_condition_avgs = rt_median_button(strcmp(group_labels, 'ADHD'), :);
nonadhd_condition_avgs = rt_median_button(strcmp(group_labels, 'nonADHD'), :);
adhd_median    = median(adhd_condition_avgs, 1, 'omitnan');
adhd_sem       = std(adhd_condition_avgs, 1, 'omitnan')./ sqrt(size(adhd_condition_avgs, 1));
nonAdhd_median = median(nonadhd_condition_avgs, 1, 'omitnan');
nonAdhd_sem    = std(nonadhd_condition_avgs, 1, 'omitnan')./ sqrt(size(nonadhd_condition_avgs, 1));

plotADHDnonADHDVariance('Median Button Press RT',...
    adhd_condition_avgs, adhd_median, adhd_sem, 'ADHD', 'northeast', ...
    nonadhd_condition_avgs, nonAdhd_median, nonAdhd_sem, 'nonADHD', 'northeast', ...
    ids, condition_labels, 'Median RT (s)', group_labels, unique_conditions, color_map, color_map_individual, ...
    fullfile(comp_results_fix, '11_rt_button_median_allinone_median.png'));

%% rt median eye
adhd_condition_avgs = rt_median_eye(strcmp(group_labels, 'ADHD'), :);
nonadhd_condition_avgs = rt_median_eye(strcmp(group_labels, 'nonADHD'), :);
adhd_median    = median(adhd_condition_avgs, 1, 'omitnan');
adhd_sem       = std(adhd_condition_avgs, 1, 'omitnan')./ sqrt(size(adhd_condition_avgs, 1));
nonAdhd_median = median(nonadhd_condition_avgs, 1, 'omitnan');
nonAdhd_sem    = std(nonadhd_condition_avgs, 1, 'omitnan')./ sqrt(size(nonadhd_condition_avgs, 1));

plotADHDnonADHDVariance('Median Gaze RT',...
    adhd_condition_avgs, adhd_median, adhd_sem, 'ADHD', 'northeast', ...
    nonadhd_condition_avgs, nonAdhd_median, nonAdhd_sem, 'nonADHD', 'northeast', ...
    ids, condition_labels, 'Median RT (s)', group_labels, unique_conditions, color_map, color_map_individual, ...
    fullfile(comp_results_fix, '11_rt_eye_median_allinone_median.png'));


%% rt diff median
adhd_condition_avgs    = rt_diff_cond_median(strcmp(group_labels, 'ADHD'), :);
nonadhd_condition_avgs = rt_diff_cond_median(strcmp(group_labels, 'nonADHD'), :);
adhd_median    = median(adhd_condition_avgs, 1, 'omitnan');
adhd_sem       = std(adhd_condition_avgs, 1, 'omitnan')./ sqrt(size(adhd_condition_avgs, 1));
nonAdhd_median = median(nonadhd_condition_avgs, 1, 'omitnan');
nonAdhd_sem    = std(nonadhd_condition_avgs, 1, 'omitnan')./ sqrt(size(nonadhd_condition_avgs, 1));

plotADHDnonADHDVariance('Median Time between Gaze Arriving at Target and Button Press',...
    adhd_condition_avgs, adhd_median, adhd_sem, 'ADHD', 'northeast', ...
    nonadhd_condition_avgs, nonAdhd_median, nonAdhd_sem, 'nonADHD', 'northeast', ...
    ids, condition_labels, 'time (ms)', group_labels, unique_conditions, color_map, color_map_individual, ...
    fullfile(comp_results_fix, '11_rt_diff_allinone_median.png'));

         
%% rtv button
adhd_condition_avgs    = rtv_button(strcmp(group_labels, 'ADHD'), :);
nonadhd_condition_avgs = rtv_button(strcmp(group_labels, 'nonADHD'), :);
adhd_median    = median(adhd_condition_avgs, 1, 'omitnan');
adhd_sem       = std(adhd_condition_avgs, 1, 'omitnan')./ sqrt(size(adhd_condition_avgs, 1));
nonAdhd_median = median(nonadhd_condition_avgs, 1, 'omitnan');
nonAdhd_sem    = std(nonadhd_condition_avgs, 1, 'omitnan')./ sqrt(size(nonadhd_condition_avgs, 1));

% % of trials gaze actually arrived at target


plotADHDnonADHDVariance(strcat('Button Press RT variability with Median per Condition [', measure, ']'),...
    adhd_condition_avgs, adhd_median, adhd_sem, 'ADHD', 'northeast', ...
    nonadhd_condition_avgs, nonAdhd_median, nonAdhd_sem, 'nonADHD', 'northeast', ...
    ids, condition_labels, 'RTV (ms)', group_labels, unique_conditions, color_map, color_map_individual, ...
    fullfile(comp_results_fix, strcat('11_rt_button_', measure, '_var_allinone_median.png')));

%% rtv eye
adhd_condition_avgs = rtv_eye(strcmp(group_labels, 'ADHD'), :);
nonadhd_condition_avgs = rtv_eye(strcmp(group_labels, 'nonADHD'), :);
adhd_median    = median(adhd_condition_avgs, 1, 'omitnan');
adhd_sem       = std(adhd_condition_avgs, 1, 'omitnan')./ sqrt(size(adhd_condition_avgs, 1));
nonAdhd_median = median(nonadhd_condition_avgs, 1, 'omitnan');
nonAdhd_sem    = std(nonadhd_condition_avgs, 1, 'omitnan')./ sqrt(size(nonadhd_condition_avgs, 1));

plotADHDnonADHDVariance(strcat('Gaze RT variability with Median per Condition [', measure, ']'),...
    adhd_condition_avgs, adhd_median, adhd_sem, 'ADHD', 'northeast', ...
    nonadhd_condition_avgs, nonAdhd_median, nonAdhd_sem, 'nonADHD', 'northeast', ...
    ids, condition_labels, 'RTV (ms)', group_labels, unique_conditions, color_map, color_map_individual, ...
    fullfile(comp_results_fix, strcat('11_rt_eye_', measure, '_var_allinone_median.png')));




% %% rt mean button
% adhd_condition_avgs = rt_mean_button(strcmp(group_labels, 'ADHD'), :);
% nonadhd_condition_avgs = rt_mean_button(strcmp(group_labels, 'nonADHD'), :);
% adhd_median = median(adhd_condition_avgs, 1, 'omitnan');
% adhd_sem = std(adhd_condition_avgs, 0, 1, 'omitnan')./ sqrt(size(adhd_condition_avgs, 1));
% nonAdhd_median = median(nonadhd_condition_avgs, 1, 'omitnan');
% nonAdhd_sem = std(nonadhd_condition_avgs, 0, 1, 'omitnan')./ sqrt(size(nonadhd_condition_avgs, 1));
% 
% plotADHDnonADHDVariance('Mean Button Press RT',...
%     adhd_condition_avgs, adhd_median, adhd_sem, 'ADHD', 'northeast', ...
%     nonadhd_condition_avgs, nonAdhd_median, nonAdhd_sem, 'nonADHD', 'northeast', ...
%     ids, condition_labels, 'Mean RT (s)', group_labels, unique_conditions, color_map, color_map_individual, ...
%     fullfile(comp_results_fix, '11_rt_button_mean_allinone_median.png'));
% 
% 
% %% rt mean eye
% adhd_condition_avgs = rt_mean_eye(strcmp(group_labels, 'ADHD'), :);
% nonadhd_condition_avgs = rt_mean_eye(strcmp(group_labels, 'nonADHD'), :);
% adhd_median = median(adhd_condition_avgs, 1, 'omitnan');
% adhd_sem = std(adhd_condition_avgs, 0, 1, 'omitnan')./ sqrt(size(adhd_condition_avgs, 1));
% nonAdhd_median = median(nonadhd_condition_avgs, 1, 'omitnan');
% nonAdhd_sem = std(nonadhd_condition_avgs, 0, 1, 'omitnan')./ sqrt(size(nonadhd_condition_avgs, 1));
% 
% plotADHDnonADHDVariance('Mean Gaze RT',...
%     adhd_condition_avgs, adhd_median, adhd_sem, 'ADHD', 'northeast', ...
%     nonadhd_condition_avgs, nonAdhd_median, nonAdhd_sem, 'nonADHD', 'northeast', ...
%     ids, condition_labels, 'Mean RT (s)', group_labels, unique_conditions, color_map, color_map_individual, ...
%     fullfile(comp_results_fix, '11_rt_eye_mean_allinone_median.png'));


% Plot RT variability per condition and group with error bars
%     plotRTV_LinePlot(rt_variability, conditions, measure, color_map, comparison_results_folder);
%
%     % Display the RT variability data
%     disp('Reaction Time Variability (Standard Deviation or Variance) by Condition and Group:');
%     disp(rt_variability);
end

function plotRTV_LinePlot(rt_variability, conditions, measure, color_map, comparison_results_folder)
% This function plots the reaction time variability (RTV) per condition and group using line plots with error bars.

% Initialize data for plotting
adhd_button_rtv_mean = zeros(1, length(conditions));
adhd_button_rtv_sem = zeros(1, length(conditions));
adhd_eye_rtv_mean = zeros(1, length(conditions));
adhd_eye_rtv_sem = zeros(1, length(conditions));

nonadhd_button_rtv_mean = zeros(1, length(conditions));
nonadhd_button_rtv_sem = zeros(1, length(conditions));
nonadhd_eye_rtv_mean = zeros(1, length(conditions));
nonadhd_eye_rtv_sem = zeros(1, length(conditions));

% Extract the data for each condition and calculate mean and SEM
for c = 1:length(conditions)
    condition = conditions{c};
    
    % ADHD group
    adhd_button_rtv_mean(c) = mean(rt_variability.(condition).ADHD.RT_ButtonPress_var, 'omitnan');
    adhd_eye_rtv_mean(c) = mean(rt_variability.(condition).ADHD.RT_Eye_var, 'omitnan');
    adhd_button_rtv_sem(c) = std(rt_variability.(condition).ADHD.RT_ButtonPress_var, 'omitnan') / sqrt(length(rt_variability.(condition).ADHD.RT_ButtonPress_var));
    adhd_eye_rtv_sem(c) = std(rt_variability.(condition).ADHD.RT_Eye_var, 'omitnan') / sqrt(length(rt_variability.(condition).ADHD.RT_Eye_var));
    
    % non-ADHD group
    nonadhd_button_rtv_mean(c) = mean(rt_variability.(condition).nonADHD.RT_ButtonPress_var, 'omitnan');
    nonadhd_eye_rtv_mean(c) = mean(rt_variability.(condition).nonADHD.RT_Eye_var, 'omitnan');
    nonadhd_button_rtv_sem(c) = std(rt_variability.(condition).nonADHD.RT_ButtonPress_var, 'omitnan') / sqrt(length(rt_variability.(condition).nonADHD.RT_ButtonPress_var));
    nonadhd_eye_rtv_sem(c) = std(rt_variability.(condition).nonADHD.RT_Eye_var, 'omitnan') / sqrt(length(rt_variability.(condition).nonADHD.RT_Eye_var));
end

% Create the plot
figure;

% Plot RT_ButtonPress variability
subplot(2, 1, 1);
hold on;
errorbar(1:length(conditions), adhd_button_rtv_mean, adhd_button_rtv_sem, '-o', 'Color', color_map('ADHD'), 'DisplayName', 'ADHD');
errorbar(1:length(conditions), nonadhd_button_rtv_mean, nonadhd_button_rtv_sem, '-o', 'Color', color_map('nonADHD'), 'DisplayName', 'non-ADHD');
set(gca, 'XTick', 1:length(conditions), 'XTickLabel', conditions);
xlabel('Condition');
ylabel(strcat('RT ButtonPress Variability: ', measure));
legend('Location', 'best');
title('Reaction Time Variability (ButtonPress) by Condition and Group');
grid on;

% Plot RT_Eye variability
subplot(2, 1, 2);
hold on;
errorbar(1:length(conditions), adhd_eye_rtv_mean, adhd_eye_rtv_sem, '-o', 'Color', color_map('ADHD'), 'DisplayName', 'ADHD');
errorbar(1:length(conditions), nonadhd_eye_rtv_mean, nonadhd_eye_rtv_sem, '-o', 'Color', color_map('nonADHD'), 'DisplayName', 'non-ADHD');
set(gca, 'XTick', 1:length(conditions), 'XTickLabel', conditions);
xlabel('Condition');
ylabel(strcat('RT Eye Variability: ', measure));
legend('Location', 'best');
title('Reaction Time Variability (Eye) by Condition and Group');

hold off;
grid on;

% Adjust the layout
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
% Save the figure if safe is set to 1
%     if safe == 1
%         set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf, fullfile(comparison_results_folder, 'RTV_group_line.png'));
%     end
end