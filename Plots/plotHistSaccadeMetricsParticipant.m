function plotHistSaccadeMetricsParticipant(trial_metrics, cutData, conditions,color_map, sr, safename_1, safename_2)

    %% All-in-one metrics
    figure;
    t = tiledlayout(2, 2, 'TileSpacing', 'Compact');
    num_sacc = trial_metrics.n_sac;
    nexttile; histogram(num_sacc, 'FaceColor', 'g', 'BinWidth', 1);
    title('Total Saccades per Trial'); xlabel('Number'); ylabel('Frequency');

    mean_dist_sac = trial_metrics.mean_dist_sac;
    nexttile; histogram(mean_dist_sac, 'FaceColor', 'g', 'BinWidth', 50);
    title('Mean Distance Saccades per Trial'); xlabel('Distance (mm)'); ylabel('Frequency');

    first_tss = cellfun(@(x) x(1), trial_metrics.ts_start);  % Extract the first number
    first_tss_time = first_tss .* sr;  % Scale by sampling rate
    nexttile; histogram(first_tss_time, 'FaceColor', 'g', 'BinWidth', 0.02);
    title('Start Time First Saccade'); xlabel('Time from Stimulation Onset (s)'); ylabel('Frequency');

    sacc_after_first_target = cellfun(@(x) x(1), trial_metrics.saccades_after_target);
    nexttile; histogram(sacc_after_first_target, 'FaceColor', 'g', 'BinWidth', 1);
    title('Saccades After Target Found'); xlabel('Number'); ylabel('Frequency');

    % Save the figure
    saveas(gcf, [safename_1, '.svg']);  % Save as SVG
    saveas(gcf, [safename_1, '.png']);  % Save as PNG

%     %% Metrics per condition
%     % Initialize variables to store global min/max for each metric
%     num_sacc_limits = [inf, -inf];
%     mean_dist_limits = [inf, -inf];
%     first_tss_time_limits = [inf, -inf];
%     sacc_after_target_limits = [inf, -inf];
% 
%     % Step 1: Calculate global axis limits for all conditions
%     for i = 1:length(conditions)
%         condition = conditions{i};  % Current condition
%         condition_trials_idx = strcmp(cutData.Condition, condition); % Filter trials for the current condition
% 
%         % Metrics for the current condition
%         num_sacc_condition = trial_metrics.n_sac(condition_trials_idx);
%         mean_dist_sac_condition = trial_metrics.mean_dist_sac(condition_trials_idx);
%         first_tss_condition = cellfun(@(x) x(1), trial_metrics.ts_start(condition_trials_idx));
%         first_tss_time_condition = first_tss_condition .* sr;  % Scale by sampling rate
%         sacc_after_target_condition = cellfun(@(x) x(1), trial_metrics.saccades_after_target(condition_trials_idx));
% 
%         % Update global min/max for each metric
%         num_sacc_limits = [min(num_sacc_limits(1), min(num_sacc_condition)), ...
%             max(num_sacc_limits(2), max(num_sacc_condition))];
% 
%         mean_dist_limits = [min(mean_dist_limits(1), min(mean_dist_sac_condition)), ...
%             max(mean_dist_limits(2), max(mean_dist_sac_condition))];
% 
%         first_tss_time_limits = [min(first_tss_time_limits(1), min(first_tss_time_condition)), ...
%             max(first_tss_time_limits(2), max(first_tss_time_condition))];
% 
%         sacc_after_target_limits = [min(sacc_after_target_limits(1), min(sacc_after_target_condition)), ...
%             max(sacc_after_target_limits(2), max(sacc_after_target_condition))];
%     end
% 
%     % Step 2: Plot histograms with consistent axis limits
%     figure;
%     tiledlayout(length(conditions), 4, 'TileSpacing', 'Compact');  % Layout for all conditions
% 
%     for i = 1:length(conditions)
%         condition = conditions{i};  % Current condition
%         condition_trials_idx = strcmp(cutData.Condition, condition); % Filter trials for the current condition
% 
%         % Metrics for the current condition
%         num_sacc_condition = trial_metrics.n_sac(condition_trials_idx);
%         mean_dist_sac_condition = trial_metrics.mean_dist_sac(condition_trials_idx);
%         first_tss_condition = cellfun(@(x) x(1), trial_metrics.ts_start(condition_trials_idx));
%         first_tss_time_condition = first_tss_condition .* sr;  % Scale by sampling rate
%         sacc_after_target_condition = cellfun(@(x) x(1), trial_metrics.saccades_after_target(condition_trials_idx));
% 
%         % Plot histograms for each metric
%         % Saccades per trial
%         nexttile;
%         histogram(num_sacc_condition, 'FaceColor', color_map(condition), 'BinWidth', 0.5);
%         title(sprintf('Saccades per Trial (%s)', condition));
%         xlabel('Number');
%         ylabel('Frequency');
%         xlim(num_sacc_limits);  % Apply consistent x-axis limits
%         ylim([0, inf]);         % Ensure consistent y-axis
% 
%         % Mean distance saccades per trial
%         nexttile;
%         histogram(mean_dist_sac_condition, 'FaceColor', color_map(condition), 'BinWidth', 50);
%         title(sprintf('Mean Distance Saccades (%s)', condition));
%         xlabel('Distance (mm)');
%         ylabel('Frequency');
%         xlim(mean_dist_limits);  % Apply consistent x-axis limits
%         ylim([0, inf]);          % Ensure consistent y-axis
% 
%         % Start time of first saccade
%         nexttile;
%         histogram(first_tss_time_condition, 'FaceColor', color_map(condition), 'BinWidth', 0.02);
%         title(sprintf('Start Time First Saccade (%s)', condition));
%         xlabel('Time from Stimulation Onset (s)');
%         ylabel('Frequency');
%         xlim(first_tss_time_limits);  % Apply consistent x-axis limits
%         ylim([0, inf]);               % Ensure consistent y-axis
% 
%         % Saccades after the target found
%         nexttile;
%         histogram(sacc_after_target_condition, 'FaceColor', color_map(condition), 'BinWidth', 0.5);
%         title(sprintf('Saccades After Target (%s)', condition));
%         xlabel('Number');
%         ylabel('Frequency');
%         xlim(sacc_after_target_limits);  % Apply consistent x-axis limits
%         ylim([0, inf]);                  % Ensure consistent y-axis
%     end
% 
%     % Save the figure
%     saveas(gcf, [safename_2, '.svg']);  % Save as SVG
%     saveas(gcf, [safename_2, '.png']);  % Save as PNG
end
