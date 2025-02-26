function rts = getAndPlotRTV_ParticipantLines(data, conditions,condition_labels, measure, groups, color_map, comparison_results_folder)
% This function computes the variability (std dev or variance) of reaction times (RTs)
% for each condition and group, and then plots the results using individual participant lines.
% Additionally, it computes and stores the mean, median, and SEM of RTs.
% measure: 'std' for standard deviation, 'var' for variance.

% Initialize a struct to store the results
rts = struct();

% Loop through each condition
for c = 1:length(conditions)
    condition = conditions{c};
    
    % Initialize arrays to store RT variability for each group
    rt_var_ADHD_button = [];
    rt_var_ADHD_eye = [];
    rt_var_nonADHD_button = [];
    rt_var_nonADHD_eye = [];
    
    % Loop through each participant
    for participant = 1:length(data)
        if ismember(condition, data(participant).Condition)
            % Get the group for this participant (ADHD or non-ADHD)
            group = data(participant).group;
            
            % Find the trials that match the current condition
            condition_trials = strcmp(data(participant).Condition, condition);
            rt_button_all = data(participant).rt(condition_trials);
            rt_eye_all = data(participant).rt_eye(condition_trials);
            
            % Calculate variability (std or var) for RT_ButtonPress and RT_Eye
            if strcmp(measure, 'std')
                rt_button_var = std(rt_button_all, 'omitnan');
                rt_eye_var = std(rt_eye_all, 'omitnan');
            elseif strcmp(measure, 'var')
                rt_button_var = var(rt_button_all, 'omitnan');
                rt_eye_var = var(rt_eye_all, 'omitnan');
            else
                error('Invalid measure. Use "std" for standard deviation or "var" for variance.');
            end
            
            % Calculate mean, median, and SEM for RT_ButtonPress and RT_Eye
            rt_button_mean = mean(rt_button_all, 'omitnan');
            rt_eye_mean = mean(rt_eye_all, 'omitnan');
            rt_button_median = median(rt_button_all, 'omitnan');
            rt_eye_median = median(rt_eye_all, 'omitnan');
            rt_button_sem = std(rt_button_all, 'omitnan') / sqrt(sum(~isnan(rt_button_all)));
            rt_eye_sem = std(rt_eye_all, 'omitnan') / sqrt(sum(~isnan(rt_eye_all)));
            
            % Store the RT variability for the respective group and condition
            if strcmp(group, 'ADHD')
                rt_var_ADHD_button = [rt_var_ADHD_button; rt_button_var];
                rt_var_ADHD_eye = [rt_var_ADHD_eye; rt_eye_var];
            elseif strcmp(group, 'nonADHD')
                rt_var_nonADHD_button = [rt_var_nonADHD_button; rt_button_var];
                rt_var_nonADHD_eye = [rt_var_nonADHD_eye; rt_eye_var];
            end
            
            % Store all computed values per participant
            rts(participant).id = data(participant).id;
            rts(participant).group = group;
            rts(participant).(condition).RT_ButtonPress_var = rt_button_var;
            rts(participant).(condition).RT_Eye_var = rt_eye_var;
            rts(participant).(condition).RT_ButtonPress_mean = rt_button_mean;
            rts(participant).(condition).RT_Eye_mean = rt_eye_mean;
            rts(participant).(condition).RT_ButtonPress_median = rt_button_median;
            rts(participant).(condition).RT_Eye_median = rt_eye_median;
            rts(participant).(condition).RT_ButtonPress_sem = rt_button_sem;
            rts(participant).(condition).RT_Eye_sem = rt_eye_sem;         
        end
        % mean rt normalized by easy conditions
        rts(participant).a.RT_Eye_norm           = data(participant).nRTa;
        rts(participant).a.RT_Button_norm        = data(participant).nRTa_eye;
        rts(participant).a_simple.RT_Eye_norm    = data(participant).nRTasimple;
        rts(participant).a_simple.RT_Button_norm = data(participant).nRTasimple_eye;
        rts(participant).b.RT_Eye_norm           = data(participant).nRTb;
        rts(participant).b.RT_Button_norm        = data(participant).nRTb_eye;
        rts(participant).b_simple.RT_Eye_norm    = data(participant).nRTbsimple;
        rts(participant).b_simple.RT_Button_norm = data(participant).nRTbsimple_eye;
    end
end


    
%% RT Visualization: 
num_conditions = length(conditions);
jitter_amount = 0.15; % Reduced jitter for better spacing
% Separate ADHD and non-ADHD indices
groups = {rts.group};
adhd_indices = strcmp(groups, 'ADHD');
nonadhd_indices = strcmp(groups, 'nonADHD');

%% Mean and Median: Eye RT - Button RT - Diff [Log Y-Axis]
set(gcf, 'Position', [100, 100, 1000, 600]); % Resize the figure window (x, y, width, height)
    % Subplot 1: Mean Eye RT
    subplot(2, 3, 1); title('Mean Eye RT');
    [me_min, me_max] = plotRTGroupComparison(rts, conditions, condition_labels, groups, color_map, 'RT_Eye_mean', 'Mean EyeRT log(ms)');

    % Subplot 2: Mean Button RT
    subplot(2, 3, 2); title('Mean Button RT');
    [mb_min, mb_max] = plotRTGroupComparison(rts, conditions, condition_labels, groups, color_map, 'RT_ButtonPress_mean', 'Mean ButtonRT log(ms)');

    % Subplot 3: Mean Button RT - Mean Eye RT
    subplot(2, 3, 3); title('Mean RT Difference');
    [md_min, md_max, adhd_handle, nonadhd_handle] = plotRTDifference(rts, conditions, condition_labels, groups, color_map, 'RT_ButtonPress_mean', 'RT_Eye_mean', 'ButtonRT-EyeRT log(ms)');
    lg = legend([adhd_handle, nonadhd_handle], {'ADHD', 'nonADHD'}, 'Location', 'northeast');
    lg.Box = 'off';
    % Subplot 4: Median Eye RT
    subplot(2, 3, 4); title('Median Eye RT');
    [ke_min, ke_max] = plotRTGroupComparison(rts, conditions, condition_labels, groups, color_map, 'RT_Eye_median', 'Median Eye RT log(ms)');

    % Subplot 5: Median Button RT
    subplot(2, 3, 5); title('Median Button RT');
    [kb_min, kb_max] = plotRTGroupComparison(rts, conditions, condition_labels, groups, color_map, 'RT_ButtonPress_median', 'Median Button RT log(ms)');

    % Subplot 6: Median Button RT - Mean Eye RT
    subplot(2, 3, 6); title('Median RT Difference');
    [kd_min, kd_max, ~, ~] = plotRTDifference(rts, conditions, condition_labels, groups, color_map, 'RT_ButtonPress_median', 'RT_Eye_mean', 'ButtonRT-EyeRT log(ms)');
    
    % Set consistent y-limits for each column
    % eye
    new_lims = [min([me_min, ke_min, mb_min, kb_min]), max([me_max, ke_max, mb_max, kb_max])];
    subplot(2, 3, 1); ylim(new_lims);
    subplot(2, 3, 4); ylim(new_lims);
    % button
    subplot(2, 3, 2); ylim(new_lims);
    subplot(2, 3, 5); ylim(new_lims);
    % diff
    subplot(2, 3, 3); ylim([min(md_min, kd_min), max(md_max, kd_max)]);
    subplot(2, 3, 6); ylim([min(md_min, kd_min), max(md_max, kd_max)]);
    safe = 1;
    % Resize and Save the Figure
    if safe == 1
        set(gcf, 'Position', [100, 100, 1000, 600]); % Resize the figure window (x, y, width, height)
        saveas(gcf, strcat(comparison_results_folder, '\03_rt_individual_group_mean_median_diff.png')); % Save the figure
    end

%% TODO Mean and Normalized Mean: Eye RT - Button RT - Diff [Log Y-Axis]


%% heatmaps
% figure;
% % Subplot 1: RTV for ButtonPress
% subplot(1, 3, 1);
% % Extract RTV values for ButtonPress
% rtv_matrix_button = cell2mat(arrayfun(@(p) ...
%     arrayfun(@(c) p.(conditions{c}).RT_ButtonPress_var, 1:length(conditions), 'UniformOutput', true), ...
%     rts, 'UniformOutput', false)');
% imagesc(rtv_matrix_button);
% colorbar;
% xticks(1:length(conditions));xticklabels(conditions);
% yticks(1:length(rts));yticklabels({rts.id});
% xlabel('Condition');ylabel('Participant');
% title('RTV Heatmap (ButtonPress)');
% 
% % Subplot 2: RTV for Eye
% subplot(1, 3, 2);
% % Extract RTV values for Eye
% rtv_matrix_eye = cell2mat(arrayfun(@(p) ...
%     arrayfun(@(c) p.(conditions{c}).RT_Eye_var, 1:length(conditions), 'UniformOutput', true), ...
%     rts, 'UniformOutput', false)');
% imagesc(rtv_matrix_eye);
% colorbar;
% xticks(1:length(conditions));xticklabels(conditions);
% yticks(1:length(rts));yticklabels({rts.id});
% xlabel('Condition');ylabel('Participant');
% title('RTV Heatmap (Eye)');
% 
% % Subplot 3: Time Difference (ButtonPress - Eye)
% subplot(1, 3, 3);
% % Calculate time differences
% time_diff_matrix = cell2mat(arrayfun(@(p) ...
%     arrayfun(@(c) p.(conditions{c}).RT_ButtonPress_var - p.(conditions{c}).RT_Eye_var, 1:length(conditions), 'UniformOutput', true), ...
%     rts, 'UniformOutput', false)');
% imagesc(time_diff_matrix);
% colorbar;
% xticks(1:length(conditions));xticklabels(conditions);
% yticks(1:length(rts));yticklabels({rts.id});
% xlabel('Condition');ylabel('Participant');
% title('Time Difference Heatmap (ButtonPress - Eye)');


end




%% Function to plot a single RT measure
function [y_min, y_max] = plotRTGroupComparison(rts, conditions, condition_labels, groups, color_map, rt_field, y_label)
    num_conditions = length(conditions);

    % Separate ADHD and non-ADHD indices
    adhd_indices = strcmp(groups, 'ADHD');
    nonadhd_indices = strcmp(groups, 'nonADHD');

    hold on;

    for c = 1:num_conditions
        condition = conditions{c};

        % Extract data
        adhd_values = arrayfun(@(p) p.(condition).(rt_field), rts(adhd_indices), 'UniformOutput', true);
        nonadhd_values = arrayfun(@(p) p.(condition).(rt_field), rts(nonadhd_indices), 'UniformOutput', true);

        % Compute means and SEMs
        adhd_mean = mean(adhd_values, 'omitnan');
        adhd_sem = std(adhd_values, 'omitnan') / sqrt(sum(~isnan(adhd_values)));

        nonadhd_mean = mean(nonadhd_values, 'omitnan');
        nonadhd_sem = std(nonadhd_values, 'omitnan') / sqrt(sum(~isnan(nonadhd_values)));

        % Plot means with error bars
        errorbar(c - 0.2, adhd_mean, adhd_sem, '.', 'Color', color_map('ADHD'), ...
            'MarkerFaceColor', color_map('ADHD'), 'LineWidth', 1.5);
        errorbar(c + 0.2, nonadhd_mean, nonadhd_sem, '.', 'Color', color_map('nonADHD'), ...
            'MarkerFaceColor', color_map('nonADHD'), 'LineWidth', 1.5);

        % Scatter individual points with jitter
        scatter(c - 0.2 + (rand(size(adhd_values)) - 0.5) * 0.15, adhd_values, 20, ...
            color_map('ADHD'), 'filled', 'MarkerEdgeColor', color_map('ADHD'), 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.3);
        scatter(c + 0.2 + (rand(size(nonadhd_values)) - 0.5) * 0.15, nonadhd_values, 20, ...
            color_map('nonADHD'), 'filled', 'MarkerEdgeColor', color_map('nonADHD'), 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.3);
    end

    % Customize plot
    xticks(1:num_conditions);
    xticklabels(condition_labels);
    ylabel(y_label);
    set(gca, 'YScale', 'log'); % Logarithmic scale for better visibility
    grid off;
    hold off;
    
    % Get y-axis limits
    y_limits = ylim();
    y_min = y_limits(1);
    y_max = y_limits(2);
end

%% Function to plot RT differences
function [y_min, y_max, adhd_handle, nonadhd_handle] = plotRTDifference(rts, conditions, condition_labels, groups, color_map, rt_field1, rt_field2, y_label)
    num_conditions = length(conditions);

    % Separate ADHD and non-ADHD indices
    adhd_indices = strcmp(groups, 'ADHD');
    nonadhd_indices = strcmp(groups, 'nonADHD');

    hold on;

    for c = 1:num_conditions
        condition = conditions{c};

        % Extract data
        adhd_values1 = arrayfun(@(p) p.(condition).(rt_field1), rts(adhd_indices), 'UniformOutput', true);
        adhd_values2 = arrayfun(@(p) p.(condition).(rt_field2), rts(adhd_indices), 'UniformOutput', true);
        nonadhd_values1 = arrayfun(@(p) p.(condition).(rt_field1), rts(nonadhd_indices), 'UniformOutput', true);
        nonadhd_values2 = arrayfun(@(p) p.(condition).(rt_field2), rts(nonadhd_indices), 'UniformOutput', true);

        % Compute differences
        adhd_diff = adhd_values1 - adhd_values2;
        nonadhd_diff = nonadhd_values1 - nonadhd_values2;

        
        %% --- for neg log values 
        % Apply a log transformation with sign preservation
        transform = @(x) sign(x) .* log10(1 + abs(x));

        % Transform your data
        adhd_diff = transform(adhd_diff);
        nonadhd_diff = transform(nonadhd_diff);

        %% ---
        
        
        
        % Compute means and SEMs
        adhd_mean_diff = mean(adhd_diff, 'omitnan');
        adhd_sem_diff = std(adhd_diff, 'omitnan') / sqrt(sum(~isnan(adhd_diff)));

        nonadhd_mean_diff = mean(nonadhd_diff, 'omitnan');
        nonadhd_sem_diff = std(nonadhd_diff, 'omitnan') / sqrt(sum(~isnan(nonadhd_diff)));

        % Plot differences with error bars
        if c == 1 % only return labels for legend once
            adhd_handle = errorbar(c - 0.2, adhd_mean_diff, adhd_sem_diff, '.', 'Color', color_map('ADHD'), ...
                    'MarkerFaceColor', color_map('ADHD'), 'LineWidth', 1.5);
            nonadhd_handle = errorbar(c + 0.2, nonadhd_mean_diff, nonadhd_sem_diff, '.', 'Color', color_map('nonADHD'), ...
                    'MarkerFaceColor', color_map('nonADHD'), 'LineWidth', 1.5);
        else
        	errorbar(c - 0.2, adhd_mean_diff, adhd_sem_diff, '.', 'Color', color_map('ADHD'), ...
                    'MarkerFaceColor', color_map('ADHD'), 'LineWidth', 1.5);
        	errorbar(c + 0.2, nonadhd_mean_diff, nonadhd_sem_diff, '.', 'Color', color_map('nonADHD'), ...
                    'MarkerFaceColor', color_map('nonADHD'), 'LineWidth', 1.5);
        end
        % Scatter individual differences with jitter
        scatter(c - 0.2 + (rand(size(adhd_diff)) - 0.5) * 0.15, adhd_diff, 20, ...
            color_map('ADHD'), 'filled', 'MarkerEdgeColor', color_map('ADHD'), 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.3);
        scatter(c + 0.2 + (rand(size(nonadhd_diff)) - 0.5) * 0.15, nonadhd_diff, 20, ...
            color_map('nonADHD'), 'filled', 'MarkerEdgeColor', color_map('nonADHD'), 'MarkerFaceAlpha', 0.3, 'MarkerEdgeAlpha', 0.3);
    end

    % Customize plot
    xticks(1:num_conditions);
    xticklabels(condition_labels);
    set(gca, 'YTick', transform([-100, -10, -1, 0, 1, 10, 100]));
    %set(gca, 'YTickLabel', {'-100', '-10', '-1', '0', '1', '10', '100'});
    % get nice latex format 10x power labels
    yticks = [-2, -1, 0, 1, 2]; % Logarithmic scale powers 
    set(gca, 'YTick', 10.^yticks, 'YTickLabel', arrayfun(@(x) ['$10^{' num2str(x) '}$'], yticks, 'UniformOutput', false));
    set(gca, 'TickLabelInterpreter', 'latex'); % Enable LaTeX interpreter for proper rendering
    
    % Match MATLAB default font size and style  
    ax = gca;% Access the current axes
    x_label_handle = ax.XLabel;% Get the XLabel properties
    set(gca, 'FontName', x_label_handle.FontName, 'FontSize', x_label_handle.FontSize);
    ylabel(y_label);
    %set(gca, 'YScale', 'log'); % Logarithmic scale for better visibility
    
    grid off; hold off;
    
    % Get y-axis limits
    y_limits = ylim();
    y_min = y_limits(1);
    y_max = y_limits(2);
end








function plotRTV_ParticipantLines(rt_variability, conditions, measure, color_map ,comparison_results_folder)
% This function plots the reaction time variability (RTV) for individual participants
% with lines colored according to their group.

% Create figure for ButtonPress variability
figure;

% Plot RT_ButtonPress variability for each participant
subplot(2, 1, 1);
hold on;

% Loop through participants and plot lines for RT_ButtonPress
for participant = 1:length(rt_variability)
    participant_group = rt_variability(participant).group;
    rt_button_var = [];
    
    for c = 1:length(conditions)
        condition = conditions{c};
        rt_button_var = [rt_button_var, rt_variability(participant).(condition).RT_ButtonPress_var];
    end
    
    % Plot the line for this participant, colored by group
    if strcmp(participant_group, 'ADHD')
        plot(1:length(conditions), rt_button_var, '-o', 'Color', color_map('ADHD'), 'DisplayName', 'ADHD', 'LineWidth', 1);
    else
        plot(1:length(conditions), rt_button_var, '-o', 'Color', color_map('nonADHD'), 'DisplayName', 'non-ADHD', 'LineWidth', 1);
    end
end

% Customize plot
set(gca, 'XTick', 1:length(conditions), 'XTickLabel', conditions);
xlabel('Condition');
ylabel(strcat('RT ButtonPress Variability: ', measure));
title('Reaction Time Variability (ButtonPress) per Participant');
grid on;
legend('Location', 'best');

% Create figure for Eye variability
subplot(2, 1, 2);
hold on;

% Plot RT_Eye variability for each participant
for participant = 1:length(rt_variability)
    participant_group = rt_variability(participant).group;
    rt_eye_var = [];
    
    for c = 1:length(conditions)
        condition = conditions{c};
        rt_eye_var = [rt_eye_var, rt_variability(participant).(condition).RT_Eye_var];
    end
    
    % Plot the line for this participant, colored by group
    if strcmp(participant_group, 'ADHD')
        plot(1:length(conditions), rt_eye_var, '-o', 'Color', color_map('ADHD'), 'DisplayName', 'ADHD', 'LineWidth', 1);
    else
        plot(1:length(conditions), rt_eye_var, '-o', 'Color', color_map('nonADHD'), 'DisplayName', 'non-ADHD', 'LineWidth', 1);
    end
end

% Customize plot
set(gca, 'XTick', 1:length(conditions), 'XTickLabel', conditions);
xlabel('Condition');
ylabel(strcat('RT Eye Variability: ', measure));
title('Reaction Time Variability (Eye) per Participant');

hold off;

grid on;

% Adjust the layout
%set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
% Save the figure if safe is set to 1
%     if safe == 1
%         set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf, fullfile(comparison_results_folder, 'RTV_participant_line.png'));
%     end
end
