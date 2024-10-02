function rt_variability = getRTvariability(data, conditions, groups, measure, color_map, comparison_results_folder)
        % This function computes the variability (std dev or variance) of reaction times (RTs)
    % for each condition and group, and then plots the results using line plots with error bars.
    % measure: 'std' for standard deviation, 'var' for variance.
    
    % Initialize a struct to store the results
    rt_variability = struct();
    
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
                % Get the group for this participant (ADHD or nonADHD)
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
                
                % Store the RT variability for the respective group and condition
                if strcmp(group, 'ADHD')
                    rt_var_ADHD_button = [rt_var_ADHD_button; rt_button_var];
                    rt_var_ADHD_eye = [rt_var_ADHD_eye; rt_eye_var];
                elseif strcmp(group, 'nonADHD')
                    rt_var_nonADHD_button = [rt_var_nonADHD_button; rt_button_var];
                    rt_var_nonADHD_eye = [rt_var_nonADHD_eye; rt_eye_var];
                end
            end
        end
        
        % Store the variability data in the struct
        rt_variability.(condition).ADHD.RT_ButtonPress_var = rt_var_ADHD_button;
        rt_variability.(condition).ADHD.RT_Eye_var = rt_var_ADHD_eye;
        rt_variability.(condition).nonADHD.RT_ButtonPress_var = rt_var_nonADHD_button;
        rt_variability.(condition).nonADHD.RT_Eye_var = rt_var_nonADHD_eye;
    end
    
    % Plot RT variability per condition and group with error bars
    plotRTV_LinePlot(rt_variability, conditions, measure, color_map, comparison_results_folder);
    
    % Display the RT variability data
    disp('Reaction Time Variability (Standard Deviation or Variance) by Condition and Group:');
    disp(rt_variability);
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