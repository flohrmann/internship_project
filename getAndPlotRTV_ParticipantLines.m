function rt_variability = getAndPlotRTV_ParticipantLines(data, conditions, measure, groups, color_map, comparison_results_folder)
    % This function computes the variability (std dev or variance) of reaction times (RTs)
    % for each condition and group, and then plots the results using individual participant lines.
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
                
                % Store the RT variability for the respective group and condition
                if strcmp(group, 'ADHD')
                    rt_var_ADHD_button = [rt_var_ADHD_button; rt_button_var];
                    rt_var_ADHD_eye = [rt_var_ADHD_eye; rt_eye_var];
                elseif strcmp(group, 'non-ADHD')
                    rt_var_nonADHD_button = [rt_var_nonADHD_button; rt_button_var];
                    rt_var_nonADHD_eye = [rt_var_nonADHD_eye; rt_eye_var];
                end
                
                % Store the RT variability per participant for individual lines
                rt_variability(participant).id = data(participant).id;
                rt_variability(participant).group = group;
                rt_variability(participant).(condition).RT_ButtonPress_var = rt_button_var;
                rt_variability(participant).(condition).RT_Eye_var = rt_eye_var;
            end
        end
    end
    
    % Plot RT variability per participant and group with lines
    plotRTV_ParticipantLines(rt_variability, conditions,measure, color_map, comparison_results_folder);
    
    % Display the RT variability data
    disp('Reaction Time Variability (Standard Deviation or Variance) by Condition and Group for each participant:');
    disp(rt_variability);
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
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);    
        % Save the figure if safe is set to 1
%     if safe == 1
%         set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
     saveas(gcf, fullfile(comparison_results_folder, 'RTV_participant_line.png'));
%     end
end
