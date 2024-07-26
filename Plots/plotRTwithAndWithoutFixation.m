function plotRTwithAndWithoutFixation(trial_results, eye_rt, lookedAtFixation, analysis_folder)
    % Compute fixation summary
    fixationSummary = any(lookedAtFixation, 2);
    fixationSummary = double(fixationSummary); % Convert logical array to double for 1 and 0 output

    % Extract unique conditions from trial_results
    uniqueConditions = categorical(trial_results.Condition);
    uniqueConditions = categories(uniqueConditions);
    
    % Initialize figure for subplots
    figure;
    
    % Loop through each condition and create subplots
    for c = 1:length(uniqueConditions)
        condition = uniqueConditions{c};
        
        % Initialize arrays to store RTs and trial indices for the current condition
        RTs = [];
        trialIndices = [];
        colors = [];
        noFixationCount = 0;
        totalTrials = 0;

        % Loop through each trial in trial_results
        for trial = 1:size(trial_results, 1)
            % Extract the current trial data
            current_data = trial_results(trial, :);

            % Check if the current trial matches the condition
            if strcmp(current_data.Condition, condition)
                totalTrials = totalTrials + 1;

                % Extract the RT for the current trial
                rt_matlab = eye_rt.RTmatlab(trial);

                % Store the RT and trial index
                RTs = [RTs; rt_matlab];
                trialIndices = [trialIndices; trial];

                % Determine the color based on fixation status
                if fixationSummary(trial) == 1
                    colors = [colors; [0 1 0]]; % Green for fixation
                else
                    colors = [colors; [1 0 0]]; % Red for no fixation
                    noFixationCount = noFixationCount + 1;
                end
            end
        end

        % Calculate the percentage of trials with no fixation
        percentageNoFixation = (noFixationCount / totalTrials) * 100;

        % Determine the subplot position
        subplot(2, 2, c);
        
        % Plot RTs with colors indicating fixation status for the current condition
        hold on;
        scatter(trialIndices(fixationSummary(trialIndices) == 1), RTs(fixationSummary(trialIndices) == 1), 50, 'g', 'filled', 'DisplayName', 'Fixation');
        scatter(trialIndices(fixationSummary(trialIndices) == 0), RTs(fixationSummary(trialIndices) == 0), 50, 'r', 'filled', 'DisplayName', 'No Fixation');
        title(['RT per Trial - Condition: ', condition, ' (No Fixation: ', num2str(percentageNoFixation, '%.1f'), '%)']);
        xlabel('Trial');
        ylabel('Reaction Time (s)');
        grid on;
        
        % Add legend
        legend show;
        hold off;
    end

    % Save the combined figure
    saveas(gcf, fullfile(analysis_folder, 'RT_per_trial_with_fixation_status_by_condition.png'));
end
