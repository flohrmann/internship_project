function condition_trials = createTrialsByCondition_new(n_trials, trial_data, conditions)
    % Initialize table to store data for each condition
    condition_trials = [];

    % Determine the number of selected conditions
    % round in case of not divideable
    conditions_split = strsplit(conditions, ',');
    num_selected_conditions = length(conditions_split);
    trials_per_condition = ceil(n_trials / num_selected_conditions);

    % Loop through each selected condition
    for i = 1:num_selected_conditions
        condition_name = conditions_split{i};
        
        % Extract trial data for the current condition
        start_index = (i-1) * trials_per_condition + 1;
        end_index = start_index + trials_per_condition -1;
        condition_trials_data = trial_data(start_index:end_index, :);

        % Determine angles for the target bars based on condition
        if contains(condition_name, {'a_simple', 'a'})
            for j = 1:trials_per_condition
                % Select first half of trials to have a target bar with 135 degrees
                % and the other half to have a target bar with 45 degrees
                if j <= trials_per_condition/2
                    % Add trial data for the current condition
                    trial = condition_trials_data(j, :);
                    trial.TrialMatrix{1}(trial.TrialMatrix{1} == 0) = 45; % All bars are 45 degrees from vertical
                    trial.TrialMatrix{1}(trial.TrialMatrix{1} == 1) = 135; % Stimulus angle is opposite
                    trial.Condition = {condition_name};
                    trial.TargetAngle = 135;
                    trial.NormalAngle = 45;
                    condition_trials = [condition_trials; trial];
                else
                    trial = condition_trials_data(j, :);
                    trial.TrialMatrix{1}(trial.TrialMatrix{1} == 0) = 135; % All bars are 45 degrees from vertical
                    trial.TrialMatrix{1}(trial.TrialMatrix{1} == 1) = 45; % Stimulus angle is opposite
                    trial.Condition = {condition_name};
                    trial.TargetAngle = 45;
                    trial.NormalAngle = 135;
                    condition_trials = [condition_trials; trial];
                end
            end
        else
            % For B and B_simple, generate 1 of four possible 20-degree rotations
            %target_angles = [-20, 20, 70, 110];
            for j = 1:trials_per_condition
                % Determine the target angle based on the condition index
                if j <= ceil(trials_per_condition*0.25)%8
                    angle = 45;
                    target_angle = angle + 20;
                elseif j <= ceil(trials_per_condition*0.5)%15
                    angle = 45;
                    target_angle = angle - 20;
                elseif j <= ceil(trials_per_condition * 0.75)%23
                    angle = 135;
                    target_angle = angle + 20;
                else
                    angle = 135;
                    target_angle = angle - 20;
                end
                % Add trial data for the current condition
                trial = condition_trials_data(j, :);
                trial.TrialMatrix{1}(trial.TrialMatrix{1} == 1) = target_angle;
                trial.TrialMatrix{1}(trial.TrialMatrix{1} == 0) = angle;
                trial.Condition = {condition_name};
                trial.TargetAngle = target_angle;
                trial.NormalAngle = angle;
                condition_trials = [condition_trials; trial];
            end
        end
        % Store condition data
        %conditionData.(genvarname(condition_name)) = condition_trials;
    end
end