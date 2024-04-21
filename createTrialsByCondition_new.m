function condition_trials = createTrialsByCondition_new(NumberOfBlocks, NTrialsEachCondition, trial_data, conditions)
% Gets trial matrix and conditions
% Returns trial matrix filled with angles of conditions
% Bsp. 2 3 2 3 2 3  + a_simple_1 = 45 45 45  45 45 45
%      2 3 1 3 2 3                 45 45 135 45 45 45

    % Initialize an empty table to store data for each condition
    condition_trials = table;

    % Get the names of the conditions
    condition_names = fieldnames(conditions);

    % Loop through each condition
    for i = 1:numel(condition_names)
        condition_name = condition_names{i};
        trials_per_condition = NTrialsEachCondition(i) * NumberOfBlocks;
        % Extract trial data for the current condition
        start_index = (i-1) * trials_per_condition + 1;
        end_index = start_index + trials_per_condition - 1;
        condition_trials_data = trial_data(start_index:end_index, :);

        % Determine sub-conditions or partitions if any
        partitions = size(conditions.(condition_name), 1);

        % Assign angles to the TrialMatrix based on the current condition
        for j = 1:trials_per_condition
            partition_index = ceil(j / (trials_per_condition / partitions));
            angles = conditions.(condition_name)(partition_index, :);

             % Assign each trial
            trial = condition_trials_data(j, :);
            if iscell(angles)  % Checking if the angles are in a cell array (for complex conditions with tuples)
                for idx = 1:numel(trial.TrialMatrix{1})
                    tuple_idx = trial.TrialMatrix{1}(idx);
                    angle_set = angles{tuple_idx};  % Access the correct tuple for target/distractors
                    % Now, handle both angles
                    % Assuming each angle set in the tuple has exactly two angles
                    trial.TrialMatrix{1}{idx} = angle_set;
                end
            else
                for idx = 1:numel(trial.TrialMatrix{1})
                    trial.TrialMatrix{1}(idx) = angles(trial.TrialMatrix{1}(idx));
                end
            end

        end
            % Set additional trial information
            trial.Condition = {condition_name};
            trial.ConditionType = {strcat(condition_name, '_', num2str(partition_index))};
            trial.TargetAngle = angles(1);  % Target 
            trial.DistractorAngle_1 = angles(2);  % Distractor 1 
            trial.DistractorAngle_2 = angles(3);  % Distractor 1 

            % Append to the main table
            condition_trials = [condition_trials; trial];
   end
end

% 
% function condition_trials = createTrialsByCondition(NumberOfBlocks, NTrialsEachCondition, trial_data, conditions)
%     % Initialize table to store data for each condition
%     condition_trials = [];
% 
%     % Get the names of the conditions
%     condition_names = fieldnames(conditions);
%     % Number of conditions
%     num_conditions = numel(condition_names);
% 
%     % Loop through each condition
%     for i = 1:num_conditions
%         condition_name = condition_names{i};
%         trials_per_condition = NTrialsEachCondition(i)*NumberOfBlocks;
%         % Extract trial data for the current condition
%         start_index = (i-1) * trials_per_condition + 1;
%         end_index = start_index + trials_per_condition -1;
%         condition_trials_data = trial_data(start_index:end_index, :);
% 
%         % Fill in angles for the bars based on condition
%         if strcmp(condition_name, 'a_simple')
%             for j = 1:trials_per_condition
%                 if j <= trials_per_condition/2
%                     % Add trial data for the current condition
%                     trial = condition_trials_data(j, :);
%                     trial.TrialMatrix{1}(trial.TrialMatrix{1} == 1) = conditions.a_simple(1); % stim
%                     trial.TrialMatrix{1}(trial.TrialMatrix{1} == 2) = conditions.a_simple(2); % distr 1
%                     trial.TrialMatrix{1}(trial.TrialMatrix{1} == 3) = conditions.a_simple(3); % distr 2
% 
%                     trial.Condition = {condition_name};
%                     trial.ConditionType = {strcat(condition_name,'_1')};
%                     trial.TargetAngle = conditions.a_simple(1);
%                     trial.NormalAngle = conditions.a_simple(2);
%                     condition_trials = [condition_trials; trial];
%                 else
%                     trial = condition_trials_data(j, :);
%                     trial.TrialMatrix{1}(trial.TrialMatrix{1} == 1) = conditions.a_simple(4);
%                     trial.TrialMatrix{1}(trial.TrialMatrix{1} == 2) = conditions.a_simple(5);
%                     trial.TrialMatrix{1}(trial.TrialMatrix{1} == 3) = conditions.a_simple(6);
%                     trial.Condition = {condition_name};
%                     trial.ConditionType = {strcat(condition_name,'_2')};
%                     trial.TargetAngle = conditions.a_simple(4);
%                     trial.NormalAngle = conditions.a_simple(5);
%                     condition_trials = [condition_trials; trial];
%                 end
%             end
%         elseif strcmp(condition_name, 'a')
%             for j = 1:trials_per_condition
%                 if j <= trials_per_condition/4
%                     % Add trial data for the current condition
%                     trial = condition_trials_data(j, :);
%                     trial.TrialMatrix{1}(trial.TrialMatrix{1} == 1) = conditions.a_simple(1);
%                     trial.TrialMatrix{1}(trial.TrialMatrix{1} == 2) = conditions.a_simple(2);
%                     trial.TrialMatrix{1}(trial.TrialMatrix{1} == 3) = conditions.a_simple(3);
% 
%                     trial.Condition = {condition_name};
%                     trial.ConditionType = {strcat(condition_name,'_1')};
%                     trial.TargetAngle = conditions.a_simple(1);
%                     trial.NormalAngle = conditions.a_simple(2);
%                     condition_trials = [condition_trials; trial];
%                 else
%                     trial = condition_trials_data(j, :);
%                     trial.TrialMatrix{1}(trial.TrialMatrix{1} == 1) = conditions.a_simple(4);
%                     trial.TrialMatrix{1}(trial.TrialMatrix{1} == 2) = conditions.a_simple(5);
%                     trial.TrialMatrix{1}(trial.TrialMatrix{1} == 3) = conditions.a_simple(6);
%                     trial.Condition = {condition_name};
%                     trial.ConditionType = {strcat(condition_name,'_2')};
%                     trial.TargetAngle = conditions.a_simple(4);
%                     trial.NormalAngle = conditions.a_simple(5);
%                     condition_trials = [condition_trials; trial];
%                 end
%             end
% 
%         elseif strcmp(condition_name, 'b_simple')
%             for j = 1:trials_per_condition
%                 if j <= ceil(trials_per_condition*0.25)%8
%                     target_angle = conditions.b_simple(1);
%                     angle_2 = conditions.b_simple(2);
%                     angle_3 = conditions.b_simple(3);
%                     trial.ConditionType = {strcat(condition_name,'_1')};
%                 elseif j <= ceil(trials_per_condition*0.5)%15
%                     target_angle = conditions.b_simple(4);
%                     angle_2 = conditions.b_simple(5);
%                     angle_3 = conditions.b_simple(6);
%                     trial.ConditionType = {strcat(condition_name,'_2')};
%                 elseif j <= ceil(trials_per_condition * 0.75)%23
%                     target_angle = conditions.b_simple(7);
%                     angle_2 = conditions.b_simple(8);
%                     angle_3 = conditions.b_simple(9);
%                     trial.ConditionType = {strcat(condition_name,'_3')};
%                 else
%                     target_angle = conditions.b_simple(10);
%                     angle_2 = conditions.b_simple(11);
%                     angle_3 = conditions.b_simple(12);
%                     trial.ConditionType = {strcat(condition_name,'_4')};
%                 end
%                 % Add trial data for the current condition
%                 trial = condition_trials_data(j, :);
%                 trial.TrialMatrix{1}(trial.TrialMatrix{1} == 1) = target_angle;
%                 trial.TrialMatrix{1}(trial.TrialMatrix{1} == 2) = angle_2;
%                 trial.TrialMatrix{1}(trial.TrialMatrix{1} == 3) = angle_3;
% 
%                 trial.Condition = {condition_name};
%                 trial.TargetAngle = target_angle;
%                 trial.NormalAngle = angle_2;
%                 condition_trials = [condition_trials; trial];
%             end
%         end
%         % Store condition data
%         %conditionData.(genvarname(condition_name)) = condition_trials;
%     end
% end