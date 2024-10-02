function [R_condition, P_condition] = correlationByCondition(data, conditions, comparison_results_folder)
    numConditions = length(conditions);

    % Initialize to store results for each condition
    R_condition = cell(numConditions, 1);
    P_condition = cell(numConditions, 1);

    % Loop through each condition
    for i = 1:numConditions
        condition = conditions{i};

        % Initialize empty arrays to store the data for this condition
        rt_button = [];
        rt_eye = [];
        accuracy = [];

        % Loop through participants
        for participant = 1:length(data)
            % Loop through trials to filter data by condition
            for trial = 1:length(data(participant).Condition)
                if strcmp(data(participant).Condition{trial}, condition)
                    % Append the data for the current trial
                    rt_button = [rt_button; data(participant).rt(trial)];
                    rt_eye = [rt_eye; data(participant).rt_eye(trial)];
                    accuracy = [accuracy; data(participant).accuracy(trial)];
                end
            end
        end

        % Handle NaNs: remove rows where any data is NaN
        valid_idx = ~isnan(rt_button) & ~isnan(rt_eye) & ~isnan(accuracy);
        rt_button = rt_button(valid_idx);
        rt_eye = rt_eye(valid_idx);
        accuracy = accuracy(valid_idx);

        % Perform correlation analysis
        [R, P] = corr([rt_button, rt_eye, accuracy]);

        % Store results for this condition
        R_condition{i} = R;
        P_condition{i} = P;

        % Display results for the condition
        disp(['Condition: ', condition]);
        disp('Correlation Matrix (R):');
        disp(R);
        disp('P-values Matrix (P):');
        disp(P);
    end

    % Save the results
    results_struct.R_condition = R_condition;
    results_struct.P_condition = P_condition;
    save(fullfile(comparison_results_folder, 'correlation_by_condition_results.mat'), 'results_struct');

    disp('Condition-specific correlation results saved.');
end
