function [R_group, P_group] = correlationByGroup(data, groups, comparison_results_folder)
    numGroups = length(groups);

    % Initialize to store results for each group
    R_group = cell(numGroups, 1);
    P_group = cell(numGroups, 1);

    % Loop through each group
    for i = 1:numGroups
        group = groups{i};

        % Initialize empty arrays to store the data for this group
        rt_button = [];
        rt_eye = [];
        accuracy = [];

        % Loop through participants
        for participant = 1:length(data)
            if strcmp(data(participant).group, group)
                % Loop through trials to gather data for this participant
                rt_button = [rt_button; data(participant).rt];
                rt_eye = [rt_eye; data(participant).rt_eye];
                accuracy = [accuracy; data(participant).accuracy];
            end
        end

        % Handle NaNs: remove rows where any data is NaN
        valid_idx = ~isnan(rt_button) & ~isnan(rt_eye) & ~isnan(accuracy);
        rt_button = rt_button(valid_idx);
        rt_eye = rt_eye(valid_idx);
        accuracy = accuracy(valid_idx);

        % Perform correlation analysis
        [R, P] = corr([rt_button, rt_eye, accuracy]);

        % Store results for this group
        R_group{i} = R;
        P_group{i} = P;

        % Display results for the group
        disp(['Group: ', group]);
        disp('Correlation Matrix (R):');
        disp(R);
        disp('P-values Matrix (P):');
        disp(P);
    end

    % Save the results
    results_struct.R_group = R_group;
    results_struct.P_group = P_group;
    save(fullfile(comparison_results_folder, 'correlation_by_group_results.mat'), 'results_struct');

    disp('Group-specific correlation results saved.');
end
