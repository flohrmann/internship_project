function checkAccuracyDistributionByCondition(data, conditions)
    numConditions = length(conditions);

    % Loop through each condition
    for i = 1:numConditions
        condition = conditions{i};

        % Initialize counters for correct and incorrect trials
        num_correct = 0;
        num_incorrect = 0;

        % Loop through participants
        for participant = 1:length(data)
            % Loop through trials to check accuracy for the current condition
            for trial = 1:length(data(participant).Condition)
                if strcmp(data(participant).Condition{trial}, condition)
                    % Count the number of correct and incorrect trials
                    if data(participant).accuracy(trial) == 1
                        num_correct = num_correct + 1;
                    elseif data(participant).accuracy(trial) == 0
                        num_incorrect = num_incorrect + 1;
                    end
                end
            end
        end

        % Display the result for the current condition
        disp(['Condition: ', condition]);
        disp(['Number of correct trials (accuracy = 1): ', num2str(num_correct)]);
        disp(['Number of incorrect trials (accuracy = 0): ', num2str(num_incorrect)]);
        disp('------------------------------------------');
    end
end
