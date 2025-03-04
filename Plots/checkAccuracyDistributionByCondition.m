function checkAccuracyDistributionByCondition(data, conditions)
    numConditions = length(conditions);

    % Loop through each condition
    for i = 1:numConditions
        condition = conditions{i};

        % Initialize counters for correct and incorrect trials
        total_correct = 0;
        total_incorrect = 0;

        % Loop through participants
        for p = 1:length(data)
            % Loop through trials to check accuracy for the current condition
            num_correct = 0;
            num_incorrect = 0;
            for trial = 1:length(data(p).Condition)
                if strcmp(data(p).Condition{trial}, condition)
                    % Count the number of correct and incorrect trials
                    if data(p).accuracy(trial) == 1
                        num_correct = num_correct + 1;
                    elseif data(p).accuracy(trial) == 0
                        num_incorrect = num_incorrect + 1;
                    end
                end
            end
            disp(['ID:',num2str(data(p).id), ': ', num2str(num_incorrect), '/',num2str(length(data(p).Condition)/4) ,' incorrect trials']);
            total_correct   = total_correct + num_correct; 
            total_incorrect = total_incorrect + num_incorrect;
        end

        % Display the result for the current condition
        disp(['Condition: ', condition]);
        disp(['Number of correct trials (accuracy = 1): ', num2str(total_correct)]);
        disp(['Number of incorrect trials (accuracy = 0): ', num2str(total_incorrect)]);
        disp('------------------------------------------');
    end
end
