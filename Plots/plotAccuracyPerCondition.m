function plotAccuracyPerCondition(trial_results)
    conditions = categorical(trial_results.Condition);
    correct = trial_results.correct;
    
    uniqueConditions = categories(conditions);
    numConditions = length(uniqueConditions);
    accuracy = zeros(numConditions, 1);
    
    for i = 1:numConditions
        conditionIndices = conditions == uniqueConditions{i};
        accuracy(i) = mean(correct(conditionIndices));
    end
    
    figure;
    bar(categorical(uniqueConditions), accuracy);
    xlabel('Condition');
    ylabel('Accuracy');
    title('Accuracy per Condition');
end
