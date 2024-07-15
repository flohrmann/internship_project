function plotAccuracyPerCondition(trial_results, analysis_folder)
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
    saveas(gcf,strcat(analysis_folder, '\accuracy_per_condition.png'));
end
