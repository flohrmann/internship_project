function plotRTDistributionPerCondition(trial_results, analysis_folder)
    conditions = categorical(trial_results.Condition);
    rt = trial_results.rt;
    uniqueConditions = categories(conditions);
    numConditions = length(uniqueConditions);
    
    figure;
    hold on;
    for i = 1:numConditions
        conditionIndices = conditions == uniqueConditions{i};
        rt_condition = rt(conditionIndices);
        histogram(rt_condition, 'DisplayName', char(uniqueConditions{i}), 'Normalization', 'pdf');
    end
    hold off;
    
    xlabel('Reaction Time (s)');
    ylabel('Density');
    title('Reaction Time Distribution per Condition');
    legend('show');
    saveas(gcf,strcat(analysis_folder, '\rt_distr_condition.png'));

end
