function plotRTvsAccuracy(trial_results)
    rt = trial_results.rt;
    correct = trial_results.correct;
    
    figure;
    scatter(rt, correct);
    xlabel('Reaction Time (s)');
    ylabel('Accuracy (0 or 1)');
    title('RT vs. Accuracy');
end
