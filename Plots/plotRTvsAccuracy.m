function plotRTvsAccuracy(trial_results, analysis_folder)
    rt = trial_results.rt;
    correct = trial_results.correct;
    
    figure;
    scatter(rt, correct);
    xlabel('Reaction Time (s)');
    ylabel('Accuracy (0 or 1)');
    title('RT vs. Accuracy');
    saveas(gcf,strcat(analysis_folder, '\accuracy_vs_rt.png'));
end
