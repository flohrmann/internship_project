function ttestResults = runTTests_ADHD_vs_nonADHD(...
    mean_rt_a_adhd, mean_rt_b_adhd, ...
    mean_rt_asimple_adhd, mean_rt_bsimple_adhd, ...
    mean_rt_a_nonadhd, mean_rt_b_nonadhd, ...
    mean_rt_asimple_nonadhd, mean_rt_bsimple_nonadhd, ...
    rtv_a_adhd, rtv_a_nonadhd,...
    rtv_b_adhd,rtv_b_nonadhd, ... 
    rtv_asimple_adhd,rtv_asimple_nonadhd,... 
    rtv_bsimple_adhd,rtv_bsimple_nonadhd, ...
    mean_accuracy_adhd, mean_accuracy_nonadhd, ...
    mean_con_adhd, mean_con_nonadhd, ...
    mean_6inatt_adhd, mean_inatt_nonadhd, ...
    count_6inatt_adhd, count_6inatt_nonadhd, ...
    mean_allinatt_adhd, mean_allinatt_nonadhd)

    % Initialize results structure
    ttestResults = struct();

    % Mean RT comparison
    ttestResults.mean_rt_a = runTTest(mean_rt_a_adhd, mean_rt_a_nonadhd, 'Mean RT (a)');
    ttestResults.mean_rt_b = runTTest(mean_rt_b_adhd, mean_rt_b_nonadhd, 'Mean RT (b)');
    ttestResults.mean_rt_asimple = runTTest(mean_rt_asimple_adhd, mean_rt_asimple_nonadhd, 'Mean RT (a simple)');
    ttestResults.mean_rt_bsimple = runTTest(mean_rt_bsimple_adhd, mean_rt_bsimple_nonadhd, 'Mean RT (b simple)');
    
    % RT variability (RTV) comparison
    ttestResults.rtv_a = runTTest(rtv_a_adhd, rtv_a_nonadhd, 'RTV (a)');
    ttestResults.rtv_b = runTTest(rtv_b_adhd, rtv_b_nonadhd, 'RTV (b)');
    ttestResults.rtv_asimple = runTTest(rtv_asimple_adhd, rtv_asimple_nonadhd, 'RTV (a simple)');
    ttestResults.rtv_bsimple = runTTest(rtv_bsimple_adhd, rtv_bsimple_nonadhd, 'RTV (b simple)');
    
    % Accuracy comparison
    ttestResults.accuracy = runTTest(mean_accuracy_adhd, mean_accuracy_nonadhd, 'Mean Accuracy');
    
    % Confusion comparison
    ttestResults.confusion = runTTest(mean_con_adhd, mean_con_nonadhd, 'Mean Confusion (a/b ratio)');
    
    % Inattention comparison
    ttestResults.mean_6inatt = runTTest(mean_6inatt_adhd, mean_inatt_nonadhd, '6-Item Inattention Mean');
    ttestResults.count_6inatt = runTTest(count_6inatt_adhd, count_6inatt_nonadhd, '6-Item Inattention Count');
    ttestResults.mean_allinatt = runTTest(mean_allinatt_adhd, mean_allinatt_nonadhd, 'All-Item Inattention Mean');
    
    % Display results
    disp('T-test Results:');
    disp(ttestResults);

end

function result = runTTest(group1, group2, variableName)
    % Perform two-sample t-test between two groups
    [h, p, ci, stats] = ttest2(group1, group2);
    
    % Store results in a struct
    result = struct();
    result.hypothesisTest = h;  % 0 means no significant difference, 1 means significant difference
    result.pValue = p;          % P-value less than 0.05 indicates significant difference
    result.confidenceInterval = ci;
    result.stats = stats;       % T-statistic, degrees of freedom, etc.
    
    % Display results for each variable
    fprintf('Variable: %s\n', variableName);
    fprintf('P-Value: %.4f\n', p);
    fprintf('T-Statistic: %.4f\n\n', stats.tstat);
end
