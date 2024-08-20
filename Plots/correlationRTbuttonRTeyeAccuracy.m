function [R, P] = correlationRTbuttonRTeyeAccuracy(data, comparison_results_folder)

    rt_button_press_all = [];
    rt_eye_all = [];
    accuracy_all = [];

    % Loop through data to aggregate all observations
    for i = 1:length(data)
        rt_button_press_all = [rt_button_press_all; data(i).rt];
        rt_eye_all = [rt_eye_all; data(i).rt_eye];
        accuracy_all = [accuracy_all; data(i).accuracy];
    end

    % Handle NaNs: remove rows where any data is NaN
    valid_idx = ~isnan(rt_button_press_all) & ~isnan(rt_eye_all) & ~isnan(accuracy_all);
    rt_button_press_all = rt_button_press_all(valid_idx);
    rt_eye_all = rt_eye_all(valid_idx);
    accuracy_all = accuracy_all(valid_idx);

    % Calculate correlations
    [R, P] = corr([rt_button_press_all, rt_eye_all, accuracy_all]);

    % Display results
    disp('Correlation matrix (R):');
    disp(R);

    disp('P-values matrix (P):');
    disp(P);

    % Save results to a .mat file
    results_struct.R = R;
    results_struct.P = P;
    save(strcat(comparison_results_folder, '\correlation_RTs_accuracy_results.mat'), 'results_struct');

    disp('Correlation and P-values matrices saved to correlation_results.mat');
end
