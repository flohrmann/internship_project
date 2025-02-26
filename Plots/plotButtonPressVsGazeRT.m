function plotButtonPressVsGazeRT(trial_results, eye_rt, analysis_folder)
    % Initialize arrays to store RTmatlab and mean eye RTs
    RTmatlab = [];
    meanEyeRTs = [];

    % Loop through each trial in trial_results
    for trial = 1:size(trial_results, 1)
        % Extract the RTmatlab for the current trial
        rt_matlab = eye_rt.RTmatlab(trial);

        % Extract the eye RTs for the current trial
        rightEyeRT = eye_rt.RightEyeRT(trial);
        leftEyeRT = eye_rt.LeftEyeRT(trial);

        % Calculate the mean of the eye RTs, ignoring zeros
        if rightEyeRT == 0
            meanEyeRT = leftEyeRT;
        elseif leftEyeRT == 0
            meanEyeRT = rightEyeRT;
        else
            meanEyeRT = nanmean([rightEyeRT, leftEyeRT]);
        end

        % Store the values in the arrays
        RTmatlab = [RTmatlab; rt_matlab];
        meanEyeRTs = [meanEyeRTs; meanEyeRT];
    end

    % Plot RTmatlab vs. mean eye RTs
    figure;
    scatter(RTmatlab, meanEyeRTs, 'filled');
    title('Button Press RT vs. Gaze RT');
    xlabel('RT from Button Press (s)');
    ylabel('Mean Gaze RT (s)');
    grid on;

    % Add a regression line to the plot
    hold on;
    p = polyfit(RTmatlab, meanEyeRTs, 1); % Linear fit
    yfit = polyval(p, RTmatlab);
    plot(RTmatlab, yfit, '-k', 'LineWidth', 1.5); % Regression line in red
    hold off;

    % Save the plot
    saveas(gcf, fullfile(analysis_folder, 'button_press_vs_gaze_rt.png'));
end
