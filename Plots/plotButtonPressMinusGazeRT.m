function plotButtonPressMinusGazeRT(trial_results, eye_rt, analysis_folder)
    % Initialize arrays to store the differences and colors
    RT_differences = [];
    colors = [];

    % Loop through each trial in trial_results
    for trial = 1:size(trial_results, 1)
        % Extract the RT from the button press
        rt_button_press = eye_rt.RTmatlab(trial);

        % Extract the gaze RTs for both eyes
        rightEyeRT = eye_rt.RightEyeRT(trial);
        leftEyeRT = eye_rt.LeftEyeRT(trial);

        % Calculate the mean gaze RT (ignoring zeros)
        if rightEyeRT == 0
            gazeRT = leftEyeRT;
        elseif leftEyeRT == 0
            gazeRT = rightEyeRT;
        else
            gazeRT = mean([rightEyeRT, leftEyeRT]);
        end

        % Calculate the difference: Button press RT - Gaze RT
        RT_difference = rt_button_press - gazeRT;

        % Store the difference
        RT_differences = [RT_differences; RT_difference];

        % Determine the color based on whether the press was correct or incorrect
        if trial_results.correct(trial) % Assuming 'Correct' is a logical field in trial_results
            colors = [colors; [0 1 0]]; % Green for correct
        else
            colors = [colors; [1 0 0]]; % Red for incorrect
        end
    end

    % Plot the RT differences
    figure;
    scatter(1:length(RT_differences), RT_differences, 50, colors, 'filled', 'DisplayName', 'Correct Button Pressed');
    title('Button Press RT - Gaze RT per Trial');
    xlabel('Trial');
    ylabel('RT Difference (Button Press - Gaze RT) (s)');
    
    % Create dummy scatter plots for the legend with correct colors
    hold on;
    scatter(NaN, NaN, 50, [1 0 0], 'filled', 'DisplayName', 'Wrong Button Pressed');
    %scatter(NaN, NaN, 50, [0 1 0], 'filled', 'DisplayName', 'Correct Button Pressed');
    hold off;
    
    % Show legend
    legend('show');
    grid on;

    % Save the plot
    saveas(gcf, fullfile(analysis_folder, 'ButtonPress_minus_GazeRT_with_correctness.png'));
end
