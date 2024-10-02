function plotFixations(data, fixations, screenXpixels, screenYpixels, plot_these, analysis_folder, safe)
    % This function visualizes the gaze points and detected fixations for specified trials.
    % data: eye-tracking data structure
    % fixations: detected fixations structure
    % screenXpixels, screenYpixels: screen resolution in pixels
    % plot_these: vector of trial indices to plot
    % analysis_folder: folder where plots should be saved
    % safe: 1/0 indicating whether to save the plots
     
     
    for i = 1:length(plot_these)
        trialIdx = plot_these(i);

        % Find the fixations for this trial
        trialFixations = fixations.stimulusFixations([fixations.stimulusFixations.trial] == trialIdx);
        
        % get colourmap for fixation clusters
        numPoints = length(trialFixations.fixations);
        cm = colormap(jet(numPoints));  % 'jet' ranges from blue to red
        
        % Extract and convert gaze points to screen coordinates for the left eye
        %% change this from eye to stimulus trial or back if needed
        xGaze = (data.stimulusTrial(trialIdx).left.gazePoint.onDisplayArea(1, :)) * screenXpixels;
        yGaze = screenYpixels - (data.stimulusTrial(trialIdx).left.gazePoint.onDisplayArea(2, :)) * screenYpixels;
        %xGaze = (data.eyeTrial(trialIdx).left.gazePoint.onDisplayArea(1, :)) * screenXpixels;
        %yGaze = screenYpixels - (data.eyeTrial(trialIdx).left.gazePoint.onDisplayArea(2, :)) * screenYpixels;
        
        % Plot the gaze path
        figure;
        hold on;
        p1 = plot(xGaze, yGaze, '-.', 'LineWidth', 1.5, 'DisplayName', 'Gaze Path');
        
        % Overlay fixations
        for fixIdx = 1:length(trialFixations.fixations)
            trial = trialFixations.fixations(fixIdx);
            fixationPointsX = xGaze(trial.startIdx:trial.endIdx);
            fixationPointsY = yGaze(trial.startIdx:trial.endIdx);
            p2 = plot(fixationPointsX, fixationPointsY, '.', 'Color', cm(fixIdx,:), 'MarkerSize', 14, 'LineWidth', 2, 'DisplayName', '');%, 'DisplayName', 'Fixation');
            % Mark the center of the fixation
            p3 = plot(trial.center(1), trial.center(2), 'x', 'Color', cm(fixIdx,:), 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Fixation Center');
        end
        
        % Mark the center of the target
        targetCenters = fixations.targetCenters(:, trialIdx);
        %targetCenterY = fixations(trialIdx).targetCenters(2);
        p4 = plot(targetCenters(1), targetCenters(2), 'g*', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Target Center');
        
        % Customize the plot
        title(['Participant ', num2str(data.id), ' - Trial ', num2str(trialIdx)]);
        xlabel('X Coordinate (pixels)');
        ylabel('Y Coordinate (pixels)');
        %legend('show');
        legend([p1, p3, p4], {'Gaze Path', 'Fixation Centers', 'Target Center'}); 
        %grid on;
        
        % Set axis limits based on screen resolution
        axis([0 screenXpixels 0 screenYpixels]);  
        
        % Save the plot if "safe" is true
        if safe == 1
            set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
            saveas(gcf, fullfile(analysis_folder, ['Fixations_Participant_', num2str(data.id), '_Trial_', num2str(trialIdx), '.png']));
            %close(gcf); % Close the figure after saving
        end
    end
end
