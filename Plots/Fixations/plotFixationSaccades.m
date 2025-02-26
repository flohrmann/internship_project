function plotFixationSaccades(data, all_saccades, all_fixations, screenXpixels, screenYpixels, plot_these, analysis_folder, safe)
% This function visualizes the gaze points, fixations, and saccades for specified trials.
% data: eye-tracking data structure
% all_saccades: struct containing detected saccades
% all_fixations: struct containing detected fixations
% screenXpixels, screenYpixels: screen resolution in pixels
% plot_these: vector of trial indices to plot
% analysis_folder: folder where plots should be saved
% safe: 1/0 indicating whether to save the plots
try
    participant_id = data.id;
catch
    participant_id = all_fixations.id;
end

try
    participant_fixations = all_fixations.fixations;
catch
    participant_fixations = all_fixations; %.stimulusFixations.fixations;
end


for i = 1:length(plot_these)
    trialIdx = plot_these(i);
    xGaze = (data.stimulusTrial(trialIdx).left.gazePoint.onDisplayArea(1, :)) * screenXpixels;
    yGaze = screenYpixels - (data.stimulusTrial(trialIdx).left.gazePoint.onDisplayArea(2, :)) * screenYpixels;
    
    % Plot the gaze path
    figure;
    p1 = plot(xGaze, yGaze, '-.', 'DisplayName', 'Gaze Path', 'Color', [0.5 0.5 0.5]);
    hold on;
    try
        % Overlay fixations
        fixations = participant_fixations.stimulusFixations(trialIdx).fixations;
        for fixIdx = 1:length(fixations)
            fixationPointsX = xGaze(fixations(fixIdx).startIdx: fixations(fixIdx).endIdx);
            fixationPointsY = yGaze(fixations(fixIdx).startIdx:fixations(fixIdx).endIdx);
            plot(fixationPointsX, fixationPointsY, 'ro-', 'LineWidth', 2, 'MarkerSize', 3);
            plot(fixations(fixIdx).center(1), fixations(fixIdx).center(2), 'kx', 'MarkerSize', 5, 'LineWidth', 2);
        end
        
        % Overlay saccades
        try
            saccades = all_saccades(trialIdx).saccades;
        catch
            participant_saccades = all_saccades.saccades;
            saccades = participant_saccades(trialIdx).saccades;
        end
        
        for saccIdx = 1:length(saccades)
            %saccadeStartIdx = saccades(saccIdx).startIdx;
            %saccadeEndIdx = saccades(saccIdx).endIdx;
            %plot([xGaze(saccadeStartIdx) xGaze(saccadeEndIdx)], ...
            %    [yGaze(saccadeStartIdx) yGaze(saccadeEndIdx)], ...
            %   'g-', 'LineWidth', 2, 'DisplayName', 'Saccade');
            p3 = plot([saccades(saccIdx).startCenter(1) saccades(saccIdx).endCenter(1)], ...
                [saccades(saccIdx).startCenter(2) saccades(saccIdx).endCenter(2)], ...
                'g-.', 'LineWidth', 2, 'DisplayName', 'Saccade');
        end
    catch
    end
    targetCenters = participant_fixations.targetCenters(:, trialIdx);
    p4 = plot(targetCenters(2), targetCenters(1), 'm*', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Target Center');
    
    title(['Participant ', num2str(participant_fixations.id), ': Trial ', num2str(trialIdx), ' - Fixation Length: ', participant_fixations.length]);
    xlabel('X Coordinate (pixels)');
    ylabel('Y Coordinate (pixels)');
    %legend('show');
    axis([0 screenXpixels 0 screenYpixels]);
    try
        legend([p1, p3, p4], {'Gaze Path', 'Fixation Centers', 'Target Center'});
    catch
        legend([p1, p4], {'Gaze Path', 'Target Center'});
        
    end
    if safe == 1
        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf, fullfile(analysis_folder, ['Saccades_Participant_', num2str(participant_id), '_Trial_', num2str(trialIdx), '_FixLength_', participant_fixations.length, '.png']));
        close(gcf);
    end
end
end
