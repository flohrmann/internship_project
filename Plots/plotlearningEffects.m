function plotlearningEffects(data, screenXcm, screenYcm, screenXpixels, screenYpixels, color_map, comparison_results_folder, distance_threshold)
    % Define conditions and groups
    conditions = {'a', 'a_simple', 'b', 'b_simple'};
    groups = {'ADHD', 'nonADHD'};
    
    % Calculate screen center in pixels
    center_x = screenXpixels / 2;
    center_y = screenYpixels / 2;

    % Convert distance threshold from degrees to pixels
    distance_threshold_pixels = tand(distance_threshold) * (screenXpixels / (2 * atand(screenXcm / 2 / screenXpixels)));
    
    % Initialize variables to store reaction times by condition, group, and target location
    rt_central = struct('ADHD', [], 'nonADHD', []);
    rt_peripheral = struct('ADHD', [], 'nonADHD', []);

    % Loop through each participant and categorize trials as central or peripheral
    for i = 1:length(data)
        participant_group = data(i).group;  % Get the group of the current participant (ADHD or nonADHD)

        for j = 1:length(data(i).TargetPosition)
            % Extract the target position in pixels for each trial
            target_x = data(i).x_centers{j}(data(i).TargetPosition(j, 2), data(i).TargetPosition(j, 1));
            target_y = screenYpixels - data(i).y_centers{j}(data(i).TargetPosition(j, 2), data(i).TargetPosition(j, 1));

            % Calculate distance from the screen center to the target
            target_distance = sqrt((target_x - center_x)^2 + (target_y - center_y)^2);

            % Loop through each condition and categorize trials into central or peripheral
            for c = 1:length(conditions)
                condition = conditions{i}; 
                condition_rts = data(i).rt(strcmp(data(i).Condition, condition), :);
                condition_rt = condition_rts(j);  % Reaction time for the current condition

                if target_distance <= distance_threshold_pixels
                    % Central targets
                    rt_central.(participant_group).(conditions{c})(end+1) = condition_rt;
                else
                    % Peripheral targets
                    rt_peripheral.(participant_group).(conditions{c})(end+1) = condition_rt;
                end
            end
        end
    end

    % Calculate mean RTs for each condition, group, and target location
    central_means = structfun(@(g) structfun(@(c) mean(c, 'omitnan'), g, 'UniformOutput', false), rt_central, 'UniformOutput', false);
    peripheral_means = structfun(@(g) structfun(@(c) mean(c, 'omitnan'), g, 'UniformOutput', false), rt_peripheral, 'UniformOutput', false);

    % Set up the figure and subplots for central and peripheral targets
    figure;
    t = tiledlayout(1, 2, 'TileSpacing', 'compact');
    title(t, 'Reaction Times by Condition and Group (Central vs Peripheral Targets)', 'FontWeight', 'bold');

    % Subplot 1: Central Targets
    nexttile;
    hold on;
    for c = 1:length(conditions)
        % Plot the bars for ADHD and nonADHD within each condition
        bar(c - 0.2, central_means.ADH.(conditions{c}), 0.4, 'FaceColor', color_map('ADHD'));
        bar(c + 0.2, central_means.nonADHD.(conditions{c}), 0.4, 'FaceColor', color_map('nonADHD'));
    end
    set(gca, 'XTick', 1:length(conditions), 'XTickLabel', conditions);
    ylabel('Mean Reaction Time (RT)');
    xlabel('Condition');
    title('Central Targets');
    legend({'ADHD', 'nonADHD'}, 'Location', 'northeastoutside');
    hold off;

    % Subplot 2: Peripheral Targets
    nexttile;
    hold on;
    for c = 1:length(conditions)
        % Plot the bars for ADHD and nonADHD within each condition
        bar(c - 0.2, peripheral_means.ADH.(conditions{c}), 0.4, 'FaceColor', color_map('ADHD'));
        bar(c + 0.2, peripheral_means.nonADHD.(conditions{c}), 0.4, 'FaceColor', color_map('nonADHD'));
    end
    set(gca, 'XTick', 1:length(conditions), 'XTickLabel', conditions);
    ylabel('Mean Reaction Time (RT)');
    xlabel('Condition');
    title('Peripheral Targets');
    legend({'ADHD', 'nonADHD'}, 'Location', 'northeastoutside');
    hold off;

    % Save figure if specified
    if comparison_results_folder
        saveas(gcf, fullfile(comparison_results_folder, 'RT_by_Condition_Group_Central_Peripheral.png'));
    end
end
