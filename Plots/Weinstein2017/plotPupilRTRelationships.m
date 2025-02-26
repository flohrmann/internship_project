function plotPupilRTRelationships(avg_adhd, avg_nonadhd, conditions, color_map, folder)
condition_names = {'a', 'a simple', 'b', 'b simple'};
pupil_fields = {'a_max', 'as_max', 'b_max', 'bs_max'};
rtv_fields   = {'rtv_eye_a', 'rtv_eye_as', 'rtv_eye_b', 'rtv_eye_bs'};


%% --- per group ---
figure;
t = tiledlayout('flow' ,'TileSpacing', 'compact');
% Loop through each condition and create the subplots
for i = 1:length(conditions)
    % ADHD
    max_pupil_adhd = avg_adhd.(pupil_fields{i});
    rtv_adhd = avg_adhd.(rtv_fields{i});
    % non-ADHD
    max_pupil_nonadhd = avg_nonadhd.(pupil_fields{i});
    rtv_nonadhd = avg_nonadhd.(rtv_fields{i});
    
    nexttile;
    hold on;
    % Plot ADHD group for the current condition
    sz = 30;
    scatter(rtv_adhd, max_pupil_adhd, sz, ...
            'MarkerEdgeColor', color_map('ADHD'), ...
            'MarkerFaceColor', color_map('ADHD'), ...
            'DisplayName', 'ADHD', ...
            'Marker', 'o', ...
            'LineWidth', 2);
        
    % Plot non-ADHD group for the current condition
    sz = 70;
    scatter(rtv_nonadhd, max_pupil_nonadhd, sz,...
        'MarkerEdgeColor', color_map('nonADHD'), ...
        'MarkerFaceColor', color_map('nonADHD'), ...
        'DisplayName', 'Non-ADHD', ...
        'Marker', 'x', ...
        'LineWidth', 3);
    
    xlabel('Reaction Time Standard Deviation (RTV)');
    ylabel('Max Pupil Diameter (z-score)');
    title(condition_names(i));
    grid off;
    hold off;
end
title(t, 'Max Pupil Dilation vs. Reaction Time Variability (RTV) - ADHD vs. Non-ADHD');
lgd = legend({'ADHD', 'Non-ADHD'}, 'Location', 'northwest', 'Orientation', 'vertical');
lgd.Location = 'bestoutside';
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1])
saveas(gcf, fullfile(folder, 'RTV_Pupil_Relationships_ADHD_vs_nonADHD.png'));


%% --- per participant ---
    figure;
    t = tiledlayout('flow', 'TileSpacing', 'compact');

    % Define color gradients using summer (for ADHD) and autumn (for non-ADHD) colormaps
    num_adhd = height(avg_adhd);
    num_nonadhd = height(avg_nonadhd);

    adhd_colors = colormap(summer(num_adhd+3));      % Summer colormap for ADHD
    nonadhd_colors = colormap(autumn(num_nonadhd+3)); % Autumn colormap for non-ADHD
    
    % Loop through each condition and create the subplots
    for i = 1:length(conditions)
        % ADHD group data for the current condition
        max_pupil_adhd = avg_adhd.(pupil_fields{i});
        rtv_adhd = avg_adhd.(rtv_fields{i});

        % Non-ADHD group data for the current condition
        max_pupil_nonadhd = avg_nonadhd.(pupil_fields{i});
        rtv_nonadhd = avg_nonadhd.(rtv_fields{i});

        % Create the tile
        nexttile;
        hold on;

        % Plot ADHD group using the summer colormap
        for j = 1:length(rtv_adhd)
            current_color = adhd_colors(j, :);  % Use summer colors
            sz = 30;
            scatter(rtv_adhd(j), max_pupil_adhd(j), sz, ...
                'MarkerEdgeColor', current_color, ...
                'MarkerFaceColor', current_color, ...
                'DisplayName', ['ADHD ', num2str(avg_adhd.id(j))], ...
                'Marker', 'o', ...
                'LineWidth', 3);
        end

        % Plot non-ADHD group using the autumn colormap
        for k = 1:length(rtv_nonadhd)
            current_color = nonadhd_colors(k, :);  % Use autumn colors
            sz = 70;
            scatter(rtv_nonadhd(k), max_pupil_nonadhd(k), sz, ...
                'MarkerEdgeColor', current_color, ...
                'MarkerFaceColor', current_color, ...
                'DisplayName', ['nonADHD ', num2str(avg_nonadhd.id(k))], ...
                'Marker', 'x', ...
                'LineWidth', 3);
        end

        % Label the axes and set title for each subplot
        xlabel('Reaction Time Standard Deviation (RTV)');
        ylabel('Max Pupil Diameter (z-score)');
        title(condition_names{i});
        grid off;
        hold off;
    end

    % Global title for the entire figure
    title(t, 'Max Pupil Dilation vs. Reaction Time Variability (RTV) - ADHD vs. Non-ADHD');

    % Create legend in the last tile with participant IDs
    lgd = legend('Location', 'bestoutside');
    lgd.Layout.Tile = 'east';  % Position legend to the right side

    % Adjust and save the figure
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf, fullfile(folder, 'RTV_Pupil_Relationships_ADHD_vs_nonADHD_per_participant.png'));

end