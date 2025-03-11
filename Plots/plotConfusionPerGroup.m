function data = plotConfusionPerGroup(average, data, group_labels, ids, conditions, color_map,color_map_other, participant_map, comparison_results_folder, safe)
disp([average, ' confusion'])

    % Initialize arrays to hold differences and ratios for eye-based and button-press metrics
    diff_adhd_eye = []; diff_nonadhd_eye = [];
    ratio_adhd_eye = []; ratio_nonadhd_eye = [];
    diff_adhd_press = []; diff_nonadhd_press = [];
    ratio_adhd_press = []; ratio_nonadhd_press = [];

    % Loop through each participant and calculate differences and ratios
    for i = 1:length(data)
        if strcmp(data(i).group, 'ADHD')
            % Eye-based metrics
            diff_adhd_eye = [diff_adhd_eye; data(i).RTa_eye - data(i).RTb_eye];
            ratio_adhd_eye = [ratio_adhd_eye; data(i).RTa_eye / data(i).RTb_eye];
            data(i).confusion_diff_eye = data(i).RTa_eye - data(i).RTb_eye;
            data(i).confusion_ratio_eye = data(i).RTa_eye / data(i).RTb_eye;
            
            % Button-press metrics
            diff_adhd_press = [diff_adhd_press; data(i).RTa - data(i).RTb];
            ratio_adhd_press = [ratio_adhd_press; data(i).RTa / data(i).RTb];
            data(i).confusion_diff_press = data(i).RTa - data(i).RTb;
            data(i).confusion_ratio_press = data(i).RTa / data(i).RTb;
            
        elseif strcmp(data(i).group, 'nonADHD')
            % Eye-based metrics
            diff_nonadhd_eye = [diff_nonadhd_eye; data(i).RTa_eye - data(i).RTb_eye];
            ratio_nonadhd_eye = [ratio_nonadhd_eye; data(i).RTa_eye / data(i).RTb_eye];
            data(i).confusion_diff_eye = data(i).RTa_eye - data(i).RTb_eye;
            data(i).confusion_ratio_eye = data(i).RTa_eye / data(i).RTb_eye;
            
            % Button-press metrics
            diff_nonadhd_press = [diff_nonadhd_press; data(i).RTa - data(i).RTb];
            ratio_nonadhd_press = [ratio_nonadhd_press; data(i).RTa / data(i).RTb];
            data(i).confusion_diff_press = data(i).RTa - data(i).RTb;
            data(i).confusion_ratio_press = data(i).RTa / data(i).RTb;
        end
    end

    % Calculate means/medians and SEMs 
    if strcmp(average, 'Mean')
        avg_diff_adhd_eye    = mean(diff_adhd_eye, 'omitnan')
        avg_diff_nonadhd_eye = mean(diff_nonadhd_eye, 'omitnan')

        avg_ratio_adhd_eye    = mean(ratio_adhd_eye, 'omitnan')
        avg_ratio_nonadhd_eye = mean(ratio_nonadhd_eye, 'omitnan')

        avg_diff_adhd_press    = mean(diff_adhd_press, 'omitnan')
        avg_diff_nonadhd_press = mean(diff_nonadhd_press, 'omitnan')

        avg_ratio_adhd_press    = mean(ratio_adhd_press, 'omitnan')
        avg_ratio_nonadhd_press = mean(ratio_nonadhd_press, 'omitnan')
        
    elseif strcmp(average, 'Median')    
        avg_diff_adhd_eye    = median(diff_adhd_eye, 'omitnan')
        avg_diff_nonadhd_eye = median(diff_nonadhd_eye, 'omitnan')

        avg_ratio_adhd_eye    = median(ratio_adhd_eye, 'omitnan')
        avg_ratio_nonadhd_eye = median(ratio_nonadhd_eye, 'omitnan')

        avg_diff_adhd_press    = median(diff_adhd_press, 'omitnan')
        avg_diff_nonadhd_press = median(diff_nonadhd_press, 'omitnan')

        avg_ratio_adhd_press    = median(ratio_adhd_press, 'omitnan')
        avg_ratio_nonadhd_press = median(ratio_nonadhd_press, 'omitnan')
    else
        disp('use Mean or Median')
    end
        
    sem_diff_adhd_eye = std(diff_adhd_eye, 'omitnan') ./ sqrt(length(diff_adhd_eye))
    sem_diff_nonadhd_eye = std(diff_nonadhd_eye, 'omitnan') ./ sqrt(length(diff_nonadhd_eye))
        
    sem_ratio_adhd_eye = std(ratio_adhd_eye, 'omitnan') ./ sqrt(length(ratio_adhd_eye))
    sem_ratio_nonadhd_eye = std(ratio_nonadhd_eye, 'omitnan') ./ sqrt(length(ratio_nonadhd_eye))
    
    sem_diff_adhd_press = std(diff_adhd_press, 'omitnan') ./ sqrt(length(diff_adhd_press))
    sem_diff_nonadhd_press = std(diff_nonadhd_press, 'omitnan') ./ sqrt(length(diff_nonadhd_press))

    sem_ratio_adhd_press = std(ratio_adhd_press, 'omitnan') ./ sqrt(length(ratio_adhd_press))
    sem_ratio_nonadhd_press = std(ratio_nonadhd_press, 'omitnan') ./ sqrt(length(ratio_nonadhd_press))

    %% standard plot
    plotADHDnonADHDDiffTwoValues([average, ' confusion for Gaze and Button-Press RT'],... % sgtitle
                        [diff_adhd_eye,     diff_adhd_press,     ratio_adhd_eye,     ratio_adhd_press],...
                        [avg_diff_adhd_eye, avg_diff_adhd_press, avg_ratio_adhd_eye, avg_ratio_adhd_press],...
                        [sem_diff_adhd_eye, sem_diff_adhd_press, sem_ratio_adhd_eye, sem_ratio_adhd_press], ...
                        "ADHD" , "northeast", ...  % adhd data
                        [diff_nonadhd_eye,     diff_nonadhd_press,     ratio_nonadhd_eye,     ratio_nonadhd_press], ....
                        [avg_diff_nonadhd_eye, avg_diff_nonadhd_press, avg_ratio_nonadhd_eye, avg_ratio_nonadhd_press],...
                        [sem_diff_nonadhd_eye, sem_diff_nonadhd_press, sem_ratio_nonadhd_eye, sem_ratio_nonadhd_press], ...                        
                        "nonADHD" , "northeast", ...  % nonadhd data
                        ids,...     % all participant ids
                        {"EyeDiff","PressDiff","EyeRatio","PressRatio"}, "Diff/Ratio (s)", ...% x, y axis labels "Diff: RTa - RTb (s) | Ratio: RTa/RTb (s)"
                        group_labels, {'1', '2', '3', '4'}, color_map_other, participant_map, ...
                        fullfile(comparison_results_folder, ['12_confusion_allinone_', average, '.png']), [100, 100, 1100, 500]);
    
    plotADHDnonADHDDiffTwoValues([average, ' confusion for Gaze and Button-Press RT'],... % sgtitle
                        [ratio_adhd_eye,     ratio_adhd_press],...
                        [avg_ratio_adhd_eye, avg_ratio_adhd_press],...
                        [sem_ratio_adhd_eye, sem_ratio_adhd_press], ...
                        "ADHD" , "northeast", ...  % adhd data
                        [ratio_nonadhd_eye,     ratio_nonadhd_press], ....
                        [avg_ratio_nonadhd_eye, avg_ratio_nonadhd_press],...
                        [sem_ratio_nonadhd_eye, sem_ratio_nonadhd_press], ...                        
                        "nonADHD" , "northeast", ...  % nonadhd data
                        ids,...     % all participant ids
                        {"Gaze","Button Press"}, "Ratio RTa/RTb (s)", ...% x, y axis labels "Diff: RTa - RTb (s) | Ratio: RTa/RTb (s)"
                        group_labels, {'1',  '3'}, color_map_other, participant_map, ...
                        fullfile(comparison_results_folder, ['12_confusion_allinone_ratio', average, '.png']), [100, 100, 1000, 450]);
                    
                    
    %% Create 2x2 figure
    if strcmp(average, 'Median')
    figure;
    sgtitle(strcat('Median per Participant - ', average, ' Per Group'));
    % Subplot 1: Eye-based Difference (RTa_eye - RTb_eye)
    subplot(2, 2, 1);
    hold on;
    % Plot bars
    bar(1, avg_diff_adhd_eye, 'FaceColor', color_map('ADHD'));
    errorbar(1, avg_diff_adhd_eye, sem_diff_adhd_eye, 'k', 'LineStyle', 'none');
    bar(2, avg_diff_nonadhd_eye, 'FaceColor', color_map('nonADHD'));
    errorbar(2, avg_diff_nonadhd_eye, sem_diff_nonadhd_eye, 'k', 'LineStyle', 'none');
    xticks([1 2]);
    xticklabels({'ADHD', 'nonADHD'});
    ylabel('RTa eye - RTb eye (s)');
   % Significance stars
    [~, p_diff_eye] = ttest2(diff_adhd_eye, diff_nonadhd_eye);
    if p_diff_eye < 0.05
        y_star_pos = max([avg_diff_adhd_eye + sem_diff_adhd_eye, ...
                          avg_diff_nonadhd_eye + sem_diff_nonadhd_eye]) * 1;
        text(1.5, y_star_pos, getSignificanceStars(p_diff_eye), 'FontSize', 12, 'HorizontalAlignment', 'center');
        else
    end
    title(sprintf('Time until Gaze reached Target: Confusion (Difference)[p=%.2f]',p_diff_eye));
    hold off;

    % Subplot 2: Eye-based Ratio (RTa_eye / RTb_eye)
    subplot(2, 2, 2);
    hold on;
    % Plot bars
    bar(1, avg_ratio_adhd_eye, 'FaceColor', color_map('ADHD'));
    errorbar(1, avg_ratio_adhd_eye, sem_ratio_adhd_eye, 'k', 'LineStyle', 'none');
    bar(2, avg_ratio_nonadhd_eye, 'FaceColor', color_map('nonADHD'));
    errorbar(2, avg_ratio_nonadhd_eye, sem_ratio_nonadhd_eye, 'k', 'LineStyle', 'none');
    xticks([1 2]);
    xticklabels({'ADHD', 'nonADHD'});
    ylabel('RTa eye / RTb eye (s)');
    % Significance stars
    [~, p_ratio_eye] = ttest2(ratio_adhd_eye, ratio_nonadhd_eye);
    if p_ratio_eye < 0.05
        y_star_pos = max([avg_ratio_adhd_eye + sem_ratio_adhd_eye, ...
                          avg_ratio_nonadhd_eye + sem_ratio_nonadhd_eye]) * 1.1;
        text(1.5, y_star_pos, getSignificanceStars(p_ratio_eye), 'FontSize', 12, 'HorizontalAlignment', 'center');
    end
    title(sprintf('Time until Gaze reached Target: Confusion (Ratio)[p=%.2f]',p_ratio_eye));
    hold off;
  
    % Subplot 3: Button Press Difference (RTa - RTb)
    subplot(2, 2, 3);
    hold on;
    % Plot bars
    bar(1, avg_diff_adhd_press, 'FaceColor', color_map('ADHD'));
    errorbar(1, avg_diff_adhd_press, sem_diff_adhd_press, 'k', 'LineStyle', 'none');
    bar(2, avg_diff_nonadhd_press, 'FaceColor', color_map('nonADHD'));
    errorbar(2, avg_diff_nonadhd_press, sem_diff_nonadhd_press, 'k', 'LineStyle', 'none');
    xticks([1 2]);
    xticklabels({'ADHD', 'nonADHD'});
    ylabel('RTa - RTb (s)');

    % Significance stars
    [~, p_diff_press] = ttest2(diff_adhd_press, diff_nonadhd_press);
    if p_diff_press < 0.05
        y_star_pos = max([avg_diff_adhd_press + sem_diff_adhd_press, ...
                          avg_diff_nonadhd_press + sem_diff_nonadhd_press]) * 1.01;
        text(1.5, y_star_pos, getSignificanceStars(p_diff_press), 'FontSize', 12, 'HorizontalAlignment', 'center');
    end
    title(sprintf('Time until Button Press: Confusion (Difference)[p=%.2f]',p_diff_press));
    hold off;

    % Subplot 4: Button Press Ratio (RTa / RTb)
    subplot(2, 2, 4);
    hold on;
    % Plot bars
    bar(1, avg_ratio_adhd_press, 'FaceColor', color_map('ADHD'));
    errorbar(1, avg_ratio_adhd_press, sem_ratio_adhd_press, 'k', 'LineStyle', 'none');
    bar(2, avg_ratio_nonadhd_press, 'FaceColor', color_map('nonADHD'));
    errorbar(2, avg_ratio_nonadhd_press, sem_ratio_nonadhd_press, 'k', 'LineStyle', 'none');
    xticks([1 2]);
    xticklabels({'ADHD', 'nonADHD'});
    ylabel('RTa / RTb (s)');
    % Significance stars
    [~, p_ratio_press] = ttest2(ratio_adhd_press, ratio_nonadhd_press);
    if p_ratio_press < 0.05
        y_star_pos = max([avg_ratio_adhd_press + sem_ratio_adhd_press, ...
                          avg_ratio_nonadhd_press + sem_ratio_nonadhd_press]) * 1.1;
        text(1.5, y_star_pos, getSignificanceStars(p_ratio_press), 'FontSize', 12, 'HorizontalAlignment', 'center');
    end
    title(sprintf('Time until Button Press: Confusion (Ratio) [p=%.2f]',p_ratio_press));
    hold off;

    % Save the figure if specified
    if safe == 1
        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf, strcat(comparison_results_folder, '\12_confusion_',average, '_barplot_with_individual_lines.png'));
    end
    else
    end
end

% Helper function to determine the number of stars based on p-value
function stars = getSignificanceStars(p_value)
    if p_value < 0.001
        stars = '***';
    elseif p_value < 0.01
        stars = '**';
    elseif p_value < 0.05
        stars = '*';
    else
        stars = '';
    end
end

