function plotConfusion(data, color_map, comparison_results_folder, safe)
    % Initialize arrays to hold the differences and ratios
    diff_adhd = [];
    diff_nonadhd = [];
    ratio_adhd = [];
    ratio_nonadhd = [];

    % Loop through each participant and calculate differences and ratios
    for i = 1:length(data)
        if strcmp(data(i).group, 'ADHD')
            diff_adhd = [diff_adhd; data(i).RTa - data(i).RTb];
            ratio_adhd = [ratio_adhd; data(i).RTa / data(i).RTb];
        elseif strcmp(data(i).group, 'nonADHD')
            diff_nonadhd = [diff_nonadhd; data(i).RTa - data(i).RTb];
            ratio_nonadhd = [ratio_nonadhd; data(i).RTa / data(i).RTb];
        end
    end

    % Calculate means and SEMs
    mean_diff_adhd = mean(diff_adhd, 'omitnan');
    mean_diff_nonadhd = mean(diff_nonadhd, 'omitnan');
    sem_diff_adhd = std(diff_adhd, 'omitnan') / sqrt(length(diff_adhd));
    sem_diff_nonadhd = std(diff_nonadhd, 'omitnan') / sqrt(length(diff_nonadhd));

    mean_ratio_adhd = mean(ratio_adhd, 'omitnan');
    mean_ratio_nonadhd = mean(ratio_nonadhd, 'omitnan');
    sem_ratio_adhd = std(ratio_adhd, 'omitnan') / sqrt(length(ratio_adhd));
    sem_ratio_nonadhd = std(ratio_nonadhd, 'omitnan') / sqrt(length(ratio_nonadhd));

    % Create figure with two subplots
    figure;
    
    % Subplot 1: Difference (RTa - RTb)
    subplot(1, 2, 1);
    hold on;
    bar(1, mean_diff_adhd, 'FaceColor', color_map('ADHD'));
    errorbar(1, mean_diff_adhd, sem_diff_adhd, 'k', 'LineStyle', 'none');
    bar(2, mean_diff_nonadhd, 'FaceColor', color_map('nonADHD'));
    errorbar(2, mean_diff_nonadhd, sem_diff_nonadhd, 'k', 'LineStyle', 'none');
    xticks([1 2]);
    xticklabels({'ADHD', 'nonADHD'});
    ylabel('RTa - RTb (s)');
    title('Confusion (Difference)');
    hold off;
    grid on;
    
    % Subplot 2: Ratio (RTa / RTb)
    subplot(1, 2, 2);
    hold on;
    bar(1, mean_ratio_adhd, 'FaceColor', color_map('ADHD'));
    errorbar(1, mean_ratio_adhd, sem_ratio_adhd, 'k', 'LineStyle', 'none');
    bar(2, mean_ratio_nonadhd, 'FaceColor', color_map('nonADHD'));
    errorbar(2, mean_ratio_nonadhd, sem_ratio_nonadhd, 'k', 'LineStyle', 'none');
    xticks([1 2]);
    xticklabels({'ADHD', 'nonADHD'});
    ylabel('RTa / RTb');
    title('Confusion (Ratio)');
    hold off;
    grid on;
    
    % Save the figure if safe is set to 1
    if safe == 1
        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf, fullfile(comparison_results_folder, 'bar_mean_confusion.png'));
    end
end
