function plotMeanRTButtonPressVsEyecomparison(data, color_map, conditions, comparison_results_folder, safe)
    % Initialize empty arrays to hold the overall mean RTs and SEMs per condition for each group
    avg_rt_adhd = [];
    avg_rt_eye_adhd = [];
    avg_rt_diff_adhd = [];
    
    avg_rt_nonadhd = [];
    avg_rt_eye_nonadhd = [];
    avg_rt_diff_nonadhd = [];
    
    avg_nrt_adhd = [];
    avg_nrt_eye_adhd = [];
    avg_nrt_diff_adhd = [];
    
    avg_nrt_nonadhd = [];
    avg_nrt_eye_nonadhd = [];
    avg_nrt_diff_nonadhd = [];
    
    % Loop through each participant and aggregate the mean RT and accuracy per condition
    for i = 1:length(data)
        if strcmp(data(i).group, 'ADHD')
            avg_rt_adhd      = [avg_rt_adhd; data(i).RTa, data(i).RTasimple, data(i).RTb, data(i).RTbsimple];
            avg_rt_eye_adhd  = [avg_rt_eye_adhd; data(i).RTa_eye, data(i).RTasimple_eye, data(i).RTb_eye, data(i).RTbsimple_eye];
            avg_rt_diff_adhd = [avg_rt_diff_adhd; data(i).nRTa - data(i).nRTa_eye, data(i).nRTb - data(i).nRTb_eye, ...
                                                  data(i).nRTasimple - data(i).nRTasimple_eye, data(i).nRTbsimple - data(i).nRTbsimple_eye];
            
            avg_nrt_adhd      = [avg_nrt_adhd; data(i).nRTa, data(i).nRTasimple, data(i).nRTb, data(i).nRTbsimple];
            avg_nrt_eye_adhd  = [avg_nrt_eye_adhd; data(i).nRTa_eye, data(i).nRTasimple_eye, data(i).nRTb_eye, data(i).nRTbsimple_eye];
            avg_nrt_diff_adhd = [avg_nrt_diff_adhd; data(i).nRTa - data(i).nRTa_eye, data(i).nRTb - data(i).nRTb_eye, ...
                                                    data(i).nRTasimple - data(i).nRTasimple_eye, data(i).nRTbsimple - data(i).nRTbsimple_eye];
        elseif strcmp(data(i).group, 'nonADHD')
            avg_rt_nonadhd      = [avg_rt_nonadhd; data(i).RTa, data(i).RTasimple, data(i).RTb, data(i).RTbsimple];
            avg_rt_eye_nonadhd  = [avg_rt_eye_nonadhd; data(i).RTa_eye, data(i).RTasimple_eye, data(i).RTb_eye, data(i).RTbsimple_eye];
            avg_rt_diff_nonadhd = [avg_rt_diff_nonadhd; data(i).nRTa - data(i).nRTa_eye, data(i).nRTb - data(i).nRTb_eye, ...
                                                        data(i).nRTasimple - data(i).nRTasimple_eye, data(i).nRTbsimple - data(i).nRTbsimple_eye];
            
            avg_nrt_nonadhd      = [avg_nrt_nonadhd; data(i).nRTa, data(i).nRTasimple, data(i).nRTb, data(i).nRTbsimple];
            avg_nrt_eye_nonadhd  = [avg_nrt_eye_nonadhd; data(i).nRTa_eye, data(i).nRTasimple_eye, data(i).nRTb_eye, data(i).nRTbsimple_eye];
            avg_nrt_diff_nonadhd = [avg_nrt_diff_nonadhd; data(i).nRTa - data(i).nRTa_eye, data(i).nRTb - data(i).nRTb_eye, ...
                                                          data(i).nRTasimple - data(i).nRTasimple_eye, data(i).nRTbsimple - data(i).nRTbsimple_eye];
        end
    end
    
    % Calculate the overall median and SEM for each condition and group
    overall_mean_rt_adhd      = median(avg_rt_adhd, 1, 'omitnan');
    overall_mean_rt_eye_adhd  = median(avg_rt_eye_adhd, 1, 'omitnan');
    overall_mean_rt_diff_adhd = median(avg_rt_diff_adhd, 1, 'omitnan');
    
    overall_mean_rt_nonadhd      = median(avg_rt_nonadhd, 1, 'omitnan');
    overall_mean_rt_eye_nonadhd  = median(avg_rt_eye_nonadhd, 1, 'omitnan');
    overall_mean_rt_diff_nonadhd = median(avg_rt_diff_nonadhd, 1, 'omitnan');
    
    overall_sem_rt_adhd      = std(avg_rt_adhd, 0, 1, 'omitnan') / sqrt(size(avg_rt_adhd, 1));
    overall_sem_rt_eye_adhd  = std(avg_rt_eye_adhd, 0, 1, 'omitnan') / sqrt(size(avg_rt_eye_adhd, 1));
    overall_sem_rt_diff_adhd = std(avg_rt_diff_adhd, 0, 1, 'omitnan') / sqrt(size(avg_rt_diff_adhd, 1));
    
    overall_sem_rt_nonadhd      = std(avg_rt_nonadhd, 0, 1, 'omitnan') / sqrt(size(avg_rt_nonadhd, 1));
    overall_sem_rt_eye_nonadhd  = std(avg_rt_eye_nonadhd, 0, 1, 'omitnan') / sqrt(size(avg_rt_eye_nonadhd, 1));
    overall_sem_rt_diff_nonadhd = std(avg_rt_diff_nonadhd, 0, 1, 'omitnan') / sqrt(size(avg_rt_diff_nonadhd, 1));
    
    % Repeat for normalized data
    overall_mean_nrt_adhd = median(avg_nrt_adhd, 1, 'omitnan');
    overall_mean_nrt_eye_adhd = median(avg_nrt_eye_adhd, 1, 'omitnan');
    overall_mean_nrt_diff_adhd = median(avg_nrt_diff_adhd, 1, 'omitnan');
    
    overall_mean_nrt_nonadhd = median(avg_nrt_nonadhd, 1, 'omitnan');
    overall_mean_nrt_eye_nonadhd = median(avg_nrt_eye_nonadhd, 1, 'omitnan');
    overall_mean_nrt_diff_nonadhd = median(avg_nrt_diff_nonadhd, 1, 'omitnan');
    
    overall_sem_nrt_adhd = std(avg_nrt_adhd, 0, 1, 'omitnan') / sqrt(size(avg_nrt_adhd, 1));
    overall_sem_nrt_eye_adhd = std(avg_nrt_eye_adhd, 0, 1, 'omitnan') / sqrt(size(avg_nrt_eye_adhd, 1));
    overall_sem_nrt_diff_adhd = std(avg_nrt_diff_adhd, 0, 1, 'omitnan') / sqrt(size(avg_nrt_diff_adhd, 1));
    
    overall_sem_nrt_nonadhd = std(avg_nrt_nonadhd, 0, 1, 'omitnan') / sqrt(size(avg_nrt_nonadhd, 1));
    overall_sem_nrt_eye_nonadhd = std(avg_nrt_eye_nonadhd, 0, 1, 'omitnan') / sqrt(size(avg_nrt_eye_nonadhd, 1));
    overall_sem_nrt_diff_nonadhd = std(avg_nrt_diff_nonadhd, 0, 1, 'omitnan') / sqrt(size(avg_nrt_diff_nonadhd, 1));
    
    % Create a 3x2 subplot
    figure;
    
    % Subplot 1: Mean Button Press RTs (Non-normalized)
    subplot(3, 2, 1);
    hold on;
    errorbar(1:length(conditions), overall_mean_rt_adhd, overall_sem_rt_adhd, '-o', 'Color', color_map('ADHD'), ...
        'DisplayName', 'ADHD', 'LineWidth', 2, 'MarkerSize', 8);
    errorbar(1:length(conditions), overall_mean_rt_nonadhd, overall_sem_rt_nonadhd, '-o', 'Color', color_map('nonADHD'), ...
        'DisplayName', 'nonADHD', 'LineWidth', 2, 'MarkerSize', 8);
    hold off;
    xticks(1:length(conditions));
    xticklabels(conditions);
    ylabel('Median Button Press RT (s)');
    title('Button Press RT');
    grid on;
    legend('Orientation', 'horizontal');
    
    % Subplot 2: Normalized Button Press RTs
    subplot(3, 2, 2);
    hold on;
    errorbar(1:length(conditions), overall_mean_nrt_adhd, overall_sem_nrt_adhd, '-o', 'Color', color_map('ADHD'), ...
        'DisplayName', 'ADHD', 'LineWidth', 2, 'MarkerSize', 8);
    errorbar(1:length(conditions), overall_mean_nrt_nonadhd, overall_sem_nrt_nonadhd, '-o', 'Color', color_map('nonADHD'), ...
        'DisplayName', 'nonADHD', 'LineWidth', 2, 'MarkerSize', 8);
    hold off;
    xticks(1:length(conditions));
    xticklabels(conditions);
    ylabel('Normalized Button Press RT (s)');
    title('Normalized Button Press RT');
    grid on;
    
    % Subplot 3: Mean Eye RTs (Non-normalized)
    subplot(3, 2, 3);
    hold on;
    errorbar(1:length(conditions), overall_mean_rt_eye_adhd, overall_sem_rt_eye_adhd, '-o', 'Color', color_map('ADHD'), ...
        'DisplayName', 'ADHD', 'LineWidth', 2, 'MarkerSize', 8);
    errorbar(1:length(conditions), overall_mean_rt_eye_nonadhd, overall_sem_rt_eye_nonadhd, '-o', 'Color', color_map('nonADHD'), ...
        'DisplayName', 'nonADHD', 'LineWidth', 2, 'MarkerSize', 8);
    hold off;
    xticks(1:length(conditions));
    xticklabels(conditions);
    ylabel('Median Eye RT (s)');
    title('Eye Movement RT');
    grid on;
    
    % Subplot 4: Normalized Eye RTs
    subplot(3, 2, 4);
    hold on;
    errorbar(1:length(conditions), overall_mean_nrt_eye_adhd, overall_sem_nrt_eye_adhd, '-o', 'Color', color_map('ADHD'), ...
        'DisplayName', 'ADHD', 'LineWidth', 2, 'MarkerSize', 8);
    errorbar(1:length(conditions), overall_mean_nrt_eye_nonadhd, overall_sem_nrt_eye_nonadhd, '-o', 'Color', color_map('nonADHD'), ...
        'DisplayName', 'nonADHD', 'LineWidth', 2, 'MarkerSize', 8);
    hold off;
    xticks(1:length(conditions));
    xticklabels(conditions);
    ylabel('Normalized Eye RT (s)');
    title('Normalized Eye Movement RT');
    grid on;
    
    % Subplot 5: Mean Difference Between Button Press RT and Eye RT (Non-normalized)
    subplot(3, 2, 5);
    hold on;
    errorbar(1:length(conditions), overall_mean_rt_diff_adhd, overall_sem_rt_diff_adhd, '-o', 'Color', color_map('ADHD'), ...
        'DisplayName', 'ADHD', 'LineWidth', 2, 'MarkerSize', 8);
    errorbar(1:length(conditions), overall_mean_rt_diff_nonadhd, overall_sem_rt_diff_nonadhd, '-o', 'Color', color_map('nonADHD'), ...
        'DisplayName', 'nonADHD', 'LineWidth', 2, 'MarkerSize', 8);
    hold off;
    xticks(1:length(conditions));
    xticklabels(conditions);
    ylabel('Median Difference RT (s)');
    title('Difference: Button Press RT - Eye RT');
    grid on;
    
    % Subplot 6: Normalized Difference Between Button Press RT and Eye RT
    subplot(3, 2, 6);
    hold on;
    errorbar(1:length(conditions), overall_mean_nrt_diff_adhd, overall_sem_nrt_diff_adhd, '-o', 'Color', color_map('ADHD'), ...
        'DisplayName', 'ADHD', 'LineWidth', 2, 'MarkerSize', 8);
    errorbar(1:length(conditions), overall_mean_nrt_diff_nonadhd, overall_sem_nrt_diff_nonadhd, '-o', 'Color', color_map('nonADHD'), ...
        'DisplayName', 'nonADHD', 'LineWidth', 2, 'MarkerSize', 8);
    hold off;
    xticks(1:length(conditions));
    xticklabels(conditions);
    ylabel('Normalized Difference RT (s)');
    title('Normalized Difference: Button Press RT - Eye RT');
    grid on;
    
    % Manually set the y-axis limits to be the same for corresponding plots
    % Get all subplot handles
    h1 = subplot(3, 2, 1); h2 = subplot(3, 2, 2);
    h3 = subplot(3, 2, 3); h4 = subplot(3, 2, 4);
    h5 = subplot(3, 2, 5); h6 = subplot(3, 2, 6);
    
    % Set y-limits for Button Press RTs (Subplot 1 and 2)
    all_ys_button_press = [get(h1, 'YLim'), get(h2, 'YLim')];
    set([h1, h2], 'YLim', [min(all_ys_button_press), max(all_ys_button_press)]);
    
    % Set y-limits for Eye RTs (Subplot 3 and 4)
    all_ys_eye = [get(h3, 'YLim'), get(h4, 'YLim')];
    set([h3, h4], 'YLim', [min(all_ys_eye), max(all_ys_eye)]);
    
    % Set y-limits for Difference RTs (Subplot 5 and 6)
    all_ys_diff = [get(h5, 'YLim'), get(h6, 'YLim')];
    set([h5, h6], 'YLim', [min(all_ys_diff), max(all_ys_diff)]);
    
    % Save the figure if safe is set to 1
    if safe == 1
        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf, fullfile(comparison_results_folder, 'mean_of_means_RT_per_condition_group.png'));
    end
end
