function plotCosineSimilarity(trial_results, trial_metrics, saccades, fixations, screenXpixels, screenYpixels, conditions, condition_labels, color_map, safe, analysis_folder)

num_trials = length(saccades);
x_s_start = []; y_s_start = [];
x_s_end = []; y_s_end = [];
x_target = []; y_target = [];
for a=1:num_trials
    if ~isempty(saccades(a).saccades)
        x_s_start = [x_s_start; saccades(a).saccades(1).startCenter(1)];
        y_s_start = [y_s_start; saccades(a).saccades(1).startCenter(2)];
        
        x_s_end = [x_s_end; saccades(a).saccades(1).endCenter(1)];
        y_s_end = [y_s_end; saccades(a).saccades(1).endCenter(2)];
        
        x_target = [x_target; fixations.targetCenters(1, a)];
        y_target = [y_target; fixations.targetCenters(2, a)];
    else
        x_s_start = [x_s_start; NaN];
        y_s_start = [y_s_start; NaN];
        
        x_s_end = [x_s_end; NaN];
        y_s_end = [y_s_end; NaN];
        
        x_target = [x_target; NaN];
        y_target = [y_target; NaN];
    end
end

x_s_end_centered  = x_s_end - x_s_start;
y_s_end_centered  = y_s_end - y_s_start;
x_target_centered = x_target - x_s_start;
y_target_centered = y_target - y_s_start;

s_end_vector =  [x_s_end_centered(:), y_s_end_centered(:)];
target_vector = [x_target_centered(:), y_target_centered(:)];


trialConditions = trial_results.Condition;
cos_sims_cond = zeros(length(conditions), length(trialConditions)/4);


%% per condition all trials
figure; grid off;
for c = 1:length(conditions)
    subplot(2, 2, c); hold on;
    condition_idx = strcmp(trialConditions, conditions{c});
    
    s_end_current  = s_end_vector(condition_idx, :);
    target_current = target_vector(condition_idx, :);
    
    cosSim_current = dot(s_end_current, target_current, 2) ./ (vecnorm(s_end_current, 2, 2) .* vecnorm(target_current, 2, 2));
    %participant_means(c) = nanmean(cosSim_current);
    cos_sims_cond(c, :) = cosSim_current;
    
    bar(c, cosSim_current, 'FaceColor', color_map(conditions{c}), 'DisplayName', conditions{c});
    xlabel('Trial'); ylabel('Cosine Similarity');
    xticks('');
    title(condition_labels{c});
end
hold off;
sgtitle(strcat('ID ' , num2str(trial_metrics.id), ': Cosine Similarity Between Actual and Optimal Paths'));

if safe == 1
    set(gcf, 'Position', [50, 50, 1000, 700]); % Resize the figure window (x, y, width, height)
    %set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf, fullfile(analysis_folder, 'cosine_similarity.png'));
end


%% cosine similarity, scatter plot, per condition, mean, median, sem, lines
mean_values = nan(1, length(conditions));
median_values = nan(1, length(conditions));
sem_values = nan(1, length(conditions));

figure;hold on;grid on;
for c = 1:length(conditions)
    valid_data = cos_sims_cond(c, ~isnan(cos_sims_cond(c, :)));
    x_jitter = c + 0.1 * (rand(size(valid_data)) - 0.5);
    scatter(x_jitter, valid_data, 30, 'filled', ...
        'MarkerFaceColor', color_map(conditions{c}), 'MarkerFaceAlpha', 0.5, ...
        'MarkerEdgeColor', color_map(conditions{c}), 'MarkerEdgeAlpha', 0.5);
    
    mean_values(c) = nanmean(valid_data);
    median_values(c) = nanmedian(valid_data);
    
    sem_values(c) = std(valid_data);% / sqrt(length(valid_data)); % Standard error of the mean
    
    errorbar(c, mean_values(c), sem_values(c), 'Color', 'k', ...
        'LineWidth', 1.5, 'CapSize', 10); % Error bar in matching color
end

% Plot mean and median lines
mean_line = plot(1:length(conditions), mean_values, '-', 'LineWidth', 2, 'Color', [0 0 0 0.9], ...
    'DisplayName', 'Mean');
median_line = plot(1:length(conditions), median_values, '--', 'LineWidth', 2, 'Color', [0 0 0 0.6], ...
    'DisplayName', 'Median');
% Scatter plot for means
scatter(1:length(conditions), mean_values, 70, 'o', ...
    'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.9, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.9, ...
    'DisplayName', 'Mean');

% Scatter plot for medians
scatter(1:length(conditions), median_values, 70, 'o', ...
    'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.6, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.6, ...
    'DisplayName', 'Median');

ylabel('Cosine Similarity');
title(strcat('ID ' , num2str(trial_metrics.id), ': Cosine Similarity by Condition with individual Trial Values'));
xticks(1:length(conditions));xticklabels(condition_labels);xlim([0.5, length(conditions) + 0.5]);
%ylim([-1.02, 1.02]);
legend([mean_line, median_line], {'Mean', 'Median'}, 'Location', 'northwest');
legend('boxoff')

if safe == 1
    set(gcf, 'Position', [50, 50, 1000, 700]);
    saveas(gcf, fullfile(analysis_folder, 'cosine_similarity_colored_scatter_mean_median_sem.png'));
    
    
    cosSim_all = dot(s_end_vector, target_vector, 2) ./ (vecnorm(s_end_vector, 2, 2) .* vecnorm(target_vector, 2, 2));
    cosSim.mean_values = mean_values; 
    cosSim.median_values = median_values; 
    cosSim.sem_values = sem_values; 
    cosSim.cosSim = cosSim_all;
    save(fullfile(analysis_folder, '\cosSim.mat'), 'cosSim');
end





%%  quiver plot for actual + optimal path in one
% cosSim_all = dot(s_end_vector, target_vector, 2) ./ (vecnorm(s_end_vector, 2, 2) .* vecnorm(target_vector, 2, 2));
% figure;        hold on;        grid on;
% for a = 1:num_trials
%     if ~isnan(cosSim_all(a))
%         % Plot start, end, and target points
%         scatter([x_s_start(a), x_s_end(a), x_target(a)], [y_s_start(a), y_s_end(a), y_target(a)], 100, 'filled');
%         % Plot vectors
%         quiver(x_s_start(a), y_s_start(a), x_s_end_centered(a), y_s_end_centered(a), 0, 'r', 'LineWidth', 1.5); % Actual path
%         quiver(x_s_start(a), y_s_start(a), x_target_centered(a), y_target_centered(a), 0, 'g', 'LineWidth', 1.5); % Optimal path
%         % Label cosine similarity
%         text(x_s_end(a), y_s_end(a), sprintf('%.2f', cosSim_all(a)), 'Color', 'b', 'FontSize', 10);
%     end
% end
% xlabel('screen x axis pixels');  ylabel('screen y axis pixels');
% title(strcat('ID ' , num2str(trial_metrics.id), ': Saccade Paths and Cosine Similarity to Optimal Paths'));
% legend('Points', 'Actual Path', 'Optimal Path', 'Location', 'best');


%% quiver two plots one optimal one actual
figure; %t = tiledlayout(1, 2, 'TileSpacing', 'Compact');
colors = lines(num_trials);  % Generates a colormap with as many colors as trials
% Tile 1: Actual saccade data
nexttile; hold on;
xlim([0, screenXpixels]);
ylim([0, screenYpixels]);
for i = 1:num_trials
    quiver(x_s_start(i), y_s_start(i), x_s_end(i)-x_s_start(i), y_s_end(i)-y_s_start(i), 0, 'Color', colors(i, :), 'LineWidth', 1.5); % Plot vector
end
title('First Saccade Per Trial Vectors');
xlabel('screen x pixels');ylabel('screen y pixels');hold off;
% Tile 2: Optimal saccade data
nexttile;hold on;
xlim([0, screenXpixels]);
ylim([0, screenYpixels]);
for i = 1:num_trials
    quiver(x_s_start(i), y_s_start(i), x_target(i)-x_s_start(i), y_target(i)-y_s_start(i), 0, 'Color', colors(i, :), 'LineWidth', 1.5); % Plot vector
end
title('Optimal Saccade Vectors');
xlabel('screen x pixels');ylabel('screen y pixels');hold off;
sgtitle(strcat('ID ' , num2str(trial_metrics.id), ': Comparison of Actual and Optimal Saccade Vectors'));
if safe == 1
    set(gcf, 'Position', [50, 50, 1000, 700]);
    saveas(gcf, fullfile(analysis_folder, 'cosine_similarity_quiver_2_plots.png'));
end

%%  quiver two plots one optimal one actual - 0 centered
figure; 
colors = lines(num_trials);  % Generates a colormap with as many colors as trials
% Tile 1: Actual saccade data
nexttile; hold on; 
xlim([-screenXpixels/2, screenXpixels/2]);
ylim([-screenYpixels/2, screenYpixels/2]);
for i = 1:num_trials
    quiver(0, 0, x_s_end(i)-x_s_start(i), y_s_end(i)-y_s_start(i), 0, 'Color', colors(i, :), 'LineWidth', 1.5); % Plot vector
end
title('Actual First Saccade per Trial Vectors');xlabel('screen x pixels');ylabel('screen y pixels');hold off;
% Tile 2: Optimal saccade data
nexttile;hold on; 
xlim([-screenXpixels/2, screenXpixels/2]);
ylim([-screenYpixels/2, screenYpixels/2]);
for i = 1:num_trials
    quiver(0, 0, x_target(i)-x_s_start(i), y_target(i)-y_s_start(i), 0, 'Color', colors(i, :), 'LineWidth', 1.5); % Plot vector
end
title('Optimal Saccade Vectors');xlabel('screen x pixels');ylabel('screen y pixels');hold off;

sgtitle(strcat('ID ' , num2str(trial_metrics.id), ': Comparison of Actual and Optimal Saccade Vectors'));


if safe == 1
    set(gcf, 'Position', [50, 50, 1000, 700]);
    saveas(gcf, fullfile(analysis_folder, 'cosine_similarity_quiver_2_plots_0_centered.png'));
end
