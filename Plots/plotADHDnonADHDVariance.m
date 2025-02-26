function plotADHDnonADHDVariance(bigtitle,... % sgtitle
                        data_1, group_1_medians, group_1_sems, label_1, legend_1_place, ...  % adhd data
                        data_2, group_2_medians, group_2_sems, label_2, legend_2_place, ...  % nonadhd data
                        ids,...     % all participant ids
                        x_labels, y_label, ...% x, y axis labels
                        group_labels, conditions, color_map, participant_map, ...
                        safe_name)
% data_1&2 : [participants x conditions]
            
                    
ymax = max(max(data_1(:)), max(data_2(:)))+ max(max(group_1_sems(:)), max(group_2_sems(:))); % Maximum individual median value
ymin = min(min(data_1(:)), min(data_2(:)))- max(max(group_1_sems(:)), max(group_2_sems(:)))/2; % Maximum individual median value

sharedYLimits = [ymin, ymax];


%%
figure; sgtitle(bigtitle);
% ADHD Plot
subplot(1, 3, 1);
plotStandardTwoGroupsEachData(data_1, group_1_medians, group_1_sems,label_1, y_label, legend_1_place, ids, group_labels, conditions, x_labels, color_map, participant_map, sharedYLimits);

% nonADHD Plot
subplot(1, 3, 2);
plotStandardTwoGroupsEachData(data_2, group_2_medians, group_2_sems, label_2, y_label, legend_2_place, ids, group_labels, conditions, x_labels, color_map, participant_map, sharedYLimits);

% Variance Difference Plot
subplot(1, 3, 3); 
plotVarianceDifference(data_1, data_2, conditions, x_labels, color_map);
  
set(gcf, 'Position', [50, 50, 1400, 700]); % Resize the figure window (x, y, width, height)
saveas(gcf, safe_name);


%% 
% figure; sgtitle(bigtitle);
% subplot(1, 2, 1);
% plotGroups(data_1, label_1, ...
%            data_2, label_2, ...
%            fixationStats, group_labels, conditions, color_map, participant_map, x_labels)
% subplot(1, 2, 2);
% plotVarianceDifference(data_1, data_2, conditions, x_labels, color_map);
% set(gcf, 'Position', [50, 50, 1400, 700]); % Resize the figure window (x, y, width, height)
% saveas(gcf, safe_name);

end



function plotStandardTwoGroupsEachData(data, medians, stds, group, y_label, position, ids, group_labels, conditions, labels, color_map, participant_map, sharedYLimits)
    % Plot group medians and scatter individual participant data
    numConditions = length(conditions);
    hold on;
    
    % Bar plot for medians with error bars
    for c = 1:numConditions
        bar(c, medians(c), 'FaceColor', color_map(conditions{c}), 'EdgeColor', 'none', ...
            'FaceAlpha', 0.3, 'DisplayName', group, 'HandleVisibility', 'off');
    end
    
    current_indices = strcmp(group_labels, group);
    group_ids = ids(current_indices);

    % Scatter individual participant data
    jitterAmount = 0.4; 
    for participant = 1:length(group_ids)
        participant_id = group_ids(participant);

        %if strcmp(group_labels{participant}, group)
            info = participant_map(participant_id); % to retrieve color and marker
            current_data = data(participant,:);
            jitteredX = (1:numConditions) + jitterAmount * (rand(1, numConditions) - 0.5); % Jitter x-coordinates
            
            scatter(jitteredX, current_data, 40, ... % Use scatter for points only
                'MarkerEdgeColor', info.color, ...
                'MarkerFaceColor', info.color, ...
                'Marker', info.marker, ...
                'DisplayName', sprintf('ID%d', participant_id));
        %end
    end
    
    for c = 1:numConditions % errorbar on top for better visibility
       errorbar(c, medians(c), stds(c), '.', 'LineWidth', 1.5, 'CapSize', 10, 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    end
    xticks(1:numConditions); xticklabels(labels);
    ylabel(y_label);  ylim(sharedYLimits);
    legend('boxoff');
    legend('Location', position, 'Orientation', 'horizontal', 'NumColumns', 2); % Two-column legend
    title(group);    hold off;
end

%%
function plotVarianceDifference(adhd_medians, nonadhd_medians, conditions, labels, color_map)
    % Calculate variance metric and perform T-tests or permutation tests
    
    % Calculate variances for ADHD and nonADHD groups
    adhd_variances = var(adhd_medians, 0, 1, 'omitnan'); % Variance across participants per condition
    nonadhd_variances = var(nonadhd_medians, 0, 1, 'omitnan');
    adhd_sd = std(adhd_medians, 0, 1, 'omitnan');
    nonadhd_sd = std(nonadhd_medians, 0, 1, 'omitnan');

    % Calculate variance metric: Sq((Variance 1 + Variance 2) / 2)
    % square root((var 1 var 2) / std deviation )
    
    %% Variances Assumed Unequal: Cohen’s d(av)
    % square root of the average variance.
    % https://cran.r-project.org/web/packages/TOSTER/vignettes/SMD_calcs.html
    sq_avg_var = sqrt((adhd_variances + nonadhd_variances) / 2);
    final = (adhd_sd - nonadhd_sd) ./ sq_avg_var;
    
    %% Variances Assumed Equal: Cohen’s d
    % the denominator is simply the pooled standard deviation.
    % https://cran.r-project.org/web/packages/TOSTER/vignettes/SMD_calcs.html
    % todo 
    
    % The standardized mean difference
    % https://handbook-5-1.cochrane.org/chapter_9/9_2_3_2_the_standardized_mean_difference.htm
    
    % Perform statistical tests and store p-values
    numConditions = length(conditions);
    p_values = nan(1, numConditions); % Store p-values for annotation
    hold on;
    for c = 1:numConditions
        % Extract medians for each condition
        adhd_data = adhd_medians(:, c); % ADHD group medians for condition c
        nonadhd_data = nonadhd_medians(:, c); % nonADHD group medians for condition c
        
        % Perform a T-test (or permutation test if needed)
        [~, p] = ttest2(adhd_data, nonadhd_data);
        p_values(c) = p; % Store p-value
        
        % Plot variance metric
        b = bar(c, final(c), 'FaceColor', color_map(conditions{c}), 'EdgeColor', 'none', ...
            'FaceAlpha', 0.3, 'DisplayName', sprintf('%s Variance', conditions{c}));
        
        xtips1 = b(1).XEndPoints;
        ytips1 = b(1).YEndPoints;
        text(xtips1,ytips1, sprintf('p=%.2f', p),'HorizontalAlignment','center',...
                'VerticalAlignment','bottom')
    end

    % Customize plot
    xticks(1:numConditions); xticklabels(labels);
    ylabel('Normalized Difference');
    %legend('Location', 'Best'); legend('boxoff');
    title('Difference in Variance');
    hold off;
end





function plotGroups(data_1, label_1, ...
                    data_2, label_2, ...
                    fixationStats, group_labels, conditions, color_map, participant_map, condition_labels)
       
                
group_1_medians = median(data_1, 1, 'omitnan');
group_1_stds = std(data_1, 0, 1, 'omitnan');
group_2_medians = median(data_2, 1, 'omitnan');
group_2_stds = std(data_2, 0, 1, 'omitnan');

bar_width = 0.4; 
offset = 0.2;    % Offset between ADHD and nonADHD bars
x_positions = 1:size(conditions,2);

positions_1 = x_positions - offset;% Adjusted positions for ADHD and nonADHD bars
positions_2 = x_positions + offset;

hold on;% Overlay bars for ADHD group
bar1 = bar(positions_1, group_1_medians, bar_width, 'FaceColor', color_map(label_1), ...
           'FaceAlpha', 0.5, 'EdgeColor', 'none');
errorbar(positions_1, group_1_medians, group_1_stds, '.', 'LineWidth', 1.5, 'CapSize', 10, 'Color', '#609860');

% Overlay bars for nonADHD group
bar2 = bar(positions_2, group_2_medians, bar_width, 'FaceColor', color_map(label_2), ...
           'FaceAlpha', 0.5, 'EdgeColor', 'none');
errorbar(positions_2, group_2_medians, group_2_stds, '.', 'LineWidth', 1.5, 'CapSize', 10, 'Color', '#9F3030');


% Scatter individual participant data
    jitterAmount = 0.3; 
    for participant = 1:size(participant_map, 1)
        participant_id = fixationStats(participant).id;
        info = participant_map(participant_id); % to retrieve color and marker

        if strcmp(group_labels{participant}, label_1)
            jitteredX = positions_1 + jitterAmount * (rand(1, size(conditions,2)) - 0.5); % Jitter for Group 1
           
            data = [fixationStats(participant).a_median, fixationStats(participant).b_median, ...
                    fixationStats(participant).as_median, fixationStats(participant).bs_median];
            scatter(jitteredX, data, 50, ... % Use scatter for points only
                'MarkerEdgeColor', info.color, ...
                'MarkerFaceColor', info.color, ...
                'Marker', info.marker, ...
                'DisplayName', sprintf('ID%d', participant_id));
        else
            jitteredX = positions_2 + jitterAmount * (rand(1, size(conditions,2)) - 0.5); % Jitter for Group 2  
            data = [fixationStats(participant).a_median, fixationStats(participant).b_median, ...
                    fixationStats(participant).as_median, fixationStats(participant).bs_median];
            scatter(jitteredX, data, 50, ... % Use scatter for points only
                'MarkerEdgeColor', info.color, ...
                'MarkerFaceColor', info.color, ...
                'Marker', info.marker, ...
                'DisplayName', sprintf('ID%d', participant_id));
        end
    end

xticks(x_positions);
xticklabels(condition_labels);
ylabel('Average Reaction Time (seconds)');
title('Median Fixation Duration');
%legend([bar1, bar2], {'ADHD', 'nonADHD'}, 'Location', 'northeast', 'Box', 'off');
legend('Location', 'northeastoutside', 'Orientation', 'vertical', 'NumColumns', 1); % Two-column legend
hold off;
if safe
    %set(gcf, 'Position', [50, 50, 1200, 600]);
    saveas(gcf, fullfile(comparison_results_folder, safe_name));
end

end