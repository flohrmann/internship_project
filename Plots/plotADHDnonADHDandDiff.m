function plotADHDnonADHDandDiff(bigtitle,... % sgtitle
                        data_1, group_1_medians, group_1_sems, label_1, legend_1_place, ...  % adhd data
                        data_2, group_2_medians, group_2_sems, label_2, legend_2_place, ...  % nonadhd data
                        ids,...     % all participant ids
                        x_labels, y_label, ...% x, y axis labels
                        group_labels, conditions, color_map, participant_map, ...
                        safe_name)
                    
% data_1&2 : [participants x conditions]
ymax = max(max(data_1(:)), max(data_2(:)))+ max(max(group_1_sems(:))*1.2, max(group_2_sems(:))); % Maximum individual median value
ymin = min(min(data_1(:)), min(data_2(:)))- max(max(group_1_sems(:)), max(group_2_sems(:)))/2; % Maximum individual median value
sharedYLimits = [ymin, ymax];

%
figure; sgtitle(bigtitle);
% ADHD Plot
subplot(1, 3, 1);
plotStandardTwoGroupsEachData(data_1, group_1_medians, group_1_sems,label_1, y_label, legend_1_place, ids, group_labels, conditions, x_labels, color_map, participant_map, sharedYLimits);

% nonADHD Plot
subplot(1, 3, 2);
plotStandardTwoGroupsEachData(data_2, group_2_medians, group_2_sems, label_2, y_label, legend_2_place, ids, group_labels, conditions, x_labels, color_map, participant_map, sharedYLimits);

% Difference Plot
subplot(1, 3, 3); 
plotDifferenceAndTheirError(data_1, group_1_medians, group_1_sems, data_2, group_2_medians, group_2_sems, conditions, x_labels, color_map, y_label);
  
set(gcf, 'Position', [50, 50, 1400, 700]); % Resize the figure window (x, y, width, height)
saveas(gcf, safe_name);

end



function plotStandardTwoGroupsEachData(data, medians, sems, group, y_label, position, ids, group_labels, conditions, labels, color_map, participant_map, sharedYLimits)
    % Plot group medians and scatter individual participant data
    numConditions = length(medians);
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
       errorbar(c, medians(c), sems(c), '.', 'LineWidth', 1.5, 'CapSize', 10, 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    end
    xticks(1:numConditions); xticklabels(labels);
    ylabel(y_label);  ylim(sharedYLimits);
    legend('boxoff');
    legend('Location', position, 'Orientation', 'horizontal', 'NumColumns', 2); % Two-column legend
    title(group);    hold off;
end

%%
function plotDifferenceAndTheirError(adhd_data, adhd_avg, adhd_sem, nonadhd_data, nonadhd_avg, nonAdhd_sem, conditions, labels, color_map, y_label)
    % Calculate diff and perform T-tests todo permutation tests
    difference = adhd_avg - nonadhd_avg;
    % difference in standard error of the mean 
    difference_sem = sqrt(adhd_sem.^2 + nonAdhd_sem.^2) ;

    % Perform statistical tests and store p-values
    numConditions = size(adhd_data, 2);
    p_values = nan(1, numConditions); % Store p-values for annotation
    legendEntries = cell(1, numConditions); % Store legend labels
    
    hold on;
    for c = 1:numConditions
        % Extract medians for each condition
        c_adhd_data = adhd_data(:, c); % ADHD group medians for condition c
        c_nonadhd_data = nonadhd_data(:, c); % nonADHD group medians for condition c
        
        % Perform a T-test (or permutation test if needed)
        %p_values(c) = p; % Store p-value
        

        % Perform Wilcoxon Rank-Sum Test (Mann-Whitney U test)
        [p, h, stats] = ranksum(c_adhd_data, c_nonadhd_data);
        if h     
            % permutation test
            p_perm = permutationTestWelch(c_adhd_data, c_nonadhd_data);
            [~, p_tt] = ttest2(c_adhd_data, c_nonadhd_data);
            fprintf('Wilcoxon Rank-Sum Test Results:\n');
            fprintf('p-value: %.4f\n', p);
            fprintf(['h0 rejected (groups means diff at 5perc sign):', log2str(h),'\n']);
            fprintf('permutated p %.2f\n', p_perm);
            fprintf('ttest p %.2f\n', p_tt);
        else
        end
        bar(c, difference(c), 'FaceColor', color_map(conditions{c}), 'EdgeColor', 'none', ...
            'FaceAlpha', 0.3, 'DisplayName', sprintf('%s Variance', conditions{c}));
        
        errorbar(c, difference(c), difference_sem(c), '.', 'LineWidth', 1.5, 'CapSize', 10, 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
        legendEntries{c} = sprintf('p=%.2f', p);
    end
    title('ADHD vs. nonADHD');
    xticks(1:numConditions); xticklabels(labels);
    ylabel(strcat('Difference in',{' '}, y_label))
    
    lgd = legend(legendEntries, 'Location', 'southoutside');
    legend('boxoff'); %lgd.ItemTextAlignment = 'left';
    title(lgd, 'ADHD/nonADHD: Wilcoxon Rank-Sum Test');

    hold off;
end

function str=log2str(a)
if a
    str='true';
else
    str='false';
end
end
