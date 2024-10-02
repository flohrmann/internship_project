function plotFixationSpatialDistribution(fixationStats, group_labels, conditions, comparison_results_folder, safe)
    % Initialize containers for ADHD and nonADHD spatial data
    spatialX_adhd_conditions = [];
    spatialY_adhd_conditions = [];
    spatialX_nonadhd_conditions = [];
    spatialY_nonadhd_conditions = [];
    
    condition_adhd = {};
    condition_nonadhd = {};

    % Loop through the struct to gather spatial data for ADHD and non-ADHD groups
    for i = 1:length(fixationStats)
        for c = 1:length(conditions)
            conditionIdx = strcmp({fixationStats(i).conditions.name}, conditions{c});
            
            if any(conditionIdx)
                conditionData = fixationStats(i).conditions(conditionIdx);
                
                if strcmp(group_labels{i}, 'ADHD')
                    % Append the fixation spatial data for each condition
                    spatialX_adhd_conditions = [spatialX_adhd_conditions; conditionData.spatialDistributionX(:)];
                    spatialY_adhd_conditions = [spatialY_adhd_conditions; conditionData.spatialDistributionY(:)];
                    condition_adhd = [condition_adhd; repmat({conditions{c}}, length(conditionData.spatialDistributionX), 1)];
                elseif strcmp(group_labels{i}, 'nonADHD')
                    % Append the fixation spatial data for each condition
                    spatialX_nonadhd_conditions = [spatialX_nonadhd_conditions; conditionData.spatialDistributionX(:)];
                    spatialY_nonadhd_conditions = [spatialY_nonadhd_conditions; conditionData.spatialDistributionY(:)];
                    condition_nonadhd = [condition_nonadhd; repmat({conditions{c}}, length(conditionData.spatialDistributionX), 1)];
                end
            end
        end
    end

    % Convert the condition labels to categorical and extract unique conditions
    condition_adhd = categorical(condition_adhd);
    condition_nonadhd = categorical(condition_nonadhd);
    unique_conditions = categories(condition_adhd); % Assuming both groups have the same conditions

    % Plot spatial distribution for each condition
    figure;
    for i = 1:length(unique_conditions)
        % Get data for the current condition
        X_adhd = spatialX_adhd_conditions(condition_adhd == unique_conditions{i});
        Y_adhd = spatialY_adhd_conditions(condition_adhd == unique_conditions{i});
        X_nonadhd = spatialX_nonadhd_conditions(condition_nonadhd == unique_conditions{i});
        Y_nonadhd = spatialY_nonadhd_conditions(condition_nonadhd == unique_conditions{i});
        
        % Plot ADHD group
        subplot(2, length(unique_conditions), i);
        scatter(X_adhd, Y_adhd, 'r.');
        title(['ADHD - ', char(unique_conditions{i})]);
        xlabel('X Position (pixels)');
        ylabel('Y Position (pixels)');
        axis equal;
        grid on;
        
        % Plot non-ADHD group
        subplot(2, length(unique_conditions), i + length(unique_conditions));
        scatter(X_nonadhd, Y_nonadhd, 'b.');
        title(['non-ADHD - ', char(unique_conditions{i})]);
        xlabel('X Position (pixels)');
        ylabel('Y Position (pixels)');
        axis equal;
        grid on;
    end

    sgtitle('Spatial Distribution of Fixations by Condition and Group');

    if safe == 1
        set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        saveas(gcf, fullfile(comparison_results_folder, 'spatial_distribution_fixations.png'));
    end
end
