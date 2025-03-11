function participant_data = plotAbandonReturnTrials(is_abandon, data_struct, all_fixations, unique_conditions, condition_labels, group_labels, color_map, comparison_results_folder, safe)

%% check if target was fixated
for p=1:size(data_struct,2)
    fixations = all_fixations(p).fixations;
    %id = all_fixations(p).id;
    for trial=1:size(fixations.targetCenters, 2)
        target_x = fixations.targetCenters(1, trial); % Coordinates of the target TARGET(x,y)
        target_y = fixations.targetCenters(2, trial);
        current_fixations = fixations.stimulusFixations(trial).fixations;
        
        for f=1:size(current_fixations, 2)
            fix_x = current_fixations(f).center(1,1);
            fix_y = current_fixations(f).center(1,2);
            dist = sqrt((target_x - fix_x)^2 + (target_y - fix_y)^2);
            if dist <=is_abandon % pixels
                all_fixations(p).fixations.stimulusFixations(trial).fixations(1,f).target = 1;
            else
                all_fixations(p).fixations.stimulusFixations(trial).fixations(1,f).target = 0;
            end
        end
    end
end

%% Count Fixations per Trial & check if target abandoned
participant_data = struct(); % To store data per participant
for p = 1:size(all_fixations, 2)
    subj_fixation = all_fixations(p).fixations;
    id = all_fixations(p).id;
    current_conditions = data_struct(p).Condition;
    fixation_counts = [];  % To store fixation counts per trial
    abandonments = [];     % To store abandonment statuses per trial
    last_found = [];
    abandon_trials = [];
    for trial = 1:size(subj_fixation.stimulusFixations, 2)
        current_fixation = subj_fixation.stimulusFixations(trial).fixations;
        try% Extract target information for fixations
            target_fixations = [current_fixation.target];
        catch % if empty
            target_fixations = 0;
        end
        % Count total target fixations
        num_target_fixations = sum(target_fixations == 1); % change so it only counts when seperated by 0
        
        % Check abandonment:
        % 1. If a target fixation is followed by a non-target fixation
        % 2. If the target fixation is the last fixation
        abandon = 0;
        was_abandoned = 0;
        for f = 1:numel(target_fixations) - 1
            if target_fixations(f) == 1 && target_fixations(f + 1) == 0
                was_abandoned = was_abandoned + 1;
                abandon = 1;
            end
        end
        % Check if the last fixation was on the target
        if ~isempty(target_fixations) && target_fixations(end) == 1
            found = 1; % Last fixation is on target
        else
            found = 0; % Last fixation is not on target
        end
        if target_fixations(end) == 1
            found = 1; % Last fixation is on target, no abandonment
        end
        fixation_counts = [fixation_counts; num_target_fixations];
        abandonments = [abandonments; was_abandoned];
        last_found = [found, last_found];
        abandon_trials = [abandon_trials; abandon];
    end
    participant_data(p).id = id;
    participant_data(p).fixation_counts = fixation_counts;
    participant_data(p).abandonments = abandonments;
    participant_data(p).last_found = last_found;
    participant_data(p).conditions = current_conditions;
    participant_data(p).abandon_trials = abandon_trials;
end

%% Separate into ADHD and nonADHD groups
adhd_indices = strcmp(group_labels, 'ADHD');
nonadhd_indices = strcmp(group_labels, 'nonADHD');
num_conditions = numel(unique_conditions);
adhd_summary = struct();
nonadhd_summary = struct();

% ADHD group
adhd_trials = find(adhd_indices); % Indices of ADHD participants
for idx = 1:numel(adhd_trials)
    p = adhd_trials(idx); % Extract the participant index
    adhd_summary(p).id = participant_data(p).id; % Save participant ID
    adhd_summary(p).counts = zeros(num_conditions, 3); % [cond x categories]
    % Loop through trials for this participant
    for trial = 1:numel(participant_data(p).conditions)
        for c = 1:num_conditions
            if strcmpi(participant_data(p).conditions{trial}, unique_conditions{c})
                % Check conditions and increment counts
                if participant_data(p).abandonments(trial) == 0 && participant_data(p).last_found(trial) == 1
                    adhd_summary(p).counts(c, 1) = adhd_summary(p).counts(c, 1) + 1;
                elseif participant_data(p).abandonments(trial) == 1 && participant_data(p).last_found(trial) == 1
                    adhd_summary(p).counts(c, 2) = adhd_summary(p).counts(c, 2) + 1;
                elseif participant_data(p).abandonments(trial) == 0 && participant_data(p).last_found(trial) == 0
                    adhd_summary(p).counts(c, 3) = adhd_summary(p).counts(c, 3) + 1;
                end
            end
        end
    end
end

% nonADHD group
nonadhd_trials = find(nonadhd_indices); % Indices of nonADHD participants
for idx = 1:numel(nonadhd_trials)
    p = nonadhd_trials(idx); % Extract the participant index
    nonadhd_summary(p).id = participant_data(p).id; % Save participant ID
    nonadhd_summary(p).counts = zeros(num_conditions, 3); % [cond x categories]
    for trial = 1:numel(participant_data(p).conditions)
        for c = 1:num_conditions
            if strcmpi(participant_data(p).conditions{trial}, unique_conditions{c})
                % Check conditions and increment counts
                % found and didnt abandon
                if participant_data(p).abandonments(trial) == 0 && participant_data(p).last_found(trial) == 1
                    nonadhd_summary(p).counts(c, 1) = nonadhd_summary(p).counts(c, 1) + 1;
                    % found, abandoned and found
                elseif participant_data(p).abandonments(trial) == 1 && participant_data(p).last_found(trial) == 1
                    nonadhd_summary(p).counts(c, 2) = nonadhd_summary(p).counts(c, 2) + 1;
                    % found, abandoned and didnt find again
                elseif participant_data(p).abandonments(trial) == 0 && participant_data(p).last_found(trial) == 0
                    nonadhd_summary(p).counts(c, 3) = nonadhd_summary(p).counts(c, 3) + 1;
                end
            end
        end
    end
end
adhd_summary_struct = adhd_summary(~cellfun('isempty', {adhd_summary.id}));
nonadhd_summary_struct = nonadhd_summary(~cellfun('isempty', {nonadhd_summary.id}));



%% barplot per participant all in one
adhd_ids = {adhd_summary_struct.id}; % Get participant IDs
nonadhd_ids = {nonadhd_summary_struct.id};
num_adhd = numel(adhd_summary_struct);
num_nonadhd = numel(nonadhd_summary_struct);

adhd_means = zeros(num_adhd, 3); % Preallocate for ADHD participants
nonadhd_means = zeros(num_nonadhd, 3); % Preallocate for Non-ADHD participants

for i = 1:num_adhd% Compute means for ADHD participants
    counts = adhd_summary_struct(i).counts; % Get counts for this participant
    adhd_means(i, :) = mean(counts, 1); % Compute mean across conditions
end

for i = 1:num_nonadhd% Compute means for Non-ADHD participants
    counts = nonadhd_summary_struct(i).counts; % Get counts for this participant
    nonadhd_means(i, :) = mean(counts, 1); % Compute mean across conditions
end

y = [adhd_means; nonadhd_means]; % Combine means into one matrix
all_ids = [adhd_ids, nonadhd_ids]; % Combine participant IDs

figure;
bar(y, 'stacked');
xticks(1:size(y, 1));
xticklabels(all_ids);%xtickangle(45);
xlabel('Participants');
ylabel('Mean Number of Trials Per Condition (50 Trials)');
title('Mean Search Strategies by Participant');
legend({'Abandonment 0, Last Found 1', ...
    'Abandonment 1, Last Found 1', ...
    'Abandonment 0, Last Found 0'}, ...
    'Location', 'northeastoutside');
if safe == 1
    set(gcf, 'Position', [50, 50, 1200, 600]);
    saveas(gcf, fullfile(comparison_results_folder, '09_abandon_return_trials_abandonment_kinds_all.png'));
end


% cant be done on matlab 2018/standrechner
try
    %% figure per group subplot per participant
    legend_items_strategies = {'Found & Stayed', 'Abandon & find again', 'No fixation'};
    
    figure;% Plot for ADHD participants
    % FIX: Add extra row for legend if number of participants isnt even
    %tiledlayout(ceil(numel(adhd_summary_struct) / 2) + 1, 2, 'TileSpacing', 'compact');
    tiledlayout(ceil(numel(adhd_summary_struct) / 3), 3, 'TileSpacing', 'compact');
    for i = 1:numel(adhd_summary_struct)
        % Extract data for this participant
        counts = adhd_summary_struct(i).counts; % [conditions x strategies]
        participant_id = adhd_summary_struct(i).id;
        
        nexttile; b = bar(counts, 'stacked'); % Create subplot
        xticks(1:num_conditions); xticklabels(condition_labels);
        ylabel('Trials'); ylim([0 50]);
        title(['ID: ', num2str(participant_id)]);
        if i == numel(adhd_summary_struct)
            ax = nexttile(1);
            lg  = legend(legend_items_strategies, 'Orientation','Vertical','Box','off', 'Location', 'northwest');
            lg.Layout.Tile = i+1;
            sgtitle('Target Fixation Abandonments per Trial (ADHD)');
        end
    end
    if safe == 1
        set(gcf, 'Position', [50, 50, 1000, 1000]);
        saveas(gcf, fullfile(comparison_results_folder, '09_abandon_return_trials_abandonment_kinds_adhd.png'));
    end
    
    figure;% Plot for Non-ADHD participants
    tiledlayout(ceil(numel(nonadhd_summary_struct) / 3), 3, 'TileSpacing', 'compact');
    for i = 1:numel(nonadhd_summary_struct)
        counts = nonadhd_summary_struct(i).counts; % [conditions x strategies]
        participant_id = nonadhd_summary_struct(i).id;
        
        nexttile; b = bar(counts, 'stacked'); % Create subplot
        xticks(1:num_conditions); xticklabels(condition_labels);
        ylabel('Trials'); ylim([0 50]);
        title(['ID: ', num2str(participant_id)]);
        if i == numel(nonadhd_summary_struct)
            ax = nexttile(1);
            lg  = legend(legend_items_strategies, 'Orientation','Vertical','Box','off', 'Location', 'northwest');
            lg.Layout.Tile = i+1;
            sgtitle('Target Fixation Abandonments per Trial (nonADHD)');
        end
    end
    if safe == 1
        set(gcf, 'Position', [50, 50, 1000, 1000]);
        saveas(gcf, fullfile(comparison_results_folder, '09_abandon_return_trials_abandonment_kinds_nonadhd.png'));
    end
    
    
    %% average per adhd/nonadhd (not shown: multiple abandonments)
    nonadhd_counts = zeros(num_conditions, 3, numel(nonadhd_summary_struct)); % [conditions x strategies x participants]
    adhd_counts = zeros(num_conditions, 3, numel(adhd_summary_struct));

    for i = 1:numel(adhd_summary_struct) % ADHD group
        for c = 1:num_conditions % data for this participant & condition
            adhd_counts(c, :, i) = adhd_summary_struct(i).counts(c, :);
        end
    end
    avg_adhd = median(adhd_counts, 3);

    for i = 1:numel(nonadhd_summary_struct) % ADHD group
        for c = 1:num_conditions % data for this participant & condition
            nonadhd_counts(c, :, i) = nonadhd_summary_struct(i).counts(c, :);
        end
    end
    avg_nonadhd = median(nonadhd_counts, 3);

 
    % Combine averages for plotting and add spacing
    spacer = nan(1, 3); % Create a spacer row
    group_averages = [];
    for c = 1:num_conditions
        group_averages = [group_averages; avg_adhd(c, :); avg_nonadhd(c, :); spacer];
    end
    group_averages(end,:) = []; % remove last spacer to safe space in plot
    
    figure; bar(group_averages, 'stacked');
    num_groups = num_conditions;
    xticks(1.5:3:(num_groups * 3));
    xticklabels(repmat({''}, 1, num_conditions)); % Temporarily remove group labels); % xlabel('Condition');
    ylabel('Average Number of Trials');
    title('Average Search Strategies per Condition by Group');
    legend(legend_items_strategies, 'Location', 'northeast', 'Box','off');
    for i = 1:num_conditions% Add group labels below x-axis % Add condition labels slightly below ADHD/Non-ADHD
        text(1 + (i - 1) * 3, -1, 'ADHD', 'HorizontalAlignment', 'center'); % ADHD label
        text(2 + (i - 1) * 3, -1, 'nonADHD', 'HorizontalAlignment', 'center'); % Non-ADHD label
        text(1.5 + (i - 1) * 3, -2.6, condition_labels{i}, 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold'); % Condition label
    end
    hold off;
    if safe == 1
        %set(gcf, 'Position', [50, 50, 1200, 600]);
        set(gcf, 'Position', [50, 50, 900, 600]);
        saveas(gcf, fullfile(comparison_results_folder, '09_abandon_return_trials_num_abandonment_kinds_group.png'));
    end
    
    
    %% bin number of abandonments per participant and condition
    
    max_abandonments = 10;
    max_abandonments = 4;
    
    legend_items_count = {'1x Abandonment', '2x Abandonments', '3x Abandonments', '4x Abandonments', '5x Abandonments', ...
        '6x Abandonments', '7x Abandonments', '8x Abandonments', '9x Abandonments', '>= 10x Abandonments'};
    legend_items_count = {'1x Abandonment', '2x Abandonments', '3x Abandonments', '>= 4x Abandonments'};
    
    colormap_matrix = [linspace(0, 1, max_abandonments)', zeros(max_abandonments, 1), linspace(1, 0, max_abandonments)'];
    
    figure;% ADHD Participants Plot
    tiledlayout(ceil(numel(adhd_summary_struct) / 3), 3, 'TileSpacing', 'compact');
    %tiledlayout(ceil(numel(adhd_summary_struct) / 2) + 1, 2, 'TileSpacing', 'compact'); % Add extra row for legend
    for i = 1:numel(adhd_summary_struct)
        participant_id = adhd_summary_struct(i).id;
        counts = zeros(num_conditions, max_abandonments); % Preallocate [conditions x abandonments]
        for c = 1:num_conditions    % Loop through conditions
            cond_trials = strcmpi(participant_data(adhd_trials(i)).conditions, unique_conditions{c}); % Trials for this condition
            for x = 1:max_abandonments
                counts(c, x) = sum(sum(participant_data(adhd_trials(i)).abandonments(cond_trials) == x));
            end
        end
        nexttile; b = bar(counts, 'stacked');
        for k = 1:max_abandonments
            b(k).FaceColor = 'flat';
            b(k).CData = repmat(colormap_matrix(k, :), num_conditions, 1); % Apply gradient color
        end
        xticks(1:num_conditions); xticklabels(condition_labels);
        ylabel('Trials'); %xlabel('Condition');
        ylim([0 50]);
        title(['ID: ', num2str(participant_id)]);
        if i == numel(adhd_summary_struct)
            ax = nexttile(1);
            lg  = legend(legend_items_count, 'Orientation','Vertical','NumColumns',2, 'Box','off');
            lg.Layout.Tile = i+1;
            sgtitle('Target Fixation Abandonments per Trial (ADHD)');
        end
    end
    if safe == 1
        set(gcf, 'Position', [50, 50, 1000, 600]); % Resize the figure window (x, y, width, height)
        saveas(gcf, fullfile(comparison_results_folder, strcat('09_abandon_return_trials_',num2str(max_abandonments),'_abandonments_adhd_individual.png')));
    end
    
    figure;% Non-ADHD Participants Plot
    tiledlayout(ceil(numel(nonadhd_summary_struct) / 3), 3, 'TileSpacing', 'compact');
    %tiledlayout(ceil(numel(nonadhd_summary_struct) / 2) + 1, 2, 'TileSpacing', 'compact'); % Add extra row for legend
    for i = 1:numel(nonadhd_summary_struct)
        participant_id = nonadhd_summary_struct(i).id;
        counts = zeros(num_conditions, max_abandonments); % Preallocate [conditions x abandonments]
        for c = 1:num_conditions
            cond_trials = strcmpi(participant_data(nonadhd_trials(i)).conditions, unique_conditions{c}); % Trials for this condition
            for x = 1:max_abandonments
                counts(c, x) = sum(sum(participant_data(nonadhd_trials(i)).abandonments(cond_trials) == x));
            end
        end
        nexttile; b2 = bar(counts, 'stacked');
        for k = 1:max_abandonments
            b2(k).FaceColor = 'flat';
            b2(k).CData = repmat(colormap_matrix(k, :), num_conditions, 1);
        end
        xticks(1:num_conditions); xticklabels(condition_labels);
        ylabel('Trials'); % xlabel('Condition');
        ylim([0 50]);
        title(['ID: ', num2str(participant_id)]);
        if i == numel(nonadhd_summary_struct)
            ax = nexttile(1);
            lg  = legend(legend_items_count, 'Orientation','Vertical','NumColumns',2, 'Box','off');
            lg.Layout.Tile = i+1;
            sgtitle('Target Fixation Abandonments per Trial (nonADHD)');
        end
    end
    if safe == 1
        set(gcf, 'Position', [50, 50, 1000, 600]); % set(gcf, 'Position', [50, 50, 1200, 600]);
        saveas(gcf, fullfile(comparison_results_folder, strcat('09_abandon_return_trials_',num2str(max_abandonments),'_abandonments_nonadhd_individual.png')));
    end
    
    %% amount of abandonments of target mean per group
    group_adhd_means = nan(num_conditions, max_abandonments); % ADHD means [conditions x abandonments]
    group_nonadhd_means = nan(num_conditions, max_abandonments); % Non-ADHD means [conditions x abandonments]
    % means for ADHD group
    adhd_participant_counts = nan(numel(adhd_summary_struct), num_conditions, max_abandonments); % Storage for all participants
    for i = 1:numel(adhd_summary_struct)
        for c = 1:num_conditions
            cond_trials = strcmpi(participant_data(adhd_trials(i)).conditions, unique_conditions{c}); % Trials for this condition
            for x = 1:max_abandonments
                adhd_participant_counts(i, c, x) = sum(participant_data(adhd_trials(i)).abandonments(cond_trials) == x);
            end
        end
    end
    group_adhd_means = squeeze(nanmean(adhd_participant_counts, 1)); % Average across participants
    % means for Non-ADHD group
    nonadhd_participant_counts = nan(numel(nonadhd_summary_struct), num_conditions, max_abandonments); % Storage for all participants
    for i = 1:numel(nonadhd_summary_struct)
        for c = 1:num_conditions
            cond_trials = strcmpi(participant_data(nonadhd_trials(i)).conditions, unique_conditions{c}); % Trials for this condition
            for x = 1:max_abandonments
                nonadhd_participant_counts(i, c, x) = sum(participant_data(nonadhd_trials(i)).abandonments(cond_trials) == x);
            end
        end
    end
    group_nonadhd_means = squeeze(nanmean(nonadhd_participant_counts, 1)); % Average across participants
    
    spacer = nan(1, max_abandonments); combined_means = [];
    for c = 1:num_conditions
        combined_means = [combined_means; group_adhd_means(c, :); group_nonadhd_means(c, :); spacer];
    end
    combined_means(end,:) = []; % remove last spacer to safe space in plot
    
    figure;tiledlayout(1, 1, 'TileSpacing', 'compact'); nexttile;
    b = bar(combined_means, 'stacked');
    for k = 1:max_abandonments
        b(k).FaceColor = 'flat';
        b(k).CData = repmat(colormap_matrix(k, :), size(combined_means, 1), 1); % Apply gradient color
    end
    xticks(1.5:3:(2 * num_conditions * 3));
    xticklabels(repmat({''}, 1, num_conditions*2)); % Temporarily remove group labels); % xlabel('Condition');
    ylabel('Average Number of Trials');
    title('Average Target Fixation Abandonments per Condition by Group');
    legend(legend_items_count, 'Location', 'northeast', 'Box', 'off');
    for i = 1:num_conditions% Add group labels below x-axis % Add condition labels slightly below ADHD/Non-ADHD
        text(1 + (i - 1) * 3, -1, 'ADHD', 'HorizontalAlignment', 'center'); % ADHD label
        text(2 + (i - 1) * 3, -1, 'nonADHD', 'HorizontalAlignment', 'center'); % Non-ADHD label
        text(1.5 + (i - 1) * 3, -2.6, condition_labels{i}, 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold'); % Condition label
    end
    if safe == 1
        set(gcf, 'Position', [50, 50, 1000, 600]);
        saveas(gcf, fullfile(comparison_results_folder, strcat('09_abandon_return_trials_',num2str(max_abandonments),'_abandonments_group_means.png')));
    end
    
    
catch
end

