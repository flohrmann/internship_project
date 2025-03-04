function plotAARtrialsRThandeye(participant_data, data_struct, ids, unique_conditions, condition_labels, group_labels, color_map, color_map_individual, comparison_results_folder, safe)


rt_struct = struct();
for p = 1:size(participant_data, 2)
    subj_aarts = participant_data(p).abandon_trials;
    subj_rt_eye = data_struct(p).rt_eye;
    subj_rt_hand = data_struct(p).rt;
    id = participant_data(p).id;
    current_conditions = participant_data(p).conditions;
    aart_eye = [];
    nonaart_eye = [];
    aart_hand = [];
    nonaart_hand = [];
    aart_diff = [];
    nonaart_diff = [];
    
    for trial = 1:size(subj_aarts, 1)
        if subj_aarts(trial) == 1
            try
                aart_eye = [aart_eye; subj_rt_eye(trial)];
                aart_hand = [aart_hand; subj_rt_hand(trial)];
                aart_diff = [aart_diff; (subj_rt_hand(trial) - subj_rt_eye(trial))];
                nonaart_eye = [nonaart_eye; NaN];
                nonaart_hand = [nonaart_hand; NaN];
                nonaart_diff = [nonaart_diff; NaN];
                
            catch % no valid rt for this trial
                aart_eye = [aart_eye; NaN];
                aart_hand = [aart_hand; NaN];
                aart_diff = [aart_diff; NaN];
                nonaart_eye = [nonaart_eye; NaN];
                nonaart_hand = [nonaart_hand; NaN];
                nonaart_diff = [nonaart_diff; NaN];
            end
        else % not a aar trial
            try
                aart_eye = [aart_eye; NaN];
                aart_hand = [aart_hand; NaN];
                aart_diff = [aart_diff; NaN];
                nonaart_eye = [nonaart_eye; subj_rt_eye(trial)];
                nonaart_hand = [nonaart_hand; subj_rt_hand(trial)];
                nonaart_diff = [nonaart_diff; (subj_rt_hand(trial) - subj_rt_eye(trial))];
            catch % no valid rt for this trial
                aart_eye = [aart_eye; NaN];
                aart_hand = [aart_hand; NaN];
                aart_diff = [aart_diff; NaN];
                nonaart_eye = [nonaart_eye; NaN];
                nonaart_hand = [nonaart_hand; NaN];
                nonaart_diff = [nonaart_diff; NaN];
            end
        end
    end
    rt_struct(p).id = id;
    rt_struct(p).aart_eye = aart_eye;
    rt_struct(p).nonaart_eye = nonaart_eye;
    rt_struct(p).aart_hand = aart_hand;
    rt_struct(p).nonaart_hand = nonaart_hand;
    rt_struct(p).aart_diff = aart_diff;
    rt_struct(p).nonaart_diff = nonaart_diff;
    rt_struct(p).conditions = current_conditions;
end

%% for all 4 conditions
nc = [2, numel(unique_conditions)];

for which_condition = 1:2
    num_conditions = nc(which_condition);
    
    nonadhd_indices = find(strcmp(group_labels, 'nonADHD'));
    num_nonadhd = numel(nonadhd_indices);
    adhd_indices = find(strcmp(group_labels, 'ADHD'));
    num_adhd = numel(adhd_indices);
    
    adhd_summary = struct();
    nonadhd_summary = struct();
    
    adhd_trials = find(adhd_indices); % Indices of ADHD participants
    nonadhd_trials = find(nonadhd_indices); % Indices of nonADHD participants
    
    
    % Initialize storage for participant means
    participant_avgs = struct();
    for p = 1:numel(rt_struct)
        participant_id = rt_struct(p).id;
        rt_hand        = rt_struct(p).aart_hand;
        n_rt_hand      = rt_struct(p).nonaart_hand;
        rt_eye         = rt_struct(p).aart_eye;
        n_rt_eye       = rt_struct(p).nonaart_eye;
        conditions     = rt_struct(p).conditions;
        aart_diff      = rt_struct(p).aart_diff;
        nonaart_diff    = rt_struct(p).nonaart_diff;
        avg_hand   = nan(num_conditions, 1); avg_eye = nan(num_conditions, 1);
        n_avg_hand = nan(num_conditions, 1); n_avg_eye = nan(num_conditions, 1);
        avg_diff = nan(num_conditions, 1); n_avg_diff = nan(num_conditions, 1);
        
        
        if num_conditions == 2
            for c = 1:num_conditions
                % Find trials for a/b only    
                %condition_trials_1 = strcmp(conditions, unique_conditions{c});
                %condition_trials_2 = strcmp(conditions, unique_conditions{c+1});
                %condition_trials  = condition_trials_1 +condition_trials_2;
                condition_trials= strcmp(conditions, unique_conditions{c});
                avg_hand(c)     = nanmedian(rt_hand(condition_trials));
                avg_eye(c)      = nanmedian(rt_eye(condition_trials));
                n_avg_hand(c)   = nanmedian(n_rt_hand(condition_trials));
                n_avg_eye(c)    = nanmedian(n_rt_eye(condition_trials));
                avg_diff(c)     = nanmedian(aart_diff(condition_trials));
                n_avg_diff(c)   = nanmedian(nonaart_diff(condition_trials));
            end
        else
            for c = 1:num_conditions
                % Find trials for this condition
                condition_trials= strcmp(conditions, unique_conditions{c});
                avg_hand(c)     = nanmedian(rt_hand(condition_trials));
                avg_eye(c)      = nanmedian(rt_eye(condition_trials));
                n_avg_hand(c)   = nanmedian(n_rt_hand(condition_trials));
                n_avg_eye(c)    = nanmedian(n_rt_eye(condition_trials));
                avg_diff(c)     = nanmedian(aart_diff(condition_trials));
                n_avg_diff(c)   = nanmedian(nonaart_diff(condition_trials));
            end
            
        end
        

        participant_avgs(p).id         = participant_id;
        participant_avgs(p).avg_hand   = avg_hand;
        participant_avgs(p).avg_eye    = avg_eye;
        participant_avgs(p).n_avg_hand = n_avg_hand;
        participant_avgs(p).n_avg_eye  = n_avg_eye;
        participant_avgs(p).avg_diff   = avg_diff;
        participant_avgs(p).n_avg_diff = n_avg_diff;
    end
    
    % Group data for ADHD and nonADHD
    adhd_avg_hand = [];     adhd_avg_eye = [];      adhd_avg_diff = [];
    nonadhd_avg_hand = [];  nonadhd_avg_eye = [];   nonadhd_avg_diff   = [];
    adhd_avg_n_hand = [];   adhd_avg_n_eye = [];    adhd_avg_n_diff = [];
    nonadhd_avg_n_hand = [];nonadhd_avg_n_eye = []; nonadhd_avg_n_diff = [];
    
    for p = 1:numel(participant_avgs)
        if ~(sum(ismember(adhd_indices, p))==0) % if its in the list of adhd subjects
            adhd_avg_hand   = [adhd_avg_hand; participant_avgs(p).avg_hand'];
            adhd_avg_eye    = [adhd_avg_eye; participant_avgs(p).avg_eye'];
            adhd_avg_diff   = [adhd_avg_diff; participant_avgs(p).avg_diff'];
            adhd_avg_n_hand = [adhd_avg_n_hand; participant_avgs(p).n_avg_hand'];
            adhd_avg_n_eye  = [adhd_avg_n_eye; participant_avgs(p).n_avg_eye'];
            adhd_avg_n_diff = [adhd_avg_n_diff; participant_avgs(p).n_avg_diff'];
        elseif ~(sum(ismember(nonadhd_indices, p))==0) % if its in the list of adhd subjects
            nonadhd_avg_hand = [nonadhd_avg_hand; participant_avgs(p).avg_hand'];
            nonadhd_avg_eye = [nonadhd_avg_eye; participant_avgs(p).avg_eye'];
            nonadhd_avg_diff   = [nonadhd_avg_diff; participant_avgs(p).avg_diff'];
            nonadhd_avg_n_hand = [nonadhd_avg_n_hand; participant_avgs(p).n_avg_hand'];
            nonadhd_avg_n_eye = [nonadhd_avg_n_eye; participant_avgs(p).n_avg_eye'];
            nonadhd_avg_n_diff = [nonadhd_avg_n_diff; participant_avgs(p).n_avg_diff'];
        end
    end
    
    
    %% AAR percentages per participant
    aar_percentages_adhd        = zeros(numel(adhd_indices), num_conditions);
    nonaar_percentages_adhd     = zeros(numel(adhd_indices), num_conditions);
    aar_percentages_nonadhd     = zeros(numel(nonadhd_indices), num_conditions);
    nonaar_percentages_nonadhd  = zeros(numel(nonadhd_indices), num_conditions);
    
    for p = 1:numel(adhd_indices)
        participant_idx = adhd_indices(p);
        conditions = rt_struct(participant_idx).conditions;
        aars = participant_data(participant_idx).abandon_trials; % AAR flags
        for c = 1:num_conditions
            condition_trials = strcmp(conditions, unique_conditions{c});
            aar_trials       = sum(aars(condition_trials));
            total_trials     = size(condition_trials, 1) / 4;
            aar_percentages_adhd(p, c)    = (aar_trials / total_trials) * 100;
            nonaar_percentages_adhd(p, c) = ((total_trials - aar_trials) / total_trials) * 100;
        end
    end
    
    for p = 1:numel(nonadhd_indices)
        participant_idx = nonadhd_indices(p);
        conditions      = rt_struct(participant_idx).conditions;
        aars            = participant_data(participant_idx).abandon_trials;
        for c = 1:num_conditions
            condition_trials = strcmp(conditions, unique_conditions{c});
            aar_trials       = sum(aars(condition_trials));
            total_trials     = size(condition_trials, 1) / 4;
            aar_percentages_nonadhd(p, c)    = (aar_trials / total_trials) * 100;
            nonaar_percentages_nonadhd(p, c) = ((total_trials - aar_trials) / total_trials) * 100;
        end
    end
    
    aar_medians_adhd    = median(aar_percentages_adhd, 1, 'omitnan') % Median per condition
    aar_sem_adhd        = std(aar_percentages_adhd, 0, 1, 'omitnan')./ sqrt(size(aar_percentages_adhd, 1))
    aar_medians_nonadhd = median(aar_percentages_nonadhd, 1, 'omitnan')
    aar_sem_nonadhd     = std(aar_percentages_nonadhd, 0, 1, 'omitnan')./ sqrt(size(aar_percentages_nonadhd, 1))
    
    
    plotADHDnonADHDandDiff('Median Abandon Return Trial Rate Comparison',...
        aar_percentages_adhd, aar_medians_adhd, aar_sem_adhd, 'ADHD', 'northeast', ...
        aar_percentages_nonadhd, aar_medians_nonadhd, aar_sem_nonadhd, 'nonADHD', 'northeast', ...
        ids, condition_labels, '% of trials', group_labels, unique_conditions, color_map, color_map_individual, ...
        fullfile(comparison_results_folder, '09_aar_trials_percentage_allinone_median.png'));
    
    %% RT hand and eye in AAR /non AAR trials per participant/condition
    % aar trials
    median_adhd_hand      = nanmedian(adhd_avg_hand, 1);
    median_adhd_eye       = nanmedian(adhd_avg_eye, 1);
    median_nonadhd_hand   = nanmedian(nonadhd_avg_hand, 1);
    median_nonadhd_eye    = nanmedian(nonadhd_avg_eye, 1);
    se_adhd_hand        = nanstd(adhd_avg_hand, 0, 1) ./ sqrt(size(adhd_avg_hand, 1));
    se_adhd_eye         = nanstd(adhd_avg_eye, 1) ./ sqrt(size(adhd_avg_eye, 1));
    se_nonadhd_hand     = nanstd(nonadhd_avg_hand, 1) ./ sqrt(size(nonadhd_avg_hand, 1));
    se_nonadhd_eye      = nanstd(nonadhd_avg_eye, 1) ./ sqrt(size(nonadhd_avg_eye, 1));
    % non aar trials
    median_adhd_n_hand    = nanmedian(adhd_avg_n_hand, 1);
    median_adhd_n_eye     = nanmedian(adhd_avg_n_eye, 1);
    median_nonadhd_n_hand = nanmedian(nonadhd_avg_n_hand, 1);
    median_nonadhd_n_eye  = nanmedian(nonadhd_avg_n_eye, 1);
    se_adhd_n_hand      = nanstd(adhd_avg_n_hand, 1) ./ sqrt(size(adhd_avg_n_hand, 1));
    se_adhd_n_eye       = nanstd(adhd_avg_n_eye, 1) ./ sqrt(size(adhd_avg_n_eye, 1));
    se_nonadhd_n_hand   = nanstd(nonadhd_avg_n_hand, 1) ./ sqrt(size(nonadhd_avg_n_hand, 1));
    se_nonadhd_n_eye    = nanstd(nonadhd_avg_n_eye, 1) ./ sqrt(size(nonadhd_avg_n_eye, 1));
    
    % hand
    plotADHDnonADHDandDiff('Median RT of Button Press in Arrive Abandon Return Trials',...
        adhd_avg_hand, median_adhd_hand, se_adhd_hand, 'ADHD', 'northwest', ...
        nonadhd_avg_hand, median_nonadhd_hand, se_nonadhd_hand, 'nonADHD', 'northeast', ...
        ids, condition_labels, 'Median RT (s)', group_labels, unique_conditions, color_map, color_map_individual, ...
        fullfile(comparison_results_folder, strcat('09_aar_trials_rt_hand_allinone_median', num2str(num_conditions),'.png')));
    
    plotADHDnonADHDandDiff('Median RT of Button Press in non-AAR-Trials',...
        adhd_avg_n_hand, median_adhd_n_hand, se_adhd_n_hand, 'ADHD', 'northeast', ...
        nonadhd_avg_n_hand, median_nonadhd_n_hand, se_nonadhd_n_hand, 'nonADHD', 'northeast', ...
        ids, condition_labels, 'Median RT (s)', group_labels, unique_conditions, color_map, color_map_individual, ...
        fullfile(comparison_results_folder, strcat('09_nonaar_trials_rt_hand_allinone_median', num2str(num_conditions),'.png')));
    
    
    
    
    % eye
    plotADHDnonADHDandDiff('Median Gaze Arrival at Target in Arrive Abandon Return Trials',...
        adhd_avg_eye, median_adhd_eye, se_adhd_eye, 'ADHD', 'northwest', ...
        nonadhd_avg_eye, median_nonadhd_eye, se_nonadhd_eye, 'nonADHD', 'northeast', ...
        ids, condition_labels, 'Median RT (s)', group_labels, unique_conditions, color_map, color_map_individual, ...
        fullfile(comparison_results_folder, strcat('09_aar_trials_rt_eye_allinone_median', num2str(num_conditions),'.png')));
    
    plotADHDnonADHDandDiff('Median Gaze Arrival at Target in non-AAR-Trials',...
        adhd_avg_n_eye, median_adhd_n_eye, se_adhd_n_eye, 'ADHD', 'northwest', ...
        nonadhd_avg_n_eye, median_nonadhd_n_eye, se_nonadhd_n_eye, 'nonADHD', 'northeast', ...
        ids, condition_labels, 'Median RT (s)', group_labels, unique_conditions, color_map, color_map_individual, ...
        fullfile(comparison_results_folder, strcat('09_nonaar_trials_rt_eye_allinone_median', num2str(num_conditions),'.png')));
    
    % time between eye rt and button rt
    median_nonadhd_n_diff  = nanmedian(nonadhd_avg_n_diff, 1);
    se_nonadhd_n_diff      = nanstd(nonadhd_avg_n_diff, 1) ./ sqrt(size(nonadhd_avg_n_diff, 1));
    median_adhd_n_diff  = nanmedian(adhd_avg_n_diff, 1);
    se_adhd_n_diff      = nanstd(adhd_avg_n_diff, 1) ./ sqrt(size(adhd_avg_n_diff, 1));
    plotADHDnonADHDandDiff('Time between Gaze arriving at Target and Button Press in non-AAR-Trials',...
        adhd_avg_n_diff, median_adhd_n_diff, se_adhd_n_diff, 'ADHD', 'northeast', ...
        nonadhd_avg_n_diff, median_nonadhd_n_diff, se_nonadhd_n_diff, 'nonADHD', 'northeast', ...
        ids, condition_labels, 'Difference in median RT (s)', group_labels, unique_conditions, color_map, color_map_individual, ...
        fullfile(comparison_results_folder, strcat('09_nonaar_trials_rt_hand__eye_diff_allinone_median', num2str(num_conditions),'.png')));
    
    median_nonadhd_diff  = nanmedian(nonadhd_avg_diff, 1);
    se_nonadhd_diff      = nanstd(nonadhd_avg_diff, 1) ./ sqrt(size(nonadhd_avg_diff, 1));
    median_adhd_diff  = nanmedian(adhd_avg_diff, 1);
    se_adhd_diff      = nanstd(adhd_avg_diff, 1) ./ sqrt(size(adhd_avg_diff, 1));
    plotADHDnonADHDandDiff('Time between Gaze arriving at Target and Button Press in AAR-Trials',...
        adhd_avg_diff, median_adhd_diff, se_adhd_diff, 'ADHD', 'northeast', ...
        nonadhd_avg_diff, median_nonadhd_diff, se_nonadhd_diff, 'nonADHD', 'northeast', ...
        ids, condition_labels, 'Difference in median RT (s)', group_labels, unique_conditions, color_map, color_map_individual, ...
        fullfile(comparison_results_folder, strcat('09_aar_trials_rt_hand_eye_diff_allinone_median', num2str(num_conditions),'.png')));
    
    
    
    %% per particpant grouped by adhd/nonadhd
    if num_conditions == 4
    try
        %adhd
        figure; % Create a new figure
        %num_participants = numel(aar_rt_struct);
        tiledlayout(ceil(num_adhd / 3), 3, 'TileSpacing', 'compact');% Number of tiles based on the number of participants
        % Loop through each participant to create individual plots
        for i = 1:num_adhd    % Extract data for this participant
            p = adhd_indices(i);
            participant_id = rt_struct(p).id;
            rt_hand = rt_struct(p).aart_hand;
            rt_eye = rt_struct(p).aart_eye;
            conditions = rt_struct(p).conditions;
            valid_indices = ~isnan(rt_hand) & ~isnan(rt_eye);    % Remove NaN values
            rt_hand = rt_hand(valid_indices);
            rt_eye = rt_eye(valid_indices);
            conditions = conditions(valid_indices);
            % Convert conditions to categorical for consistent grouping
            condition_cats = categorical(conditions, unique_conditions, unique_conditions);
            nexttile;     % Plot data for this participant
            hold on;
            scatter(condition_cats, rt_hand, '.', 'DisplayName', 'RT Hand', 'MarkerEdgeColor', 'k'); % Plot RT_hand
            scatter(condition_cats, rt_eye, 'x', 'DisplayName', 'RT Eye', 'MarkerEdgeColor', [1, 0, 1]); % Plot RT_eye
            hold off;
            xticklabels(condition_labels);    %xticks(1:4);
            ylabel('RT (ms)');  title(['ID: ', num2str(participant_id)]);    %ylim([0, max([rt_hand; rt_eye], [], 'omitnan')]); % Set Y-limit
            if i == num_adhd    % Add legend for the last tile
                legend('RT Hand', 'RT Eye', 'Location', 'best');
            end
        end
        sgtitle('Reaction Times (RT Hand and RT Eye) per Participant [ADHD]');
        if safe
            %set(gcf, 'Position', [50, 50, 1200, 600]);
            saveas(gcf, fullfile(comparison_results_folder, '09_aar_trials_rt_eye_rt_hand_overlay_bars_conditions_individual_adhd.png'));
        end
        % nonadhd
        figure; % Create a new figure
        tiledlayout(ceil(num_nonadhd / 3), 3, 'TileSpacing', 'compact');% Number of tiles based on the number of participants
        % Loop through each participant to create individual plots
        for i = 1:num_nonadhd    % Extract data for this participant
            p = nonadhd_indices(i); % Participant index
            participant_id = rt_struct(p).id;
            rt_hand = rt_struct(p).aart_hand;
            rt_eye = rt_struct(p).aart_eye;
            conditions = rt_struct(p).conditions;
            valid_indices = ~isnan(rt_hand) & ~isnan(rt_eye);    % Remove NaN values
            rt_hand = rt_hand(valid_indices);
            rt_eye = rt_eye(valid_indices);
            conditions = conditions(valid_indices);
            condition_cats = categorical(conditions, unique_conditions, unique_conditions);        % Convert conditions to categorical for consistent grouping
            nexttile; hold on;
            scatter(condition_cats, rt_hand, '.', 'DisplayName', 'RT Hand', 'MarkerEdgeColor', 'k'); % Plot RT_hand
            scatter(condition_cats, rt_eye, 'x', 'DisplayName', 'RT Eye', 'MarkerEdgeColor', [1, 0, 1]); % Plot RT_eye
            hold off;
            xticklabels(condition_labels);    %xticks(1:4);
            ylabel('RT (ms)');  title(['ID: ', num2str(participant_id)]);    %ylim([0, max([rt_hand; rt_eye], [], 'omitnan')]); % Set Y-limit
            if i == num_nonadhd    % Add legend for the last tile
                legend('RT Hand', 'RT Eye', 'Location', 'best');
            end
        end
        sgtitle('Reaction Times (RT Hand and RT Eye) per Participant [nonADHD]');
        if safe
            %set(gcf, 'Position', [50, 50, 1200, 600]);
            saveas(gcf, fullfile(comparison_results_folder, '09_aar_trials_rt_eye_rt_hand_overlay_bars_conditions_individual_nonadhd.png'));
        end
        
    catch
    end
    
    %%
    plotGroups(median_adhd_hand, median_adhd_eye, ...
        se_adhd_hand, se_adhd_eye,...
        median_nonadhd_hand, median_nonadhd_eye, ...
        se_nonadhd_hand, se_nonadhd_eye, ...
        adhd_avg_hand, adhd_avg_eye, ...
        nonadhd_avg_hand, nonadhd_avg_eye, ...
        'RThand (lighter coloured bars) and RTeye (darker coloured bars)[AARtrials]', ...
        '09_aar_trials_rt_eye_rt_hand_overlay_bars_conditions_groups_individual.png',...
        num_conditions, color_map, comparison_results_folder, condition_labels, safe)
    %% same but for non AAR trials
    plotGroups(median_adhd_n_hand, median_adhd_n_eye, ...
        se_adhd_n_hand, se_adhd_n_eye,...
        median_nonadhd_n_hand, median_nonadhd_n_eye, ...
        se_nonadhd_n_hand, se_nonadhd_n_eye, ...
        adhd_avg_n_hand, adhd_avg_n_eye, ...
        nonadhd_avg_n_hand, nonadhd_avg_n_eye, ...
        'RThand (lighter coloured bars) and RTeye (darker coloured bars)[nonAARtrials]', ...
        '09_aar_NON_trials_rt_eye_rt_hand_overlay_bars_conditions_groups_individual.png',...
        num_conditions, color_map, comparison_results_folder, condition_labels, safe)
    
    
    
    
    
    %% % of trials that were AAR
    
    % Initialize storage for AAR percentages
    aar_percentages_adhd = zeros(num_conditions, 1);
    nonaar_percentages_adhd = zeros(num_conditions, 1);
    aar_percentages_nonadhd = zeros(num_conditions, 1);
    nonaar_percentages_nonadhd = zeros(num_conditions, 1);
    
    % Count trials for ADHD participants
    for p = 1:numel(adhd_indices)
        participant_idx = adhd_indices(p);
        conditions = rt_struct(participant_idx).conditions;
        aars = participant_data(participant_idx).abandon_trials; % AAR flags
        
        for c = 1:num_conditions
            condition_trials = strcmp(conditions, unique_conditions{c});
            aar_trials  = sum(aars(condition_trials));
            total_trials = size(condition_trials, 1)/4;
            aar_percentages_adhd(c) = aar_percentages_adhd(c) + (aar_trials / total_trials) * 100;
            nonaar_percentages_adhd(c) = nonaar_percentages_adhd(c) + ((total_trials - aar_trials) / total_trials) * 100;
            
        end
    end
    
    % Average percentages across ADHD participants
    aar_percentages_adhd = aar_percentages_adhd / numel(adhd_indices);
    nonaar_percentages_adhd = nonaar_percentages_adhd / numel(adhd_indices);
    
    % Count trials for nonADHD participants
    for p = 1:numel(nonadhd_indices)
        participant_idx = nonadhd_indices(p);
        conditions = rt_struct(participant_idx).conditions;
        aars = participant_data(participant_idx).abandon_trials; % AAR flags
        
        for c = 1:num_conditions
            condition_trials = strcmp(conditions, unique_conditions{c});
            aar_trials = sum(condition_trials & aars);
            total_trials = size(condition_trials, 1)/4;
            aar_percentages_nonadhd(c) = aar_percentages_nonadhd(c) + (aar_trials / total_trials) * 100;
            nonaar_percentages_nonadhd(c) = nonaar_percentages_nonadhd(c) +((total_trials - aar_trials) / total_trials) * 100;
        end
    end
    
    % Average percentages across nonADHD participants
    aar_percentages_nonadhd = aar_percentages_nonadhd / numel(nonadhd_indices);
    nonaar_percentages_nonadhd = nonaar_percentages_nonadhd / numel(nonadhd_indices);
    
    % Combine percentages for plotting
    percentages_adhd = [aar_percentages_adhd, nonaar_percentages_adhd];
    percentages_nonadhd = [aar_percentages_nonadhd, nonaar_percentages_nonadhd];
    
    % Plotting
    figure;
    bar_width = 0.4;
    offset = 0.2;
    
    % Bar positions
    x_positions = 1:num_conditions;
    adhd_positions = x_positions - offset;
    nonadhd_positions = x_positions + offset;
    
    hold on;
    
    % ADHD bars
    bar(adhd_positions, percentages_adhd, bar_width, 'stacked', ...
        'EdgeColor', 'none', 'DisplayName', 'ADHD');
    
    % nonADHD bars
    bar(nonadhd_positions, percentages_nonadhd, bar_width, 'stacked', ...
        'EdgeColor', 'none', 'DisplayName', 'nonADHD');
    
    % Add labels and legend
    xticks(x_positions);
    xticklabels(condition_labels);
    xlabel('Conditions');
    ylabel('% of Trials');
    title('Percentage of Trials that were AAR or non-AAR by Condition');
    legend({'AAR Trials', 'non-AAR Trials'}, 'Location', 'northeastoutside', 'Box', 'off');
    hold off;
    
    % Save figure if in safe mode
    if safe
        saveas(gcf, fullfile(comparison_results_folder, '09_aar_nonaar_percentage_conditions.png'));
    end
    
    else
    end 
end

end






function plotGroups(mean_adhd_hand, mean_adhd_eye, ...
    se_adhd_hand, se_adhd_eye,...
    mean_nonadhd_hand, mean_nonadhd_eye, ...
    se_nonadhd_hand, se_nonadhd_eye, ...
    adhd_means_hand, adhd_means_eye, ...
    nonadhd_means_hand, nonadhd_means_eye, ...
    plot_title, safe_name,...
    num_conditions, color_map, comparison_results_folder, condition_labels, safe)

bar_width = 0.4; % Bar width for both RT_hand and RT_eye
offset = 0.2;    % Offset between ADHD and nonADHD bars
x_positions = 1:num_conditions;

adhd_positions = x_positions - offset;% Adjusted positions for ADHD and nonADHD bars
nonadhd_positions = x_positions + offset;

figure; hold on;% Overlay bars for ADHD group
adhd_hand_bar = bar(adhd_positions, mean_adhd_hand, bar_width, 'FaceColor', color_map('ADHD'), ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'ADHD RT Hand');
adhd_eye_bar  = bar(adhd_positions, mean_adhd_eye, bar_width, 'FaceColor', color_map('ADHD'), ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'ADHD RT Eye');
errorbar(adhd_positions, mean_adhd_hand, se_adhd_hand, '.', 'LineWidth', 1.5, 'CapSize', 10, 'Color', '#609860');
errorbar(adhd_positions, mean_adhd_eye, se_adhd_eye, '.', 'LineWidth', 1.5, 'CapSize', 10, 'Color', '#609860');

% Overlay bars for nonADHD group
nonadhd_hand_bar = bar(nonadhd_positions, mean_nonadhd_hand, bar_width, 'FaceColor', color_map('nonADHD'), ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'nonADHD RT Hand');
nonadhd_eye_bar  = bar(nonadhd_positions, mean_nonadhd_eye, bar_width, 'FaceColor', color_map('nonADHD'), ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'nonADHD RT Eye');
errorbar(nonadhd_positions, mean_nonadhd_hand, se_nonadhd_hand, '.', 'LineWidth', 1.5, 'CapSize', 10, 'Color', '#9F3030');
errorbar(nonadhd_positions, mean_nonadhd_eye, se_nonadhd_eye, '.', 'LineWidth', 1.5, 'CapSize', 10, 'Color', '#9F3030');


for c = 1:num_conditions% Scatter individual participant data on top of bars
    % Scatter ADHD participant data
    scatter(adhd_positions(c) + zeros(size(adhd_means_hand(:, c))), adhd_means_hand(:, c), ...
        10, color_map('ADHD'), 'filled', 'MarkerFaceAlpha', 1);
    scatter(adhd_positions(c) + zeros(size(adhd_means_eye(:, c))), adhd_means_eye(:, c), ...
        10, color_map('ADHD'), 'filled', 'MarkerFaceAlpha', 0.1);
    % Scatter nonADHD participant data
    scatter(nonadhd_positions(c) + zeros(size(nonadhd_means_hand(:, c))), nonadhd_means_hand(:, c), ...
        10, color_map('nonADHD'), 'filled', 'MarkerFaceAlpha', 1);
    scatter(nonadhd_positions(c) + zeros(size(nonadhd_means_eye(:, c))), nonadhd_means_eye(:, c), ...
        10, color_map('nonADHD'), 'filled', 'MarkerFaceAlpha', 0.1);
end

xticks(x_positions);
xticklabels(condition_labels);
ylabel('Average Reaction Time (seconds)');
title(plot_title);
legend([adhd_hand_bar, nonadhd_hand_bar], {'ADHD', 'nonADHD'}, 'Location', 'northeast', 'Box', 'off');

%grid on;
hold off;
if safe
    %set(gcf, 'Position', [50, 50, 1200, 600]);
    saveas(gcf, fullfile(comparison_results_folder, safe_name));
end

end
