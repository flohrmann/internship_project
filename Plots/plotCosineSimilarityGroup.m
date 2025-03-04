function plotCosineSimilarityGroup(data, group_labels, ids, all_saccades, all_fixations, ...
                          unique_conditions, condition_labels, color_map, color_map_individual,  comp_results_fix)
                      
                      
% Initialize variables
numParticipants = length(data); % Number of participants
group_cos_sims_cond = cell(length(unique_conditions), numParticipants); % Cell array to store cosine similarities for all participants

% Loop through each participant
for p = 1:numParticipants
    trial_conditions = data(p).Condition;
    saccades = all_saccades(p).saccades;
    fixations = all_fixations(p).fixations;

    num_trials = length(saccades);
    x_s_start = []; y_s_start = []; % sacade start x/y 
    x_s_end = []; y_s_end = [];     % sacade end x/y 
    x_target = []; y_target = [];   % target x/y 

    % Process each trial for the participant
    for a = 1:num_trials
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

    % Compute vectors
    x_s_end_centered  = x_s_end - x_s_start;
    y_s_end_centered  = y_s_end - y_s_start;
    x_target_centered = x_target - x_s_start;
    y_target_centered = y_target - y_s_start;

    s_end_vector = [x_s_end_centered(:), y_s_end_centered(:)];
    target_vector = [x_target_centered(:), y_target_centered(:)];
    
    % cosine similarity per condition
    for t = 1:num_trials
        s_end = s_end_vector(t, :);
        target = target_vector(t, :);

        % Ensure vectors are valid
        if all(~isnan(s_end)) && all(~isnan(target))
            cos_sim = dot(s_end, target) / (norm(s_end) * norm(target));
        else
            cos_sim = NaN; % Handle invalid trials
        end
        
        % Find condition index and store cosine similarity
        conditionIdx = find(strcmp(unique_conditions, trial_conditions{t}));
        if ~isempty(conditionIdx)
            try
                group_cos_sims_cond{conditionIdx, p} = [group_cos_sims_cond{conditionIdx, p}, cos_sim];
            catch
                group_cos_sims_cond{conditionIdx, p} = [group_cos_sims_cond{conditionIdx, p}, NaN];
            end
        end
    end
end


av_cos_sim = zeros(numParticipants, length(unique_conditions));
for p=1:numParticipants
        avg_cos_sim(p,:) = [nanmedian(group_cos_sims_cond{1,p}), nanmedian(group_cos_sims_cond{2,p}), ...
                            nanmedian(group_cos_sims_cond{3,p}), nanmedian(group_cos_sims_cond{4,p})];           
end

% get data into format for plot
adhd_indices = find(strcmp(group_labels, 'ADHD'));
nonadhd_indices = find(strcmp(group_labels, 'nonADHD'));
for idx = 1:numel(adhd_indices)
    adhd_condition_avgs(idx, :) = avg_cos_sim(adhd_indices(idx),:); % Extract the participant index & data                             
end

for idx = 1:numel(nonadhd_indices)
    nonadhd_condition_avgs(idx, :) = avg_cos_sim(nonadhd_indices(idx),:); % Extract the participant index & data                                 
end

adhd_median = median(adhd_condition_avgs, 1, 'omitnan');
adhd_sem = std(adhd_condition_avgs, 0, 1, 'omitnan')./ sqrt(size(adhd_condition_avgs, 1));
nonAdhd_median = median(nonadhd_condition_avgs, 1, 'omitnan');
nonAdhd_sem = std(nonadhd_condition_avgs, 0, 1, 'omitnan')./ sqrt(size(nonadhd_condition_avgs, 1));

plotADHDnonADHDandDiff('Cosine Similarity of 1st Saccade of Trials',...
    adhd_condition_avgs, adhd_median, adhd_sem, 'ADHD', 'southeast', ...
    nonadhd_condition_avgs, nonAdhd_median, nonAdhd_sem, 'nonADHD', 'southeast', ...
    ids, condition_labels, 'Cosine Similarity', group_labels, unique_conditions, color_map, color_map_individual, ...
    fullfile(comp_results_fix, '04_cosine_similarity_allinone_median.png'));



