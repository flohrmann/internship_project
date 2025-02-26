function condition_diam_avg = analyzePupilDataWithPreprocessing(result_table, cut_data, conditions, base, id, analysis_folder, color_map, sr)

% Initialize containers for group-level data
condition_diam_avg = struct();  % Store mean data for each condition


% Normalize the diameter data relative to the baseline first


for i = 1:length(conditions)
    condition = conditions{i};
    condition_mask = strcmp(result_table.Condition, condition);
    data_before = result_table.DataPointsBefore(condition_mask, :);
    data_after = result_table.DataPointsAfter(condition_mask, :);
    
    % Normalize data relative to the baseline
    if contains(base, "BlankScreen")
        baseline_mean = result_table.DiamBothBlank3StdFilterMean(condition_mask, :);  % Using the mean of blank screen as baseline
        baseline_median = result_table.DiamBothBlank3StdFilterMedian(condition_mask, :);
    else % 'StartOfTrial'
        baseline_mean = result_table.DiamBothStart3StdFilterMean(condition_mask, :);  % Using the mean of trial start as baseline
        baseline_median = result_table.DiamBothStart3StdFilterMedian(condition_mask, :);
    end
    
    
    % Preprocess pupil data for each trial
    processed_data_before_mean_list = [];
    processed_data_after_mean_list = [];
    processed_data_before_median_list = [];
    processed_data_after_median_list = [];
    max_mean = [];
    max_median = [];
    min_mean    = [];
    min_median  = [];
    mean_mean   = [];
    mean_median = [];
    processed_both_mean_list   = [];
    steps_before_mean_list     = [];
    steps_before_mean_struct   = struct();

    steps_after_mean_list      = [];
    steps_both_mean_list       = [];
    processed_both_median_list   = [];
    steps_before_median_list     = [];
    steps_after_median_list      = [];
    steps_both_median_list       = [];
    
    
    
    for trial = 1:size(data_before, 1)
        pupil_data_before = data_before(trial,:);
        pupil_data_after = data_after(trial,:);
        
        % Preprocess the pupil data for before and after target found
        % Ignore autofilled NaNs to buffer length for this
        if ~isnan(baseline_mean(trial,1)) && ~isnan(baseline_median(trial,1))
            [processed_before_mean, processed_after_mean, processed_both_mean, steps_before_mean, steps_after_mean, steps_both_mean]      = preprocessPupilData(pupil_data_before, pupil_data_after, baseline_mean(trial,1));
            [processed_before_median, processed_after_median, processed_both_median, steps_before_median, steps_after_median, steps_both_median]  = preprocessPupilData(pupil_data_before, pupil_data_after, baseline_median(trial,1));
            
            
            % Combine preprocessed data
            processed_data_before_mean_list   = [processed_data_before_mean_list; processed_before_mean];
            processed_data_after_mean_list    = [processed_data_after_mean_list; processed_after_mean];
            processed_data_before_median_list = [processed_data_before_median_list; processed_before_median];
            processed_data_after_median_list  = [processed_data_after_median_list; processed_after_median];
            
            % Compute the maximum pupil dilation after the stimulus onset
            max_mean    = [max_mean; max(processed_after_mean)];
            max_median  = [max_median; max(processed_after_median)];
            min_mean    = [min_mean; min(processed_after_mean)];
            min_median  = [min_median; min(processed_after_median)];
            mean_mean   = [mean_mean; mean(processed_after_mean)];
            mean_median = [mean_median; mean(processed_after_median)];
            
            
            % extra steps for plotting
            processed_both_mean_list   = [processed_both_mean_list; processed_both_mean];
            steps_before_mean_list = [steps_before_mean_list; {steps_before_mean}];
            steps_before_mean_struct.(trial) = steps_before_mean;
            steps_after_mean_list      = [steps_after_mean_list; {steps_after_mean}];
            steps_both_mean_list       = [steps_both_mean_list; {steps_both_mean}];
            
            processed_both_median_list   = [processed_both_median_list; processed_both_median];
            steps_before_median_list     = [steps_before_median_list; {steps_before_median}];
            steps_after_median_list      = [steps_after_median_list; {steps_after_median}];
            steps_both_median_list       = [steps_both_median_list; {steps_both_median}];
            
        else
        end
    end
    condition_diam_avg.(condition).beforeMean = processed_data_before_mean_list;
    condition_diam_avg.(condition).afterMean  = processed_data_after_mean_list;
    condition_diam_avg.(condition).maxMean    = max_mean;
    condition_diam_avg.(condition).minMean    = min_mean;
    condition_diam_avg.(condition).meanMean   = mean_mean;
    
    condition_diam_avg.(condition).beforeMedian = processed_data_before_median_list;
    condition_diam_avg.(condition).afterMedian  = processed_data_after_median_list;
    condition_diam_avg.(condition).maxMedian    = max_median;
    condition_diam_avg.(condition).minMedian    = min_median;
    condition_diam_avg.(condition).meanMedian   = mean_median;
    
    % extra steps for plotting
    condition_diam_avg.(condition).processed_both_mean =  processed_both_mean_list;
    condition_diam_avg.(condition).steps_before_mean   =  steps_before_mean_list;
    condition_diam_avg.(condition).steps_after_mean    =  steps_after_mean_list;
    condition_diam_avg.(condition).steps_both_mean     = steps_both_mean_list;
    
    condition_diam_avg.(condition).processed_both_median =  processed_both_median_list;
    condition_diam_avg.(condition).steps_before_median   =  steps_before_median_list;
    condition_diam_avg.(condition).steps_after_median    = steps_after_median_list;
    condition_diam_avg.(condition).steps_both_median     =  steps_both_median_list;
end

save(fullfile(analysis_folder, strcat('\pupil_processed_found_stim_', base ,'.mat')), 'condition_diam_avg');

% Plot the normalized pupil responses before and after the target is found
%plotPupilResponses(condition_diam_avg, id, conditions, analysis_folder, color_map, sr, base)
end

