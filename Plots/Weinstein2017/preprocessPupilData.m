function [processed_before, processed_after, processed_both, steps_before, steps_after, steps_both]  = preprocessPupilData(pupil_data_before, pupil_data_after, baseline)
    % Interpolate only internal NaNs for pupil_data_before
    pupil_data_before_new = interpolateInternalNaNs(pupil_data_before);

    % Interpolate only internal NaNs for pupil_data_after
    pupil_data_after_new = interpolateInternalNaNs(pupil_data_after);
    
    % Check if both lists have enough valid data points (e.g., more than 12)
    if numel(pupil_data_before_new(~isnan(pupil_data_before_new))) > 12 && ...
       numel(pupil_data_after_new(~isnan(pupil_data_after_new))) > 12
       
        % Normalize pupil_data_before by baseline
        [processed_before, steps_before] = processSingleList(pupil_data_before_new, baseline);
        
        % Normalize pupil_data_after by baseline
        [processed_after, steps_after] = processSingleList(pupil_data_after_new, baseline);
        
        % Normalize both at once by baseline
        current_pupil_data = [pupil_data_before_new, pupil_data_after_new];
        [processed_both, steps_both] = processSingleList(current_pupil_data, baseline);

    else
        % If either list is too short, disregard both
        processed_before = [];
        processed_after = [];
        processed_both = [];
        steps_before = struct();
        steps_after = struct();
        steps_both = struct();
    end
end

% Wainstein 2017
% "data of each participant were baseline-adjusted 
% and smoothed by a bandpass Butterworth filter between 0.025 Hz and 4 Hz. 
% outliers, defined as periods of pupil change (derivative function) higher than 3 standard errors from the mean were discarded
% all trials with more than 50% of missing data (due to blinks or outliers) were not considered in the analysis
% pupil timeseries was normalized by means of a z-score, separately for each trial."

% Helper function to process each list (baseline correction, filtering, z-score)
function [processed_data, processed_steps] = processSingleList(pupil_data, baseline)
    missing_data_limit = 0.5;
    outlier_threshold  = 3;
    
    nan_mask = isnan(pupil_data);  % Find NaNs

    %% Step 1: Remove NaNs for processing (Butterworth filter doesn’t work with NaNs)
    data_without_nans = pupil_data(~nan_mask);
    processed_steps.data_without_nans = data_without_nans;

    %% Step 2: Apply baseline correction
    data_corrected = data_without_nans - baseline;
    processed_steps.data_corrected = data_corrected;

    %% Step 3: Apply bandpass filter to clean noise (Butterworth filter) [skip]
    %[b, a] = butter(2, [0.025, 4] / (60 / 2), 'bandpass');  % 60Hz sampling rate
    %filtered_data = filtfilt(b, a, data_corrected);
    %processed_steps.filtered_data = filtered_data;
    filtered_data= data_corrected;
    
    %% Step 4: Detect outliers based on the derivative of the pupil data
    derivative_data = diff(filtered_data);  % Calculate derivative
    mean_derivative = mean(derivative_data, 'omitnan');
    std_derivative = std(derivative_data, 'omitnan');
    outliers = abs(derivative_data) > (mean_derivative + outlier_threshold * std_derivative);
    processed_steps.derivative_data = derivative_data;
    processed_steps.outliers = outliers;

    %% Step 5: Mark these outliers in the original data (offset by one because of diff)
    filtered_data(2:end) = filtered_data(2:end) .* ~outliers;
    outlier_mask = isnan(filtered_data);
    processed_steps.filtered_data_with_outliers_removed = filtered_data;

    % Check if more than 50% of data is missing or marked as outliers
    if mean(outlier_mask) > missing_data_limit
        processed_data = nan(size(pupil_data));
        processed_steps.final_data = processed_data;
        warning('Trial discarded due to more than 50% missing/outlier data.');
        return;
    end

    %% Step 6: Interpolate data where outliers were removed
    processed_data_no_outliers = fillmissing(filtered_data, 'spline');
    processed_steps.interpolated_data = processed_data_no_outliers;

    %% Step 7: Z-score normalization of the interpolated data
    normalized_data = zscore(processed_data_no_outliers);
    processed_steps.normalized_data = normalized_data;

%     %% Step 8: Reintroduce NaNs at the original positions
%     processed_data = nan(size(pupil_data));  % Initialize with NaNs
%     processed_data(~nan_mask) = normalized_data;  % Insert processed values back into non-NaN positions
%     processed_steps.final_data = processed_data;

    % Optiona
    % save('processed_steps.mat', 'processed_steps');
end



% Helper function to interpolate only internal NaNs (between two valid points)
function interpolated_data = interpolateInternalNaNs(pupil_data)
    % Find the first and last non-NaN values
    first_valid_idx = find(~isnan(pupil_data), 1, 'first');
    last_valid_idx = find(~isnan(pupil_data), 1, 'last');
    
    % Only interpolate NaNs that are between two valid data points
    interpolated_data = pupil_data;  % Start with the original data
    
    % Use interpolation only on the range between first_valid_idx and last_valid_idx
    interpolated_range = first_valid_idx:last_valid_idx;
    
    % Interpolate internal NaNs within the valid range
    interpolated_data(interpolated_range) = fillmissing(pupil_data(interpolated_range), 'spline');
end
