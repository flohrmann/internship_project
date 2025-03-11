function [avg_adhd, avg_nonadhd] = extractPupilData(data, outlier_threshold, usecase)

% Initialize tables
avg_adhd = table([], [], [], [], [],...
    [], [], [], [], ...
    [], [], [], [], ...
    [], [], [], [], ...
    [], [], [], [], [], [], [], [],...
    [], [], [], [], [], [], [], [],...
    'VariableNames', {'id', 'a', 'a_simple', 'b', 'b_simple', ...
    'a_max', 'as_max', 'b_max', 'bs_max', ...
    'a_min', 'as_min', 'b_min', 'bs_min', ...
    'a_mean', 'as_mean', 'b_mean', 'bs_mean', ...
    'rt_a', 'rtv_a', 'rt_as', 'rtv_as', 'rt_b', 'rtv_b', 'rt_bs', 'rtv_bs', ...
    'rt_eye_a', 'rtv_eye_a', 'rt_eye_as', 'rtv_eye_as', 'rt_eye_b', 'rtv_eye_b', 'rt_eye_bs', 'rtv_eye_bs'});

avg_nonadhd = table([], [], [], [], [], ...
    [], [], [], [], ...
    [], [], [], [], ...
    [], [], [], [], ...
    [], [], [], [], [], [], [], [],...
    [], [], [], [], [], [], [], [],...
    'VariableNames', {'id', 'a', 'a_simple', 'b', 'b_simple', ...
    'a_max', 'as_max', 'b_max', 'bs_max', ...
    'a_min', 'as_min', 'b_min', 'bs_min', ...
    'a_mean', 'as_mean', 'b_mean', 'bs_mean', ...
    'rt_a', 'rtv_a', 'rt_as', 'rtv_as', 'rt_b', 'rtv_b', 'rt_bs', 'rtv_bs', ...
    'rt_eye_a', 'rtv_eye_a', 'rt_eye_as', 'rtv_eye_as', 'rt_eye_b', 'rtv_eye_b', 'rt_eye_bs', 'rtv_eye_bs'});

% helper function to remove outliers
remove_outliers = @(x) replace_outliers_with_nan(x, outlier_threshold);

for i = 1:length(data)
        diam = data(i).pupilDiamNorm;
        id = data(i).id;
        
        % Pupil data without outlier removal
        a = nanmean([diam.a.beforeMean, diam.a.afterMean], 1);
        a_simple = nanmean([diam.a_simple.beforeMean, diam.a_simple.afterMean], 1);
        b = nanmean([diam.b.beforeMean, diam.b.afterMean], 1);
        b_simple = nanmean([diam.b_simple.beforeMean, diam.b_simple.afterMean], 1);
        
        % Calculate min, max, and mean for pupil data without outlier removal
        as_max = nanmean(diam.a_simple.maxMean, 1);
        a_max = nanmean(diam.a.maxMean, 1);
        bs_max = nanmean(diam.b_simple.maxMean, 1);
        b_max = nanmean(diam.b.maxMean, 1);
        
        as_min = nanmean(diam.a_simple.minMean, 1);
        a_min = nanmean(diam.a.minMean, 1);
        bs_min = nanmean(diam.b_simple.minMean, 1);
        b_min = nanmean(diam.b.minMean, 1);
        
        as_mean = nanmean(diam.a_simple.meanMean, 1);
        a_mean = nanmean(diam.a.meanMean, 1);
        bs_mean = nanmean(diam.b_simple.meanMean, 1);
        b_mean = nanmean(diam.b.meanMean, 1);
        
        % Reaction Time (RT) and Reaction Time Variability (RTV) with outlier removal
        rt_a = nanmean(remove_outliers(data(i).rt(strcmp(data(i).Condition, 'a'))));
        rtv_a = nanstd(remove_outliers(data(i).rt(strcmp(data(i).Condition, 'a'))));
        rt_eye_a = nanmean(remove_outliers(data(i).rt_eye(strcmp(data(i).Condition, 'a'))));
        rtv_eye_a = nanstd(remove_outliers(data(i).rt_eye(strcmp(data(i).Condition, 'a'))));
        
        rt_as = nanmean(remove_outliers(data(i).rt(strcmp(data(i).Condition, 'a_simple'))));
        rtv_as = nanstd(remove_outliers(data(i).rt(strcmp(data(i).Condition, 'a_simple'))));
        rt_eye_as = nanmean(remove_outliers(data(i).rt_eye(strcmp(data(i).Condition, 'a_simple'))));
        rtv_eye_as = nanstd(remove_outliers(data(i).rt_eye(strcmp(data(i).Condition, 'a_simple'))));
        
        rt_b = nanmean(remove_outliers(data(i).rt(strcmp(data(i).Condition, 'b'))));
        rtv_b = nanstd(remove_outliers(data(i).rt(strcmp(data(i).Condition, 'b'))));
        rt_eye_b = nanmean(remove_outliers(data(i).rt_eye(strcmp(data(i).Condition, 'b'))));
        rtv_eye_b = nanstd(remove_outliers(data(i).rt_eye(strcmp(data(i).Condition, 'b'))));
        
        rt_bs = nanmean(remove_outliers(data(i).rt(strcmp(data(i).Condition, 'b_simple'))));
        rtv_bs = nanstd(remove_outliers(data(i).rt(strcmp(data(i).Condition, 'b_simple'))));
        rt_eye_bs = nanmean(remove_outliers(data(i).rt_eye(strcmp(data(i).Condition, 'b_simple'))));
        rtv_eye_bs = nanstd(remove_outliers(data(i).rt_eye(strcmp(data(i).Condition, 'b_simple'))));
        
        % Add the data to the appropriate table
        new_row = {id, a, a_simple, b, b_simple, ...
            a_max, as_max, b_max, bs_max, ...
            a_min, as_min, b_min, bs_min, ...
            a_mean, as_mean, b_mean, bs_mean, ...
            rt_a, rtv_a, rt_as, rtv_as, rt_b, rtv_b, rt_bs, rtv_bs, ...
            rt_eye_a, rtv_eye_a, rt_eye_as, rtv_eye_as, rt_eye_b, rtv_eye_b, rt_eye_bs, rtv_eye_bs};
        
        if strcmp(usecase, 'ADHD')
            if strcmp(data(i).group, 'ADHD')
                avg_adhd = [avg_adhd; new_row];
            else
                avg_nonadhd = [avg_nonadhd; new_row];
            end
        else % MS
            if data(i).id == 1 || data(i).id == 13
                avg_adhd = [avg_adhd; new_row];
            elseif strcmp(data(i).group, 'nonADHD')
                avg_nonadhd = [avg_nonadhd; new_row];
            else
                % ignore adhd participants for now
            end
        end
end
end

% Helper function to detect and replace outliers with NaN
function cleaned_data = replace_outliers_with_nan(filtered_data, outlier_threshold)
derivative_data = diff(filtered_data); % Calculate derivative
mean_derivative = mean(derivative_data, 'omitnan');
std_derivative = std(derivative_data, 'omitnan');

% Detect outliers based on the threshold
outliers = abs(derivative_data) > (mean_derivative + outlier_threshold * std_derivative);

% Mark outliers in the original data (offset by one because of diff)
cleaned_data = filtered_data;
cleaned_data(2:end) = filtered_data(2:end) .* ~outliers;
cleaned_data(cleaned_data == 0) = NaN; % Replace outliers with NaN
end
