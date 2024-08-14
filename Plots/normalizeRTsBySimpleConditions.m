function data = normalizeRTsBySimpleConditions(data, condition_A, condition_B)
    % Function to normalize reaction times (RTs) by the mean RT of two specified conditions
    % and store the normalized RTs back into the data_struct.
    % Inputs:
    %   data_struct: struct array containing the data for each observer
    %   condition_A: string, name of the first simple condition (e.g., 'Asimple')
    %   condition_B: string, name of the second simple condition (e.g., 'Bsimple')
    % Output:
    %   data_struct: struct array with an added field 'normalized_rt' for each observer
    
    % Loop through each observer in the data_struct
    for i = 1:length(data)
        rt        = data(i).rt;
        rt_eye    = data(i).rt_eye;
        condition = data(i).Condition; 
        
        % Find the indices for condition_A and condition_B
        idx_A_simple = strcmp(condition, condition_A);
        idx_B_simple = strcmp(condition, condition_B);
        
        % Calculate the mean RT for condition_A and condition_B
        mean_rt_simple     = (nanmean(rt(idx_A_simple))     + nanmean(rt(idx_B_simple)))     / 2;
        mean_rt_eye_simple = (nanmean(rt_eye(idx_A_simple)) + nanmean(rt_eye(idx_B_simple))) / 2;

        % Normalize the RTs for this observer
        normalized_rt     = rt / mean_rt_simple;
        normalized_rt_eye = rt_eye / mean_rt_eye_simple;

        % Store the normalized RTs back into the data_struct
        data(i).normalized_rt     = normalized_rt;
        data(i).normalized_rt_eye = normalized_rt_eye;
    end
    
    % Display a message indicating that the normalization is complete
    disp('Normalization of RTs by simple conditions and storage in data_struct is complete.');
end
