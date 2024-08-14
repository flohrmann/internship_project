function data_struct = processEyeRTs(data_struct)
    % Function to process eye reaction times (RTs) in the data_struct
    % It selects the faster RT between the left and right eyes and ignores cases where RT is zero.
    %
    % Input: 
    %   data_struct - Structure array containing the fields rt_right and rt_left for eye RTs
    %
    % Output:
    %   data_struct - Structure array with an additional field rt_eye containing the processed RTs
    
    % Loop through each participant's data
    for i = 1:length(data_struct)
        % Extract right and left eye RTs
        rt_right = data_struct(i).rt_right;
        rt_left = data_struct(i).rt_left;
        
        % Initialize the processed eye RT array
        rt_eye = nan(size(rt_right));  % Initialize with NaN
        
        % Loop through each trial
        for j = 1:length(rt_right)
            % Check if both RTs are non-zero
            if rt_right(j) > 0 && rt_left(j) > 0
                % Take the faster (smaller) RT
                rt_eye(j) = min(rt_right(j), rt_left(j));
            elseif rt_right(j) > 0
                % If only right eye RT is non-zero, use it
                rt_eye(j) = rt_right(j);
            elseif rt_left(j) > 0
                % If only left eye RT is non-zero, use it
                rt_eye(j) = rt_left(j);
            else
                % If both RTs are zero, leave it as NaN
                rt_eye(j) = NaN;
            end
        end
        
        % Assign the processed RT to the data_struct
        data_struct(i).rt_eye = rt_eye;
    end
end
