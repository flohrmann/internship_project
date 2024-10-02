function compareGroupCorrelations(R_group_condition, groups, data, conditions)
    variables = {'RT_ButtonPress', 'RT_Eye', 'Accuracy'};
    numConditions = length(conditions);

    % Initialize to store z-scores and p-values for each condition and variable pair
    z_scores = cell(numConditions, 1);
    p_values = cell(numConditions, 1);
    
    % Loop through each condition
    for c = 1:numConditions
        condition = conditions{c};
        
        % Initialize matrices to store z-scores and p-values
        z_matrix = zeros(3, 3);
        p_matrix = zeros(3, 3);
        
        % Extract the correlation matrices for ADHD and non-ADHD groups
        R_ADHD = R_group_condition.ADHD{c};
        R_nonADHD = R_group_condition.nonADHD{c};
        
        % Get sample sizes for both groups
        n_ADHD = sum(strcmp({data.group}, 'ADHD') & strcmp({data.Condition}, condition));
        n_nonADHD = sum(strcmp({data.group}, 'non-ADHD') & strcmp({data.Condition}, condition));

        % Perform Fisher's r-to-z transformation for each variable pair
        for i = 1:3
            for j = i+1:3
                % Fisher r-to-z transformation
                z_ADHD = 0.5 * log((1 + R_ADHD(i, j)) / (1 - R_ADHD(i, j)));
                z_nonADHD = 0.5 * log((1 + R_nonADHD(i, j)) / (1 - R_nonADHD(i, j)));
                
                % Calculate the z-score difference between ADHD and non-ADHD
                z_diff = (z_ADHD - z_nonADHD) / sqrt(1/(n_ADHD - 3) + 1/(n_nonADHD - 3));
                
                % Compute the p-value from the z-score
                p_value = 2 * (1 - normcdf(abs(z_diff))); % Two-tailed test
                
                % Store the z-score and p-value
                z_matrix(i, j) = z_diff;
                p_matrix(i, j) = p_value;
            end
        end
        
        % Store the z-scores and p-values for this condition
        z_scores{c} = z_matrix;
        p_values{c} = p_matrix;
        
        % Display the z-scores and p-values for this condition
        disp(['Condition: ', condition]);
        disp('Z-Scores Matrix:');
        disp(z_matrix);
        disp('P-Values Matrix:');
        disp(p_matrix);
        disp('------------------------------------------');
    end
end
