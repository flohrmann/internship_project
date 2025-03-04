function groupQuestionnaire(quest_table, group_labels)
    % group_diagnostic_analysis - Analyzes diagnostic criteria fulfillment and group statistics
    % Performs statistical tests and returns results in a table.
    %
    % Inputs:
    %   quest_table   - Table containing participant responses to 18 questions
    %   group_labels  - Cell array containing 'ADHD' or 'nonADHD' labels
    %
    % Outputs:
    %   A table summarizing median, mean, SEM, significance of group differences.

    % Convert group_labels to categorical for indexing
    group_labels = categorical(group_labels);

    % Define Part A and Total Score Calculation
    partA_questions = {'Question1', 'Question2', 'Question3', 'Question4', 'Question5', 'Question6'};
    total_questions = quest_table(:, 2:end); % All question columns

    % Calculate scores
    partA_scores = sum(table2array(quest_table(:, partA_questions)), 2);
    total_scores = sum(table2array(total_questions), 2);

    % Add scores to table
    quest_table.PartAScore = partA_scores;
    quest_table.TotalScore = total_scores;

    % Apply Diagnostic Criteria
    partA_diagnosis = partA_scores >= 14;
    total_diagnosis = total_scores >= 40;
    
    quest_table.MeetsPartACriteria = partA_diagnosis;
    quest_table.MeetsTotalCriteria = total_diagnosis;
    
    % Split into groups
    ADHD_idx = group_labels == "ADHD";
    nonADHD_idx = group_labels == "nonADHD";

    % Extract group data
    partA_ADHD = partA_scores(ADHD_idx);
    partA_nonADHD = partA_scores(nonADHD_idx);
    total_ADHD = total_scores(ADHD_idx);
    total_nonADHD = total_scores(nonADHD_idx);

    % Compute statistics for each group
    stats = table;
    groups = {'ADHD', 'nonADHD'};

    for i = 1:length(groups)
        if i == 1
            partA_group = partA_ADHD;
            total_group = total_ADHD;
        else
            partA_group = partA_nonADHD;
            total_group = total_nonADHD;
        end
        
        % Compute statistics
        stats.Group{i} = groups{i};
        stats.PartA_Median(i) = median(partA_group);
        stats.PartA_Mean(i) = mean(partA_group);
        stats.PartA_SEM(i) = std(partA_group) / sqrt(length(partA_group));
        
        stats.Total_Median(i) = median(total_group);
        stats.Total_Mean(i) = mean(total_group);
        stats.Total_SEM(i) = std(total_group) / sqrt(length(total_group));
        
        % Count participants meeting criteria
        stats.Meets_PartA_Criteria(i) = sum(partA_group >= 14);
        stats.Meets_Total_Criteria(i) = sum(total_group >= 40);
    end

    % Perform statistical tests
%     % Check normality (Shapiro-Wilk test)
%     [h_partA_ADHD, p_partA_ADHD] = swtest(partA_ADHD);
%     [h_partA_nonADHD, p_partA_nonADHD] = swtest(partA_nonADHD);
%     [h_total_ADHD, p_total_ADHD] = swtest(total_ADHD);
%     [h_total_nonADHD, p_total_nonADHD] = swtest(total_nonADHD);

    % Determine appropriate test
%     if p_partA_ADHD > 0.05 && p_partA_nonADHD > 0.05
        p_partA = permutationTestWelch(partA_ADHD, partA_nonADHD);
        stats.PartA_pValue(1) = p_partA;
%     else
%         p_partA = ranksum(partA_ADHD, partA_nonADHD);
%         stats.PartA_pValue(1) = p_partA;
%     end

%     if p_total_ADHD > 0.05 && p_total_nonADHD > 0.05
        p_total = permutationTestWelch(total_ADHD, total_nonADHD);
        stats.Total_pValue(1) = p_total;
%     else
%         p_total = ranksum(total_ADHD, total_nonADHD);
%         stats.Total_pValue(1) = p_total;
%     end

    % Fill p-values for second row (nonADHD group)
    stats.PartA_pValue(2) = NaN;
    stats.Total_pValue(2) = NaN;

    % Display results
    disp('Summary of ADHD and nonADHD Group Statistics:');
    disp(stats);

    % Return results as a table
    assignin('base', 'group_stats', stats);
end

function [h, p] = swtest(x)
    % Shapiro-Wilk test for normality
    % Returns h = 1 if data is not normal (p < 0.05), h = 0 if normal
    if length(x) < 3
        h = 0; p = 1; % Too few samples for normality test
    else
        [h, p] = swtest(x); % MATLAB’s Shapiro-Wilk function (from File Exchange)
    end
end
