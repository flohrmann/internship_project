function getGroupStats(data_struct, group_labels)
    adhd_idx = strcmp(group_labels, 'ADHD');
    nonadhd_idx = strcmp(group_labels, 'nonADHD');

    %  ages
    ages = arrayfun(@(x) x.age, data_struct);
    adhd_ages = ages(adhd_idx);
    nonadhd_ages = ages(nonadhd_idx);

    stats.Age_Median = [median(adhd_ages); median(nonadhd_ages)];
    stats.Age_Mean = [mean(adhd_ages); mean(nonadhd_ages)];
    stats.Age_SEM = [std(adhd_ages)/sqrt(length(adhd_ages)); std(nonadhd_ages)/sqrt(length(nonadhd_ages))];

    % gender 
    genders = arrayfun(@(x) x.gender, data_struct, 'UniformOutput', false);
    genders = [genders{:}]; % flatten cell array
    adhd_gender = genders(adhd_idx);
    nonadhd_gender = genders(nonadhd_idx);

    [adhd_unique_genders, ~, idx] = unique(adhd_gender);
    adhd_gender_counts = accumarray(idx, 1);
    [nonadhd_unique_genders, ~, idx] = unique(nonadhd_gender);
    nonadhd_gender_counts = accumarray(idx, 1);

    % gaming 
%     gaming = arrayfun(@(x) x.hGame, data_struct, 'UniformOutput', false);
%     gaming = [gaming{:}]; 
%     adhd_gaming = gaming(adhd_idx);
%     nonadhd_gaming = gaming(nonadhd_idx);
% 
%     [adhd_unique_gaming, ~, idx] = unique(adhd_gaming);
%     adhd_gaming_counts = accumarray(idx, 1);
%     [nonadhd_unique_gaming, ~, idx] = unique(nonadhd_gaming);
%     nonadhd_gaming_counts = accumarray(idx, 1);
% 
%     
    
    
    %%
    % summary table
    summary_table = table(["ADHD"; "nonADHD"], stats.Age_Median, stats.Age_Mean, stats.Age_SEM, ...
        'VariableNames', {'Group', 'Age_Median', 'Age_Mean', 'Age_SEM'});

    disp('Summary of Group Statistics:');
    disp(summary_table);
    fprintf('\\begin{table}[h]\n\\centering\n\\caption{Age Statistics for ADHD and nonADHD Groups}\n\\label{tab:age_stats}\n');
    fprintf('\\begin{tabular}{lccc}\n\\toprule\n');
    fprintf('Group & Age Median & Age Mean & Age SEM \\\\\n\\midrule\n');
    fprintf('ADHD & %.2f & %.2f & %.2f \\\\\n', stats.Age_Median(1), stats.Age_Mean(1), stats.Age_SEM(1));
    fprintf('nonADHD & %.2f & %.2f & %.2f \\\\\n', stats.Age_Median(2), stats.Age_Mean(2), stats.Age_SEM(2));
    fprintf('\\bottomrule\n\\end{tabular}\n\\end{table}\n\n');

    
    %% 
    
    %  Gender 
    fprintf('\nGender Distribution:\n');
    gender_table = table(adhd_unique_genders', adhd_gender_counts,nonadhd_gender_counts, 'VariableNames', {'Gender', 'ADHD', 'nonADHD'});
    disp(gender_table);

%     %  Gaming 
%     fprintf('\nGaming Preferences for ADHD Group:\n');
%     gaming_table = table(adhd_unique_gaming', adhd_gaming_counts, nonadhd_gaming_counts, 'VariableNames', {'Gaming', 'ADHD', 'nonADHD'});
%     disp(gaming_table);

end
