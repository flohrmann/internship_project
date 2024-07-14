function plotRT(trial_results)
% Extract relevant data from the table
conditions = trial_results.Condition;
rt = trial_results.rt;

% Convert conditions to categorical data for coloring
unique_conditions = unique(conditions);
colors = lines(length(unique_conditions)); % Generate distinct colors

% Create a figure for the plot
figure;
hold on;

% Plot RTs colored by condition
for i = 1:length(unique_conditions)
    % Find the indices for the current condition
    condition_indices = strcmp(conditions, unique_conditions{i});
    
    % Extract the RTs for the current condition
    rt_condition = rt(condition_indices);
    
    % Plot the RTs for the current condition
    scatter(find(condition_indices), rt_condition, [], colors(i,:), 'filled');
end

% Add labels and legend
xlabel('Trial Number');
ylabel('Reaction Time (s)');
title('Reaction Time Colored by Condition');
legend(unique_conditions, 'Location', 'best');
hold off;
end