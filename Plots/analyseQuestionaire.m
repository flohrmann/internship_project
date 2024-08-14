% Define the folder paths and corresponding file names
folders = {
    'C:\Users\flohrmann\Documents\Results\1_20240714_140048'; % alex
    'C:\Users\flohrmann\Documents\Results\2_20240726_191755'; % mara adhd
    'C:\Users\flohrmann\Documents\Results\3_20240805_105213'; % tilo adhd
    'C:\Users\flohrmann\Documents\Results\4_20240811_131601'; % anu adhd
};

% Define the classification of the folders (ADHD or non-ADHD)
group_labels = {
    'non_ADHD'; % alex
    'ADHD';     % mara adhd
    'ADHD';     % tilo adhd
    'ADHD';     % anu adhd
};

% Initialize structures to hold the data
data.ADHD = table(); % Initialize as empty table to handle empty data appropriately
data.non_ADHD = table(); % Initialize as empty table to handle empty data appropriately

% Create a generic set of column names including the ID column
num_questions = 18;
column_names = [{'ID'}, arrayfun(@(x) sprintf('Question%d', x), 1:num_questions, 'UniformOutput', false)];

% Define the response mapping
response_mapping = containers.Map({'Never', 'Rarely', 'Sometimes', 'Often', 'Very Often'}, [1, 2, 3, 4, 5]);

% Load the data from each folder
for i = 1:length(folders)
    folder = folders{i};
    label = group_labels{i};
    
    % Extract the ID from the folder name (the number after the last backslash)
    [~, folder_name] = fileparts(folder);
    id = sscanf(folder_name, '%d_%*s'); % Extract the leading number (ID)
    
    % Find the CSV file that matches the pattern 'questionnaire*.csv' in the folder
    file = dir(fullfile(folder, 'questionnaire*.csv'));
    
    % Ensure that the file exists before trying to load it
    if ~isempty(file)
        filepath = fullfile(folder, file.name);
        file_data = readtable(filepath);
        
        % Replace the response words with their corresponding numeric values
        person_responses = cellfun(@(x) response_mapping(x), file_data.Response)';
        
        % Include the ID in the first column
        person_data = [id, person_responses];
        
        % Convert to a table and append to the appropriate group
        if isempty(data.(label))
            % Create the table with ID and generic column names (Question1, Question2, ...)
            data.(label) = array2table(person_data, 'VariableNames', column_names);
        else
            % Convert responses to a table and append as a new row
            person_table = array2table(person_data, 'VariableNames', column_names);
            data.(label) = [data.(label); person_table];
        end
    else
        warning('No questionnaire file found in folder: %s', folder);
    end
end

% Example: Accessing ADHD and non-ADHD data
adhd_data = data.ADHD;
non_adhd_data = data.non_ADHD;

% Display a summary of the data
fprintf('Loaded %d ADHD entries and %d non-ADHD entries.\n', height(adhd_data), height(non_adhd_data));

% Convert responses to numeric for easier analysis and plotting
adhd_numeric = adhd_data{:, 2:end};  % Skip the ID column
non_adhd_numeric = non_adhd_data{:, 2:end};  % Skip the ID column

% Plotting

%% 1. Bar Plot for response distribution for each question
figure;
for q = 1:num_questions
    subplot(3, 6, q); % Arrange in a 3x6 grid
    adhd_counts = histcounts(adhd_numeric(:, q), 1:6); % Use histcounts for numeric data
    non_adhd_counts = histcounts(non_adhd_numeric(:, q), 1:6); % Use histcounts for numeric data
    bar([adhd_counts; non_adhd_counts]', 'grouped');
    title(sprintf('Question %d', q));
    xticks(1:5);
    xticklabels({'Never', 'Rarely', 'Sometimes', 'Often', 'Very Often'});
    ylabel('Frequency');
    legend({'ADHD', 'non-ADHD'}, 'Location', 'best');
end
sgtitle('Response Distribution for Each Question');

%% 2. Box Plot for comparing response distribution by group
figure;
group_labels = [repmat({'ADHD'}, height(adhd_data), 1); repmat({'non-ADHD'}, height(non_adhd_data), 1)];
all_data = [adhd_numeric; non_adhd_numeric];
group = repelem(group_labels, num_questions, 1);
boxplot(all_data(:), group);
ylabel('Response (1=Never, 5=Very Often)');
title('Response Distribution by Group');

%% 3. Pie Chart for overall response proportions for ADHD group
figure;
adhd_counts_total = histcounts(adhd_numeric(:), 0.5:1:5.5); % Count occurrences of each response category
pie(adhd_counts_total);
%legend({'Never', 'Rarely', 'Sometimes', 'Often', 'Very Often'}, 'Location', 'bestoutside');
title('Overall Response Proportions for ADHD Group');

%% 4. Pie Chart for overall response proportions for non-ADHD group
figure;
non_adhd_counts_total = histcounts(non_adhd_numeric(:), 0.5:1:5.5); % Count occurrences of each response category
pie(non_adhd_counts_total);
%legend({'Never', 'Rarely', 'Sometimes', 'Often', 'Very Often'}, 'Location', 'bestoutside');
title('Overall Response Proportions for non-ADHD Group');

%% Answers per question with a unique color for each participant
% Combine ADHD and non-ADHD data for plotting
all_data = [adhd_data; non_adhd_data];
ids = all_data{:, 1}; % Extract the ID column
all_numeric = all_data{:, 2:end}; % Extract the numeric responses

% Get the number of participants and questions
[num_participants, num_questions] = size(all_numeric);

% Generate a color map with a unique color for each participant
colors = lines(num_participants);

% Create the plot
figure;
hold on;
for i = 1:num_participants
    plot(1:num_questions, all_numeric(i, :), '-o', 'Color', colors(i, :), 'LineWidth', 1.5);
end
hold off;

xlabel('Question Number');
ylabel('Answer');
title('Participant Answers by Question');
xticks(1:num_questions);
xlim([1 num_questions]);
yticks(1:5);
yticklabels({'Never', 'Rarely', 'Sometimes', 'Often', 'Very Often'});
legend(arrayfun(@(x) sprintf('ID %d', ids(x)), 1:num_participants, 'UniformOutput', false), ...
    'Location', 'eastoutside');
grid on;


%% Same but only points and colored by group (ADHD vs Non-ADHD)
% Assuming `adhd_data` and `non_adhd_data` are the tables with the loaded data

% Extract the numeric responses for ADHD and non-ADHD groups
adhd_numeric = adhd_data{:, 2:end};  % Skip the ID column
non_adhd_numeric = non_adhd_data{:, 2:end};  % Skip the ID column

% Get the number of participants and questions
num_questions = size(adhd_numeric, 2);

% Create the plot
figure;
hold on;

% Plot the non-ADHD data points in one color (e.g., blue)
for i = 1:size(non_adhd_numeric, 1)
    scatter(1:num_questions, non_adhd_numeric(i, :), 50, 'blue', 'filled', 'MarkerFaceAlpha', 0.5);
end

% Plot the ADHD data points in another color (e.g., red)
for i = 1:size(adhd_numeric, 1)
    scatter(1:num_questions, adhd_numeric(i, :), 50, 'red', 'filled', 'MarkerFaceAlpha', 0.5);
end
xlabel('Question Number');
ylabel('Response');
title('ADHD vs. Non-ADHD Responses by Question');
xticks(1:num_questions);
xlim([1 num_questions]);
yticks(1:5);
yticklabels({'Never', 'Rarely', 'Sometimes', 'Often', 'Very Often'});
legend({'Non-ADHD', 'ADHD'}, 'Location', 'eastoutside');
grid on;
hold off;


%% Means of answers per question per group
% Assuming `adhd_data` and `non_adhd_data` are the tables with the loaded data

% Extract the numeric responses for ADHD and non-ADHD groups
adhd_numeric = adhd_data{:, 2:end};  % Skip the ID column
non_adhd_numeric = non_adhd_data{:, 2:end};  % Skip the ID column

% Get the number of questions
num_questions = size(adhd_numeric, 2);

% Calculate the mean and standard deviation for each question
adhd_means = mean(adhd_numeric, 1, 'omitnan');
adhd_stds = std(adhd_numeric, 0, 1, 'omitnan');
non_adhd_means = mean(non_adhd_numeric, 1, 'omitnan');
non_adhd_stds = std(non_adhd_numeric, 0, 1, 'omitnan');

figure;
hold on;
% Plot the mean and standard deviation for non-ADHD group
errorbar(1:num_questions, non_adhd_means, non_adhd_stds, '-o', 'Color', 'blue', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'blue', 'MarkerSize', 8);

% Plot the mean and standard deviation for ADHD group
errorbar(1:num_questions, adhd_means, adhd_stds, '-o', 'Color', 'red', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'red', 'MarkerSize', 8);

xlabel('Question Number');
ylabel('Response');
title('Mean and Standard Deviation of Responses by Question');
xticks(1:num_questions);
xlim([1 num_questions]);
yticks(1:5);
yticklabels({'Never', 'Rarely', 'Sometimes', 'Often', 'Very Often'});
legend({'Non-ADHD', 'ADHD'}, 'Location', 'eastoutside');
grid on;
hold off;




