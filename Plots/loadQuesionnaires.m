function data = loadQuesionnaires(folders, group_labels)
% loads questionaire results from paths provided
% sepperates them by group_labels (ADHD, nonADHD)
% return struct with two tables: ADHD and nonADHD
%                                contains row per participants answers
%                                answers are converted to numbers (1 never: 5 very often)

data.ADHD = table(); % Initialize as empty table to handle empty data appropriately
data.nonADHD = table(); % Initialize as empty table to handle empty data appropriately

% Create a generic set of column names including the ID column
num_questions = 18;
column_names = [{'id'}, arrayfun(@(x) sprintf('Question%d', x), 1:num_questions, 'UniformOutput', false)];

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