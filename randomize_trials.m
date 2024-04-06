function rand_trials = randomize_trials(data, folder_name)
% Generate a random permutation of indices
num_rows = height(data);
random_indices = randperm(num_rows);

% Rearrange the rows of the table using the random indices
rand_trials = data(random_indices, :);
rand_data_file_name = fullfile(folder_name, 'rand_trials.mat');
save(rand_data_file_name, 'rand_trials');