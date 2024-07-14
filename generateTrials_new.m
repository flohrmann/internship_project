function trial_data = generateTrials_new(n_trials, n_rows, n_columns, grid_visual_angle, ec_circle, ec_min)
% Generate matrix: 23232323   where 1 is for target bar/stimulus
%                  32321232         2 is for 1st distractor
%                  23232323         3 is for 2nd distractor
% Position of 1 is generated randomly and annotated with its position
% (left/right side of screen)

% TODO implement eccentricity

% Initialize empty table to store trial data
trial_data = table();

corner = 1; % min abstand to border of screen
middle = 1; % min abstand to middle of screen

% get all possible column/row indices
possible_rows = 1+corner:n_rows-corner;
num_rows = length(possible_rows);

possible_columns = removeMiddleValues(n_columns);
num_col_idx = length(possible_columns);
possible_columns = possible_columns(1+corner:num_col_idx-corner);
% Split columns into left and right halves
%half_idx = floor(num_col_idx / 2);
left_half = possible_columns(1:(length(possible_columns)/2));
right_half = possible_columns((length(possible_columns)/2 +1):end);

% Generate a random position for the target object within the determined side of the grid
for trial = 1:n_trials
    % Determine the target side based on the trial number (so its 50/50 and changes with every trial)
    if mod(trial, 2) == 1
        % Odd trials go to the left side ('L')
        unique_point_col = left_half(randi(length(left_half)));
        target_side = 'L';
    else
        % Even trials go to the right side ('R')
        unique_point_col = right_half(randi(length(right_half)));
        target_side = 'R';
    end
    
    unique_point_row = possible_rows(randi(num_rows));
    %unique_point_row = getRandomNumberFromList(possible_rows);
    % Initialize a matrix with zeros for the grid
    %trial_matrix = zeros(n_rows, n_columns);
    
    % Initialize a matrix with 2s and 3s for the grid/distractors
    trial_matrix = 2 + mod((1:n_rows)' + (1:n_columns), 2);  % Adds 1 wherever the sum of indices is odd
    
    % Set the target object to 1
    trial_matrix(unique_point_row, unique_point_col) = 1;
    
    % Add trial data to the table
    trial_data(trial, :) = {trial_matrix, [unique_point_col, unique_point_row], target_side};
    
end
% Rename the variables in the table
trial_data.Properties.VariableNames = {'TrialMatrix', 'TargetPosition', 'TargetSide'};
end

%%
function resultList = removeMiddleValues(n_col)
% Create a list of column indices
colIndices = 1:n_col;

% Check if the number of columns is even
if mod(n_col, 2) == 0
    % If even, remove the two middle values
    middle1 = n_col / 2;
    middle2 = middle1 + 1;
    resultList = colIndices([1:middle1-1, middle2+1:end]);
else
    % If odd, remove the middle value
    middle = (n_col + 1) / 2;
    resultList = colIndices([1:middle-1, middle+1:end]);
end
end

%%
function number_list = generateNumberList(corner, middle, n_rows)
% Generate the first part of the list from 1 to y-1
part1 = 1+corner:ceil(n_rows/2-middle);
% Generate the second part of the list from y+1 to x
part2 = ceil(n_rows/2+middle):n_rows-corner;
% Concatenate the two parts to form the final list
number_list = [part1, part2];
end

%%
function random_number = getRandomNumberFromList(number_list)
% Get the length of the list
list_length = length(number_list);
% Generate a random index within the range of the list
random_index = randi(list_length);
% Select the number corresponding to the random index
random_number = number_list(random_index);
end
