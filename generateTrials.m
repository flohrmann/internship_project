%% QUESTION: how exactly to use the eccentricity? 
% can i also use coordinates? this wouldnt change with number of stims tho

function trial_data = generateTrials(n_trials, conditions, n_rows, n_columns, grid_visual_angle, ec_circle, ec_min)
% Initialize empty table to store trial data
trial_data = table();

corner = 2; % min abstand to border of screen
middle = 1; % min abstand to middle of screen

% Loop over the number of trials
for trial = 1:n_trials
    % Determine the target side based on the trial number (so its 50/50 and changes with every trial)
    
    %% QUESTION: is it okay to take the columns right at the middle? or
    % have one column of distance so its very clear if its left/right?


   if mod(trial, 2) == 1
        % Odd trials go to the left side ('L')
        unique_point_col = randi([corner, n_columns/2 - middle]);
        target_side = 'L';
    else
        % Even trials go to the right side ('R')
        unique_point_col = randi([n_columns/2 + middle, n_columns - corner]);
        target_side = 'R';
    end
    
    % Generate a random position for the target object within the determined side of the grid
    number_list = generateNumberList(corner, middle, n_rows); % neither at cornder nor middle
    unique_point_row = getRandomNumberFromList(number_list); 
    % Initialize a matrix with zeros for the grid
    trial_matrix = zeros(n_rows, n_columns);

    % Set the target object to 1
    trial_matrix(unique_point_col, unique_point_row) = 1;

    % Add trial data to the table
    trial_data.trialMatrix{trial, 1} = trial_matrix;
    trial_data.targetPosition(trial, :) = [unique_point_col, unique_point_row];
    trial_data.targetSide{trial, 1} = target_side;
end
% Rename the variables in the table
trial_data.Properties.VariableNames = {'TrialMatrix', 'TargetPosition', 'TargetSide'};
end

function number_list = generateNumberList(corner, middle, n_rows)
    % Generate the first part of the list from 1 to y-1
    part1 = corner:n_rows-middle;
    % Generate the second part of the list from y+1 to x
    part2 = n_rows+middle:n_rows-corner;
    % Concatenate the two parts to form the final list
    number_list = [part1, part2];
end

function random_number = getRandomNumberFromList(number_list)
    % Get the length of the list
    list_length = length(number_list);
    
    % Generate a random index within the range of the list
    random_index = randi(list_length);
    
    % Select the number corresponding to the random index
    random_number = number_list(random_index);
end
