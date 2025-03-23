function plotRTbyEccentricity(data, screenXpixels, screenYpixels, analysis_folder, unique_conditions)
    
    n_rows = size(data(1).angleMatrix{1,1},2);
    n_columns = size(data(1).angleMatrix{1,1},1);

    % Calculate center row and column
    centerRow = round(n_rows / 2);
    centerCol = round(n_columns / 2);
       
    innerPositions = [ % row x column
        4,4; 5,4; 6,4;
        3,5; 4,5; 5,5; 6,5; 7,5;
        3,8; 4,8; 5,8; 6,8; 7,8; 
        4,9; 5,9; 6,9
    ];

    % Initialize arrays to store RTs for inner and outer circle positions
    RT_inner_eye_all = [];
    RT_outer_eye_all = [];
    RT_inner_button_all = [];
    RT_outer_button_all = [];
   

    % Loop through each trial in trial_results
    for p=1:size(data,2)
        participant_data = data(p);

        % Initialize arrays to store RTs for inner and outer circle positions
%         RT_inner_eye = [];
%         RT_outer_eye = [];
%         RT_inner_button = [];
%         RT_outer_button = [];

        for trial = 1:size(data(p).Condition, 1)
            % Extract the target position (row and column)
            targetRow = data(p).TargetPosition(trial, 2);
            targetCol = data(p).TargetPosition(trial, 1);
            
            % Extract the RT for the current trial
            %rt_eye = participant_data.rt_eye(trial);
            %rt_button = participant_data.rt(trial);
    
            % Check if the target position is one of the inner positions
            if ismember([targetRow, targetCol], innerPositions, 'rows')
                %RT_inner_eye = [RT_inner_eye; rt_eye];
                %RT_inner_button = [RT_inner_button; rt_button];
                data(p).central(trial) = 1;
            else
                %RT_outer_eye = [RT_outer_eye; rt_eye];
                %RT_outer_button = [RT_outer_button; rt_button];
                data(p).central(trial) = 0;
            end
        end
    end

    avg_rt_button = struct();
    for p=1:size(data,2)
        avg_rt_button{1,:} = nanmedian(data(1).rt(strcmp(condition, unique_conditions{1})));



    end
    % Plot the RTs for inner vs. outer circle positions
    figure;
    hold on;
    scatter(ones(size(RT_inner_eye)), RT_inner_eye, 50, 'b', 'filled', 'DisplayName', 'Inner Positions');
    scatter(ones(size(RT_outer)) * 2, RT_outer, 50, 'r', 'filled', 'DisplayName', 'Outer Positions');
    xlim([0 3]);
    xticks([1 2]);
    xticklabels({'Inner Positions', 'Outer Positions'});
    title('Reaction Times by Eccentricity (Inner vs. Outer Positions)');
    xlabel('Eccentricity');
    ylabel('Reaction Time (s)');
    ylim([0 max([RT_inner_eye; RT_outer]) + 0.5]); % Adjust y-limits based on data
    grid on;
    legend show;
    hold off;

    % Save the plot
    saveas(gcf, fullfile(analysis_folder, 'RT_by_eccentricity_positions.png'));
end
