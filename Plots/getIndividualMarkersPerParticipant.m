function participant_map = getIndividualMarkersPerParticipant(group_labels, ids, comparison_results_folder)
    % Input:
    %   group_labels: Cell array of group labels ('nonADHD' or 'ADHD')
    %   ids: Array of participant IDs (unique keys for the map)
    % Output:
    %   participant_map: Map with IDs as keys and a struct of color and marker as values

    % Total number of participants
    numParticipants = length(group_labels);

    % Generate a rainbow gradient colormap with space between groups
    rainbow_colors = jet(numParticipants + 14); % Add extra colors for spacing
    numNonADHD = sum(strcmp(group_labels, 'nonADHD'));
    numADHD = sum(strcmp(group_labels, 'ADHD'));
    
    
    % Assign nonADHD and ADHD colors with a gap in between
    adhd_colors = rainbow_colors(10:numADHD+10, :);                 % blue to green
    nonadhd_colors = rainbow_colors(numParticipants-numNonADHD+14:end, :);   % yellow to red, skip some in between for spacing between groups

    % Define unique markers for participants
    markers = {'o', 's', '^','d', 'v', 'p', '<', 'h','>', '*'};
    numMarkers = length(markers);

    % Initialize the map
    participant_map = containers.Map('KeyType', 'double', 'ValueType', 'any');

    % Assign colors and markers for nonADHD participants
    nonadhd_index = 1;
    for i = 1:numParticipants
        if strcmp(group_labels{i}, 'nonADHD')
            id = ids(i); % Use participant ID as the map key
            participant_map(id) = struct( ...
                'color', nonadhd_colors(nonadhd_index, :), ...
                'marker', markers{mod(nonadhd_index-1, numMarkers) + 1});
            nonadhd_index = nonadhd_index + 1;
        end
    end

    % Assign colors and markers for ADHD participants
    adhd_index = 1;
    for i = 1:numParticipants
        if strcmp(group_labels{i}, 'ADHD')
            id = ids(i); % Use participant ID as the map key
            participant_map(id) = struct( ...
                'color', adhd_colors(adhd_index, :), ...
                'marker', markers{mod(adhd_index-1, numMarkers) + 1});
            adhd_index = adhd_index + 1;
        end
    end





% Example data for plotting
numObservations = 2;
data = rand(numObservations, length(ids)); % Simulated data for participants

% Plot
figure; hold on;
for i = 1:length(ids)
    participant_id = ids(i);
    info = participant_map(participant_id); % Retrieve color and marker
    plot(data(:, i), ...
        'Color', info.color, ...
        'Marker', info.marker, ...
        'LineStyle', '-', ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('Participant %d (%s)', participant_id, group_labels{i}));
end

% Add legend and labels
legend('show', 'Location', 'best');
xlabel('Observation');
ylabel('Value');
title('Participant Fake Data to test Colors and Markers');
hold off;

saveas(gcf, fullfile(comparison_results_folder, '00_individual_markers_test.png'));

