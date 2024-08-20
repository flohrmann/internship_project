function plotQParticipantAnswers(all_data, color_map, safe, comparison_results_folder)
% Function to plot participant answers by question with unique colors for each participant
%
% Inputs:
%   all_data: Combined table of ADHD and non-ADHD data
%   color_map: Color mapping for groups

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
if safe == 1
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf, fullfile(comparison_results_folder, 'quest_all_answer_participants.png'));
end
