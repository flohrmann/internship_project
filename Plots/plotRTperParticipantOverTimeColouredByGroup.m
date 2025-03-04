function plotRTperParticipantOverTimeColouredByGroup(data, color_map, safe, comparison_results_folder)

% rt over time per participant with different colors for ADHD and non-ADHD
figure;
hold on;

% Loop through each participant and plot their reaction times
for i = 1:length(data)
    % Determine color based on group
    if strcmp(data(i).group, 'ADHD')
        plot_color = color_map('ADHD');
    else
        plot_color = color_map('nonADHD');
    end
    
    % Plot reaction times for this participant
    plot(1:length(data(i).rt), data(i).rt, 'Color', plot_color, 'DisplayName', ['ID ' num2str(data(i).id)]);
end

title('Reaction Times Over Trials');
xlabel('Trial Number');
ylabel('Reaction Time');

% Update legend to show ADHD and non-ADHD groups
legend('Location', 'northeast');
hold off;

if safe == 1
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf, fullfile(comparison_results_folder, 'lines_RTbutton_allTrials_perParticipant.png'));
else
end