function createQuestionnaire(folder_name, subject_id)
    % Create a figure for the GUI
    fig = uifigure('Name', 'Questionnaire', 'Position', [100, 100, 900, 450]);

    % Define the response options
    response_options = {'Never', 'Rarely', 'Sometimes', 'Often', 'Very Often'};

    % Define the questions for Part A and Part B
    partA_questions = {
        'How often do you have trouble wrapping up the final details of a project, once the challenging parts have been done?', ...
        'How often do you have difficulty getting things in order when you have to do a task that requires organization?', ...
        'How often do you have problems remembering appointments or obligations?', ...
        'When you have a task that requires a lot of thought, how often do you avoid or delay getting started?', ...
        'How often do you fidget or squirm with your hands or feet when you have to sit down for a long time?', ...
        'How often do you feel overly active and compelled to do things, like you were driven by a motor?'
    };

    partB1_questions = {
        'How often do you make careless mistakes when you have to work on a boring or difficult project?', ...
        'How often do you have difficulty keeping your attention when you are doing boring or repetitive work?', ...
        'How often do you have difficulty concentrating on what people say to you, even when they are speaking to you directly?', ...
        'How often do you misplace or have difficulty finding things at home or at work?', ...
        'How often are you distracted by activity or noise around you?', ...
        'How often do you leave your seat in meetings or other situations in which you are expected to remain seated?'
    };

    partB2_questions = {
        'How often do you feel restless or fidgety?', ...
        'How often do you have difficulty unwinding and relaxing when you have time to yourself?', ...
        'How often do you find yourself talking too much when you are in social situations?', ...
        'When you’re in a conversation, how often do you find yourself finishing the sentences of the people you are talking to, before they can finish them themselves?', ...
        'How often do you have difficulty waiting your turn in situations when turn taking is required?', ...
        'How often do you interrupt others when they are busy?'
    };

    % Create a cell array to store the responses
    all_questions = [partA_questions, partB1_questions, partB2_questions];
    responses = cell(1, length(all_questions));
    button_groups = cell(1, length(all_questions));

    % Show the initial instructions screen
    showInstructionsScreen(fig, @continueButtonCallback);

    % Callback function for the "Continue" button
    function continueButtonCallback(~, ~)
        % Clear the figure
        clf(fig);
        % Show Part A
        showPart(fig, partA_questions, response_options, 'Part A', @nextButtonCallbackPartB1, 1);
    end

    % Callback function for the "Next" button in Part A
    function nextButtonCallbackPartB1(~, ~)
        % Store Part A responses
        storeResponses(1:length(partA_questions));
        % Clear the figure
        clf(fig);
        % Show Part B1
        showPart(fig, partB1_questions, response_options, 'Part B1', @nextButtonCallbackPartB2, length(partA_questions) + 1);
    end

    % Callback function for the "Next" button in Part B1
    function nextButtonCallbackPartB2(~, ~)
        % Store Part B1 responses
        storeResponses(length(partA_questions) + (1:length(partB1_questions)));
        % Clear the figure
        clf(fig);
        % Show Part B2
        showPart(fig, partB2_questions, response_options, 'Part B2', @submitButtonCallback, length(partA_questions) + length(partB1_questions) + 1);
    end

    % Callback function for the "Submit" button
    function submitButtonCallback(~, ~)
        % Store Part B2 responses
        storeResponses(length(partA_questions) + length(partB1_questions) + (1:length(partB2_questions)));
        % Save responses to a CSV file
        saveResponsesToCSV(responses, all_questions);
        % Close the figure
        close(fig);
        disp('Questionnaire submitted!');
    end

    function storeResponses(question_indices)
        for i = question_indices
            bg = button_groups{i};
            selectedButton = bg.SelectedObject;
            if isempty(selectedButton)
                responses{i} = '';
            else
                responses{i} = selectedButton.UserData;
            end
        end
    end

    function showInstructionsScreen(fig, buttonCallback)
        % Create a title label
        uilabel(fig, 'Text', 'Instructions', 'Position', [20, 410, 960, 30], 'FontSize', 16);

        % Create instructions label
        uilabel(fig, 'Text', ...
            'Please answer the questions rating yourself on each of the criteria shown using the scale on the right side of the page. As you answer each question, click on the box that best describes how you have felt and conducted yourself over the past 6 months.', ...
            'Position', [20, 270, 800, 150], 'FontSize', 16, 'WordWrap', 'on'); % Increased font size

        % Create the "Continue" button
        uibutton(fig, 'Text', 'Continue', 'Position', [450, 20, 100, 30], 'ButtonPushedFcn', buttonCallback);
    end

    function showPart(fig, questions, response_options, titleText, buttonCallback, start_index)
        % Create a title label
        uilabel(fig, 'Text', titleText, 'Position', [20, 410, 960, 30], 'FontSize', 16);

        % Create response option labels at the top
        for j = 1:length(response_options)
            uilabel(fig, 'Text', response_options{j}, 'Position', [460 + (j-1)*80, 370, 90, 22], 'FontSize', 14); % Adjusted position and increased font size
        end

        % Create labels and radio buttons for each question
        num_questions = length(questions);
        y_position = 340;

        for i = 1:num_questions
            uilabel(fig, 'Text', sprintf('%d. %s', start_index + i - 1, questions{i}), 'Position', [20, (y_position-15), 420, 55], 'FontSize', 15, 'WordWrap', 'on'); % Increased font size

            bg = uibuttongroup(fig, 'Position', [450, y_position, 600, 30], 'BorderType', 'none');
            button_groups{start_index + i - 1} = bg;
            for j = 1:length(response_options)
                rb = uiradiobutton(bg, 'Text', '', 'Position', [30 + (j-1)*80, 5, 90, 22], 'FontSize', 1);  % Reduce font size and adjust position for closer buttons
                rb.UserData = response_options{j};  % Store the response option in UserData
            end

            y_position = y_position - 50;  % Reduce space between questions
        end

        % Determine the button text based on the part
        buttonText = 'Next';
        if strcmp(titleText, 'Part B2')
            buttonText = 'Submit';
        end

        % Create the "Next" or "Submit" button
        uibutton(fig, 'Text', buttonText, 'Position', [450, 20, 100, 30], 'ButtonPushedFcn', buttonCallback);
    end

    function saveResponsesToCSV(responses, questions)
        % Combine questions and responses into a cell array
        data = [questions; responses]';

        % Convert to table and write to CSV
        T = cell2table(data, 'VariableNames', {'Question', 'Response'});
        writetable(T, strcat(folder_name, '\questionnaire_responses_', subject_id, '.csv'));
    end
end
