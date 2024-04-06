function [data] = startGUI()
% Function to open GUI and get user input
% Create a figure for the GUI
fig = figure('Name', 'Data Input', 'NumberTitle', 'off', 'Position', [120, 120, 400, 300]);

% Create text boxes for the user to enter today's date, ID, name, trial length, and notes
uicontrol('Style', 'text', 'String', 'Today''s Date:', 'Position', [20, 240, 100, 20]);
dateEditBox = uicontrol('Style', 'edit', 'String', datestr(now, 'yyyy-mm-dd'), 'Position', [160, 240, 150, 20]);

uicontrol('Style', 'text', 'String', 'ID:', 'Position', [20, 210, 100, 20]);
idEditBox = uicontrol('Style', 'edit', 'String', '', 'Position', [160, 210, 150, 20]);

uicontrol('Style', 'text', 'String', 'Name:', 'Position', [20, 180, 100, 20]);
nameEditBox = uicontrol('Style', 'edit', 'String', '', 'Position', [160, 180, 150, 20]);

uicontrol('Style', 'text', 'String', 'Number Trials:', 'Position', [20, 150, 100, 20]);
trialLengthEditBox = uicontrol('Style', 'edit', 'String', '120', 'Position', [160, 150, 150, 20]);

uicontrol('Style', 'text', 'String', 'Conditions:', 'Position', [20, 120, 100, 20]);
conditionsEditBox = uicontrol('Style', 'edit', 'String', 'a_simple,b_simple,a,b', 'Position', [160, 120, 150, 20]);

uicontrol('Style', 'text', 'String', 'Notes:', 'Position', [20, 90, 100, 20]);
notesEditBox = uicontrol('Style', 'edit', 'String', '', 'Position', [160, 90, 150, 20]);

% Create a button for the user to submit the data
submitButton = uicontrol('Style', 'pushbutton', 'String', 'Lets go', 'Position', [150, 30, 100, 30], 'Callback', @submitCallback);

% Wait for the user to close the figure
waitfor(fig);

%% custom callback
    function submitCallback(~, ~)
        % Get the data from the edit boxes
        data.date = get(dateEditBox, 'String');
        data.id = get(idEditBox, 'String');
        data.name = get(nameEditBox, 'String');
        data.n_trials = str2double(get(trialLengthEditBox, 'String'));
        data.conditions = get(conditionsEditBox, 'String');
        data.notes = get(notesEditBox, 'String');

        % Close the figure
        close(fig);
    end
end