%% Function to process the data received from the GUI
function saveData(data, path, name)
disp('Data entered:');
disp(data);

% Save the table in the folder
fileName = fullfile(path,name);
writetable(data, fileName); 
disp(['Data saved as ' fileName]);
end
