function results_table = questionnaireScale(data, folder)
    %% Adult ADHD Self-Report Scale (ASRS-v1.1) Symptom Checklist
    % "If four or more marks appear in the darkly shaded boxes within Part A then the patient has 
    % symptoms highly consistent with ADHD in adults and further investigation is warranted."
    % "No total score or diagnostic likelihood is utilized for the twelve questions. It has been 
    % found that the six questions in Part A are the most predictive of the disorder and are best 
    % for use as a screening instrument."

    % Extract IDs and scores
    ids = data{:, 1}; % First column with participant IDs
    scores = data{:, 2:19}; % Columns 2-19 contain the test scores
    a_scores = data{:, 2:7}; % Questionnaire Part A

    % Define thresholds for each question (grey parts from ASRS-v1.1)
    thresholds = [3, 3, 3, 4, 4, 4, ... % Part A thresholds
                  4, 4, 3, 4, 4, 3, 4, 4, 4, 3, 4, 3]; % Part B thresholds

    % Apply thresholds to determine if each score counts towards a symptom
    symptom_counts = scores >= thresholds;
    a_symptom_counts = a_scores >= thresholds(1:6);
    
    % Mean answers 
    mean_symptoms = sum(scores, 2);
    a_mean_symptoms = sum(a_scores, 2);

    % Sum the number of symptoms per participant
    sum_symptom_counts = sum(symptom_counts, 2);
    sum_a_symptom_counts = sum(a_symptom_counts, 2);

    results_table = table(ids, sum_symptom_counts, mean_symptoms, sum_a_symptom_counts, a_mean_symptoms, ...
                    'VariableNames', {'id', 'TotalSymptomsCount', 'MeanSymptoms', 'PartASymptomsCount', 'PartAMeanSymptoms'});
    disp(results_table); 
    
    % Save the results as a .mat file in the specified folder
    mat_filename = fullfile(folder, 'questionnaire_symptom_counts.mat');
    save(mat_filename, 'results_table');
