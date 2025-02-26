function all_saccades = analyseSaccadesAllWrapper(data, all_fixations, plot_these, screenXpixels, screenYpixels, safe, folders, subfolder_fixation, comparison_results_folder, fixation_duration)

fd_label = strcat(num2str(fixation_duration), 'ms');
all_saccades = struct();
try
    load(fullfile(comparison_results_folder, 'all_saccades_', fd_label, 'mat'), 'all_saccades');% load the saccades data
catch
    for participant = 1:size(data, 2)
        saccades = struct();
        for trial = 1:size(data(participant).stimulusTrial, 1)
            saccades(trial).trial = trial;
            participant_data = data(participant).stimulusTrial(trial);
            participant_fixations = all_fixations(participant).fixations.stimulusFixations(trial).fixations; %.stimulusFixations;
            stim_start = all_fixations(participant).fixations.stimulusFixations(trial).stimStart;
            saccades(trial).saccades = detectSaccadesFromFixations(participant_fixations, stim_start, screenXpixels, screenYpixels);
            saccades(trial).saccadeDurations = getSaccadeTimes(participant_fixations, participant_data); % double check times -> looks good!
        end
        all_saccades(participant).id = data(participant).id;
        all_saccades(participant).saccades = saccades;
    end
    save(strcat(comparison_results_folder, 'all_saccades_', fd_label, 'mat'), 'all_saccades');% save the saccades data
end

if size(all_saccades, 2) == size(folders, 1) % if this contains data for all participants do nothing
else % if file doesnt contain everyone, compute for everyone
    for participant = 1:size(data, 2)
        saccades = struct();
        for trial = 1:size(data(participant).stimulusTrial, 1)
            saccades(trial).trial = trial;
            participant_data = data(participant).stimulusTrial(trial);
            participant_fixations = all_fixations(participant).fixations.stimulusFixations(trial).fixations; %.stimulusFixations;
            stim_start = all_fixations(participant).fixations.stimulusFixations(trial).stimStart;
            saccades(trial).saccades = detectSaccadesFromFixations(participant_fixations, stim_start, screenXpixels, screenYpixels);
            saccades(trial).saccadeDurations = getSaccadeTimes(participant_fixations, participant_data); % double check times -> looks good!
        end
        all_saccades(participant).id = data(participant).id;
        all_saccades(participant).saccades = saccades;
    end
    save(fullfile(comparison_results_folder, 'all_saccades_', fd_label, 'mat'), 'all_saccades');% save the saccades data
end

% plot saccades
for participant = 1:size(data, 2)
    analysis_folder = strcat(folders{participant}, subfolder_fixation, '\');
    participant_data = data(participant).stimulusTrial;
    participant_fixations = all_fixations(participant).fixations.stimulusFixations;
    plotFixationSaccades(data(participant), all_saccades(participant), all_fixations(participant), screenXpixels, screenYpixels, plot_these, analysis_folder, safe);
end
close all