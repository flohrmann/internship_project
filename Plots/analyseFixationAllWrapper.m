function all_fixations = analyseFixationAllWrapper(data, plot_these, dist_threshold, screenXpixels, screenYpixels, safe, path_to_folders, folders, subfolder_fixation, comparison_results_folder, fixation_duration)
min_duration = round(fixation_duration/16); % min num of consecutive points to count as a fixation
fd_label = strcat(num2str(fixation_duration), 'ms');
all_fixations = struct();

try
    load(strcat(comparison_results_folder, 'all_fixations_', fd_label, '.mat'), 'all_fixations'); % load for all participants
catch % if file cant be loaded, compute for everyone
    for participant=1:size(data,2)
        %analysis_folder = strcat(folders{participant}, '\analysis\'); % save per participant
        analysis_folder = strcat(path_to_folders, folders{participant}, subfolder_fixation,'\'); % save per participant
        fixations = analyseFixation(data(participant).id, plot_these, data(participant), dist_threshold, min_duration, screenXpixels, screenYpixels, safe, analysis_folder,fd_label);
        all_fixations(participant).id        = data(participant).id;
        all_fixations(participant).fixations = fixations;
    end
    save(strcat(comparison_results_folder, 'all_fixations_', fd_label, '.mat'), 'all_fixations'); % save for all participants
end

if size(all_fixations, 2) == size(folders, 1) % if this contains data for all participants do nothing
else % if file doesnt contain everyone, compute for everyone
    for participant=1:size(data,2)
        %analysis_folder = strcat(folders{participant}, '\analysis\'); % save per participant
        analysis_folder = strcat(folders{participant}, subfolder_fixation,'\'); % save per participant
        fixations = analyseFixation(data(participant).id, plot_these, data(participant), dist_threshold, min_duration, screenXpixels, screenYpixels, safe, analysis_folder, fd_label);
        all_fixations(participant).id        = data(participant).id;
        all_fixations(participant).fixations = fixations;
    end
    save(strcat(comparison_results_folder, 'all_fixations_', fd_label', '.mat'), 'all_fixations'); % save for all participants
end