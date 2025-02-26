function saccades = reverseEngineerSaccades(cut_data, fixations, screenXpixels, screenYpixels, plot_these, analysis_folder, safe)

for trial = 1:size(cut_data.stimulusTrial, 1)
    saccades(trial).trial = trial;
    current_fixations = fixations.stimulusFixations(trial).fixations; %.stimulusFixations;
    current_data = cut_data.stimulusTrial(trial);
    stim_start = fixations.stimulusFixations(trial).stimStart;
    saccades(trial).saccades = detectSaccadesFromFixations(current_fixations, stim_start, screenXpixels, screenYpixels);
    saccades(trial).saccadeDurations = getSaccadeTimes(current_fixations, current_data); % double check times -> looks good!
end
%saccades.length_fixation = {fixations.length};
save(fullfile(analysis_folder, '\all_saccades.mat'), 'saccades');% save the saccades data

% plot saccades
if length(plot_these) >= 1
    plot_these = [5, 19];
    plotFixationSaccades(cut_data, saccades, fixations, screenXpixels, screenYpixels, plot_these, analysis_folder, safe);
    close all
else
end