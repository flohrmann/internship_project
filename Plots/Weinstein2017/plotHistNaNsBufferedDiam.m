function plotHistNaNsBufferedDiam(diam_around_stim_table, analysis_folder)
    figure;
    title('Distribution of Missing Data Across Trials');
    subplot(1, 2, 1);
    histogram(sum(isnan(diam_around_stim_table.DataPointsBefore), 2), 'FaceColor', [0.8, 0.4, 0.2]);
    xlabel('Number of NaNs');
    ylabel('Number of Trials');
    title('NaN Count Before Target Found');

    subplot(1, 2, 2);
    histogram(sum(isnan(diam_around_stim_table.DataPointsAfter), 2), 'FaceColor', [0.2, 0.6, 0.8]);
    xlabel('Number of NaNs');
    ylabel('Number of Trials');
    title('NaN Count After Target Found');

    saveas(gcf, fullfile(analysis_folder,'histogram_diam_nans_before_after_stim_found.png'));
end