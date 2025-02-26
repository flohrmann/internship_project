function plotHistNaNsBufferedDiamSacc(diam_around_stim_table, analysis_folder, bigtitle, name)
    figure; sgtitle(bigtitle);
    title('Distribution of Missing Data Across Trials');
    subplot(1, 2, 1);
    % Assuming `DataPointsBefore` is a cell array of numeric arrays
    nan_counts_before = cellfun(@(x) sum(isnan(x)), diam_around_stim_table.DataPointsBefore);
    histogram(nan_counts_before, 'FaceColor', [0.8, 0.4, 0.2]);
    %histogram(sum(isnan(diam_around_stim_table.DataPointsBefore{:})), 'FaceColor', [0.8, 0.4, 0.2]);
    xlabel('Number of NaNs');
    ylabel('Number of Trials');
    title('NaN Count Before Target Found');

    subplot(1, 2, 2);
    nan_counts_after = cellfun(@(x) sum(isnan(x)), diam_around_stim_table.DataPointsAfter);
    histogram(nan_counts_after, 'FaceColor', [0.2, 0.6, 0.8]);
    %histogram(sum(isnan(diam_around_stim_table.DataPointsAfter), 2), 'FaceColor', [0.2, 0.6, 0.8]);
    xlabel('Number of NaNs');
    ylabel('Number of Trials');
    title('NaN Count After Target Found');

    saveas(gcf, fullfile(analysis_folder, strcat(name,'.png')));
end