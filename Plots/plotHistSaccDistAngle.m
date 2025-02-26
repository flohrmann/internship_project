function plotHistSaccDistAngle(distances_1, angles_1, name_1, ...
                               distances_2, angles_2, name_2,...
                               bigname)
                           
figure;
t = tiledlayout(2, 2, 'TileSpacing', 'Compact');

% Histogram for distances: Actual
nexttile;
histogram(distances_1, 'FaceColor', 'r', 'BinWidth', 50);
title(name_1);
xlabel('Distance');
ylabel('Frequency');

% Histogram for distances: Optimal
nexttile;
histogram(distances_2, 'FaceColor', 'b', 'BinWidth', 50);
title(name_2);
xlabel('Distance');
ylabel('Frequency');

% Histogram for angles: Actual
nexttile;
histogram(angles_1, 'FaceColor', 'r', 'BinWidth', 10);
title(name_1);
xlabel('Angle (degrees)');
ylabel('Frequency');

% Histogram for angles: Optimal
nexttile;
histogram(angles_2, 'FaceColor', 'b', 'BinWidth', 10);
title(name_2);
xlabel('Angle (degrees)');
ylabel('Frequency');

% Add a shared title for the layout
title(t,bigname);