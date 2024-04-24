% Start Tobii Pro SDK Operations
tobii = EyeTrackingOperations();

% Find connected eye trackers
found_eyetrackers = tobii.find_all_eyetrackers();
if isempty(found_eyetrackers)
    error('No eye trackers found. Make sure the eye tracker is connected and try again.');
end
my_eyetracker = found_eyetrackers(1);
eyetracker = my_eyetracker;

% Display eye tracker details
disp(["Address: ", my_eyetracker.Address]);
disp(["Model: ", my_eyetracker.Model]);
disp(["Name (It's OK if this is empty): ", my_eyetracker.Name]);
disp(["Serial number: ", my_eyetracker.SerialNumber]);

%% Calibration


eyetracker_address = my_eyetracker.Address;


% Check if a specific eye tracker address has been provided, and if so,
% try to locate it and return the corresponding eye tracker object.
if (~strcmp(eyetracker_address, ''))
    try
        eyetracker = Tobii.get_eyetracker(eyetracker_address);
    catch ME
        if (strcmp(ME.identifier,'EyeTrackerGet:error204'))
            fprintf('Unable to connect to specified eye tracker.\n');
            return
        end
    end
else
    try
        eyetrackers = Tobii.find_all_eyetrackers();
        eyetracker = eyetrackers(1);
    catch
        fprintf("Unable to find any connected eye trackers.\n");
        return
    end
end

eyetracker.get_gaze_data();
% check if gaze data collection seems to be working
% TODO before commit: pause
pause(1);
result = eyetracker.get_gaze_data();

if isa(result,'StreamError')
    fprintf('Gaze data stream Error: %s\n',string(result.Error.value));
    fprintf('Source: %s\n',string(result.Source.value));
    fprintf('SystemTimeStamp: %d\n',result.SystemTimeStamp);
    fprintf('Message: %s\n',result.Message);
    return
end

% Start calibration procedure (while gaze recording is running
% in the background)
calib = ScreenBasedCalibration(eyetracker);

try
    calib.enter_calibration_mode();
catch ME
    if (strcmp(ME.identifier,'EnterCalibrationMode:error210'))
        fprintf('The previous calibration was not completed!\n');
        calib.leave_calibration_mode()
        fprintf('Calibration is restarted\n');
        calib.enter_calibration_mode()
    else
        fprintf('%s\n', ME.identifier);
        return
    end
end
% Define the points on screen we should calibrate at.
% The coordinates are normalized, i.e. [0.0, 0.0] is the upper left corner and
% [1.0, 1.0] is the lower right corner.
points_to_collect = [[0.1,0.1];[0.1,0.9];[0.5,0.5];[0.9,0.1];[0.9,0.9]];

% When collecting data a point should be presented on the screen in the
% appropriate position. We do not provide this here, but you may wish to
% use Matlab's own 'imshow' function for this, or a 3rd party package
% like Psychtoolbox for this purpose.
for i=1:size(points_to_collect,1)
    collect_result = calib.collect_data(points_to_collect(i,:));
    fprintf('Point [%.2f,%.2f] Collect Result: %d\n',points_to_collect(i,:),collect_result.value);
end

calibration_result = calib.compute_and_apply();
fprintf('Calibration Status: %d\n',calibration_result.Status.value.value);

% After analysing the calibration result one might want to re-calibrate
% some of the points (this would need to be done in an interactive fashion
% however - you would most likely want to write your own logic for handling
% this, or find a 3rd party package that 'wraps' this Tobii package and handles
% things for you).
points_to_recalibrate = [[0.1,0.1];[0.1,0.9]];

for i=1:size(points_to_recalibrate,1)
    calib.discard_data(points_to_recalibrate(i,:));
    collect_result = calib.collect_data(points_to_recalibrate(i,:));
    fprintf('Point [%.2f,%.2f] Collect Result: %d\n',points_to_recalibrate(i,:),collect_result.value);
end

calibration_result = calib.compute_and_apply();
fprintf('Calibration Status: %d\n',calibration_result.Status.value.value);

% Check if the calibration result was successful, meaning a calibration could
% be performed. Note that depending on your requirements, the calibration
% might still not be of high enough quality. You would need to also run a
% validation of the calibration to evaluate this.
if calibration_result.Status == CalibrationStatus.Success
    points = calibration_result.CalibrationPoints;

    number_points = size(points,2);
    % plot the calibration points and associated data
    for i=1:number_points
        plot(points(i).PositionOnDisplayArea(1),points(i).PositionOnDisplayArea(2),'ok','LineWidth',10);
        mapping_size = size(points(i).RightEye,2);
        set(gca, 'YDir', 'reverse');
        axis([-0.2 1.2 -0.2 1.2])
        hold on;
        for j=1:mapping_size
            if points(i).LeftEye(j).Validity == CalibrationEyeValidity.ValidAndUsed
                plot(points(i).LeftEye(j).PositionOnDisplayArea(1), points(i).LeftEye(j).PositionOnDisplayArea(2),'-xr','LineWidth',3);
            end
            if points(i).RightEye(j).Validity == CalibrationEyeValidity.ValidAndUsed
                plot(points(i).RightEye(j).PositionOnDisplayArea(1),points(i).RightEye(j).PositionOnDisplayArea(2),'xb','LineWidth',3);
            end
        end

    end
end

calib.leave_calibration_mode();

% Now that calibration is finished, get the gaze data recorded in the meantime.
gaze_data = eyetracker.get_gaze_data();

eyetracker.stop_gaze_data();

n_collected = size(gaze_data, 1);
fprintf('Collected %d data points\n', n_collected);

% Store the relevant data to N-by-1 arrays to simplify saving below
% (please see the reference guide documentation on GazeData for
% information on what other properties/data are available)
left_eye_xs = zeros(n_collected, 1);
right_eye_ys = zeros(n_collected, 1);
left_eye_xs = zeros(n_collected, 1);
right_eye_ys = zeros(n_collected, 1);
sample_times_from_start = zeros(n_collected, 1);

first_t = gaze_data(1).SystemTimeStamp;
for i=1:length(gaze_data)
    single_data = gaze_data(i);
    left_gp = single_data.LeftEye.GazePoint;
    right_gp = single_data.RightEye.GazePoint;
    left_eye_xs(i) = left_gp.OnDisplayArea(1);
    left_eye_ys(i) = left_gp.OnDisplayArea(2);
    right_eye_xs(i) = right_gp.OnDisplayArea(1);
    right_eye_ys(i) = right_gp.OnDisplayArea(2);
    sample_time_from_start_us = single_data.SystemTimeStamp - first_t;
    % convert from microseconds to seconds
    sample_time_from_start_s = double(sample_time_from_start_us) / 1000000;
    sample_times_from_start(i) = sample_time_from_start_s;
end

output_matrix = [
    sample_times_from_start(:), ...
    left_eye_xs(:), ...
    left_eye_ys(:), ...
    right_eye_xs(:), ...
    right_eye_ys(:)
];

% Save the data to a CSV format file
% (first writing column headers to file)
header_lines = ["time_seconds", "left_x", "left_y", "right_x", "right_y"];
writematrix(header_lines, 'gaze_data.csv');
writematrix(output_matrix, 'gaze_data.csv', 'WriteMode', 'append');



%% Start collecting gaze data
% my_eyetracker.get_gaze_data();
% 
% pause(1); % Collect data for a second
% 
% % Attempt to retrieve and display the latest gaze data
% gaze_data = my_eyetracker.get_gaze_data();
% if isempty(gaze_data)
%     error('No gaze data retrieved. Check the eye tracker connection and settings.');
% end
% latest_gaze_data = gaze_data(end);
% 
% % Display gaze data for both eyes
% disp('SystemRequestTimeStamp: ' + string(latest_gaze_data.SystemTimeStamp));
% disp('DeviceTimeStamp: ' + string(latest_gaze_data.DeviceTimeStamp));
% 
% display_eye_data('Left Eye', latest_gaze_data.LeftEye);
% display_eye_data('Right Eye', latest_gaze_data.RightEye);
% 
% % Stop collecting gaze data
% my_eyetracker.stop_gaze_data();
% 
% % Helper function to display eye data
% function display_eye_data(eye_label, eye_data)
%     fprintf('%s\n', eye_label);
%     fprintf('GazePoint.OnDisplayArea: %.2f %.2f\n', eye_data.GazePoint.OnDisplayArea);
%     fprintf('GazePoint.InUserCoordinateSystem: %.2f %.2f %.2f\n', eye_data.GazePoint.InUserCoordinateSystem);
%     %fprintf('GazePoint.Validity: %s\n', char(eye_data.GazePoint.Validity.Valid));
%     fprintf('GazeOrigin.InUserCoordinateSystem: %.2f %.2f %.2f\n', eye_data.GazeOrigin.InUserCoordinateSystem);
%     fprintf('GazeOrigin.InTrackBoxCoordinateSystem: %.2f %.2f %.2f\n', eye_data.GazeOrigin.InTrackBoxCoordinateSystem);
%     %fprintf('GazeOrigin.Validity: %s\n', char(eye_data.GazeOrigin.Validity.Valid));
%     fprintf('Pupil.Diameter: %.2f\n', eye_data.Pupil.Diameter);
%     %fprintf('Pupil.Validity: %s\n', char(eye_data.Pupil.Validity));
% end

