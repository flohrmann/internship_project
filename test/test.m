%% safe seed
rand('state', seed4);
random_number4 = rand;

%% Get Eye Image -> not supported
clear;
Tobii = EyeTrackingOperations();
eyetrackers = Tobii.find_all_eyetrackers();
eyetracker_address = eyetrackers.Address;
% Example:
% eyetracker_address = 'tet-tcp://10.46.32.51';
try
    eyetracker = Tobii.get_eyetracker(eyetracker_address);
catch ME
    if (strcmp(ME.identifier,'EyeTrackerGet:error204'))
        fprintf('Unable to connect eye tracker.\n');
        return
    end
end
% It is possible to check if the eyetracker supports this stream
if ismember(EyeTrackerCapabilities.HasEyeImages,eyetracker.DeviceCapabilities)
    disp('Eye Image Supported');
else
    disp('Eye Image Not Supported');
end
% The first call subscribes to the stream and returns either data
% (might be empty if no data has been received yet) or any error that
% happened during the subscription.
eyetracker.get_eye_image();
pause(1);
result = eyetracker.get_eye_image();
if isa(result,'StreamError')
    fprintf('Error: %s\n',string(result.Error.value));
    fprintf('Source: %s\n',string(result.Source.value));
    fprintf('SystemTimeStamp: %d\n',result.SystemTimeStamp);
    fprintf('Message: %s\n',result.Message);
elseif isa(result,'EyeImage')
    % Collect data for 1 seconds.
    pause(1);
    % The subsequent calls return the current values in the stream buffer.
    % If a flat structure is prefered just use an extra input 'flat'.
    % i.e. eye_images = eyetracker.get_eye_image('flat');
    eye_images = eyetracker.get_eye_image();
    eyetracker.stop_eye_image();
    fprintf('Collected %d eye images\n',size(eye_images,2));
    % To select the most recent data point simply select the last value of the
    % buffer.
    latest_eye_image = eye_images(end);
    fprintf('SystemTimeStamp: %d\n',latest_eye_image.SystemTimeStamp);
    fprintf('DeviceTimeStamp: %d\n',latest_eye_image.DeviceTimeStamp);
    fprintf('CameraId: %d\n',latest_eye_image.CameraId);
    fprintf('Type: %d\n',latest_eye_image.Type.value);
    imshow(latest_eye_image.Image);
end
tracker.deInit();
