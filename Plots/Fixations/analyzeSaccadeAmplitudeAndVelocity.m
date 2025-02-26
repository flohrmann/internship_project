function saccadeStats = analyzeSaccadeAmplitudeAndVelocity(all_saccades, data, conditions)
saccadeStats = struct();

% Initialize conditions
numConditions = length(conditions);

for participant = 1:length(all_saccades)
    participant_saccades = all_saccades(participant).saccades;
    trialDistances = zeros(length(participant_saccades), 1); % Total distance per trial
    trialVelocities = cell(length(participant_saccades), 1);
    trialDistancesPerSecond = nan(length(participant_saccades), 1); % Distance per second
    conditionAmplitudes = cell(numConditions, 1);
    conditionVelocities = cell(numConditions, 1);
    conditionDistancesPerSecond = cell(numConditions, 1);
    
    % Loop over each trial
    for trial = 1:length(participant_saccades)
        saccades = participant_saccades(trial).saccades;
        saccadeCount = length(saccades);
        condition = data(participant).Condition{trial};
        conditionIdx = find(strcmp(conditions, condition));
        rt = data(participant).rt(trial); % Reaction time for this trial
        
        totalDistance = 0;
        trialVelocitiesTmp = [];
        
        % Process saccades in this trial
        for sIdx = 1:saccadeCount
            saccadeDistance = saccades(sIdx).distance;
            saccadeDuration = saccades(sIdx).duration;
            
            totalDistance = totalDistance + saccadeDistance; % Sum distances
            trialVelocitiesTmp = [trialVelocitiesTmp, saccadeDistance / saccadeDuration];
        end
        
        trialDistances(trial) = totalDistance; % Store total distance for this trial
        
        % Calculate distance per second (avoid division by zero)
        if rt > 0
            trialDistancesPerSecond(trial) = totalDistance / rt;
        else
            trialDistancesPerSecond(trial) = NaN; % Handle invalid RTs
        end
        
        % Store trial-level velocities
        if ~isempty(trialVelocitiesTmp)
            trialVelocities{trial} = nanmedian(trialVelocitiesTmp);
        else
            trialVelocities{trial} = NaN;
        end
        
        % Aggregate by condition
        conditionAmplitudes{conditionIdx} = [conditionAmplitudes{conditionIdx}, totalDistance];
        conditionVelocities{conditionIdx} = [conditionVelocities{conditionIdx}, trialVelocities{trial}];
        conditionDistancesPerSecond{conditionIdx} = [conditionDistancesPerSecond{conditionIdx}, trialDistancesPerSecond(trial)];
    end
    
    % Calculate medians per condition (aggregating over trials)
    for c = 1:numConditions
        if ~isempty(conditionAmplitudes{c})
            medianSaccadeAmplitude = nanmedian(conditionAmplitudes{c});
            medianSaccadeVelocity = nanmedian(conditionVelocities{c});
            medianDistancePerSecond = nanmedian(conditionDistancesPerSecond{c});
        else
            medianSaccadeAmplitude = NaN;
            medianSaccadeVelocity = NaN;
            medianDistancePerSecond = NaN;
        end
        
        % Store results for this condition
        saccadeStats(participant).id = all_saccades(participant).id;
        saccadeStats(participant).conditions(c).name = conditions{c};
        saccadeStats(participant).conditions(c).medianSaccadeAmplitude = medianSaccadeAmplitude;
        saccadeStats(participant).conditions(c).saccadeAmplitudes = conditionAmplitudes{c};
        saccadeStats(participant).conditions(c).medianSaccadeVelocity = medianSaccadeVelocity;
        saccadeStats(participant).conditions(c).saccadeVelocities = conditionVelocities{c};
        saccadeStats(participant).conditions(c).medianDistancePerSecond = medianDistancePerSecond;
    end
    
    % Store trial-level results
    saccadeStats(participant).trials = struct('totalDistances', trialDistances, ...
                                              'velocities', trialVelocities, ...
                                              'distancesPerSecond', trialDistancesPerSecond);
end
end
