function saccadeStats = analyzeSaccadeAmplitudeAndVelocity(all_saccades, data, conditions)
    saccadeStats = struct();
    
    % Initialize conditions
    numConditions = length(conditions);
    
    for participant = 1:length(all_saccades)
        participant_saccades = all_saccades(participant).saccades;
        conditionAmplitudes = cell(numConditions, 1);
        conditionVelocities = cell(numConditions, 1);
        conditionSaccadeCounts = zeros(numConditions, 1);
        conditionTotalAmplitudes = zeros(numConditions, 1);
        conditionTotalDurations = zeros(numConditions, 1);
        
        % Loop over each trial
        for trial = 1:length(participant_saccades)
            saccades = participant_saccades(trial).saccades;
            saccadeCount = length(saccades);
            condition = data(participant).Condition{trial};
            conditionIdx = find(strcmp(conditions, condition));

            % Sum saccade distances and durations
            for sIdx = 1:saccadeCount
                saccadeDistance = saccades(sIdx).distance;
                saccadeDuration = saccades(sIdx).duration;
                
                % Accumulate statistics per condition
                conditionTotalAmplitudes(conditionIdx) = conditionTotalAmplitudes(conditionIdx) + saccadeDistance;
                conditionTotalDurations(conditionIdx) = conditionTotalDurations(conditionIdx) + saccadeDuration;
                
                % Store individual saccade amplitudes and velocities
                conditionAmplitudes{conditionIdx} = [conditionAmplitudes{conditionIdx}, saccadeDistance];
                conditionVelocities{conditionIdx} = [conditionVelocities{conditionIdx}, saccadeDistance / saccadeDuration];
            end
            
            % Track the number of saccades per condition
            conditionSaccadeCounts(conditionIdx) = conditionSaccadeCounts(conditionIdx) + saccadeCount;
        end
        
        % Calculate averages per condition
        for c = 1:numConditions
            if conditionSaccadeCounts(c) > 0
                avgSaccadeAmplitude = conditionTotalAmplitudes(c) / conditionSaccadeCounts(c);
                avgSaccadeVelocity = conditionTotalAmplitudes(c) / conditionTotalDurations(c);
                avgSaccadeCount = conditionSaccadeCounts(c) / (length(participant_saccades)/4); % Average number of saccades per trial
            else
                avgSaccadeAmplitude = NaN;
                avgSaccadeVelocity = NaN;
                avgSaccadeCount = NaN;
            end
            
            % Store results for this condition
            saccadeStats(participant).id = all_saccades(participant).id;
            saccadeStats(participant).conditions(c).name = conditions{c};
            saccadeStats(participant).conditions(c).avgSaccadeAmplitude = avgSaccadeAmplitude;
            saccadeStats(participant).conditions(c).saccadeAmplitudes = conditionAmplitudes{c};
            saccadeStats(participant).conditions(c).avgSaccadeVelocity = avgSaccadeVelocity;
            saccadeStats(participant).conditions(c).saccadeVelocities = conditionVelocities{c};
            saccadeStats(participant).conditions(c).saccadeCount = conditionSaccadeCounts(c);
            saccadeStats(participant).conditions(c).avgSaccadeCount = avgSaccadeCount;
        end
    end
end
