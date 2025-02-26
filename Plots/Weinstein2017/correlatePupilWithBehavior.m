% Correlate pupil dilation with behavioral measures
function [correlation, p_value] = correlatePupilWithBehavior(pupil_dilation, behavior_metric)
    [correlation, p_value] = corr(pupil_dilation, behavior_metric, 'Type', 'Spearman');
end
