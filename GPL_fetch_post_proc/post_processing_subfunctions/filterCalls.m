function [figCalls] = filterCalls(figPARAMS,plot_calls,STATS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will filter out potential calls by the standard devation
% parameters.
% Written: Ian 09/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over all metrics
for n = 1:length(figPARAMS.filter.metric)

    % Unload
    metric = figPARAMS.filter.metric{n};
    range = figPARAMS.filter.stdDevRange{n};

    switch metric
        case 'Duration [cm]'
            [idx] = subFilter(range,STATS.cm.dur_s);
        case 'Duration [cm_max]'
            [idx] = subFilter(range,STATS.cm_max.dur_s);
        case 'AbsBandwidth [cm]'
            [idx] = subFilter(range,STATS.cm.abs_band_hz);
        case 'AbsBandwidth [cm_max]'
            [idx] = subFilter(range,STATS.cm_max.abs_band_hz);
        case 'NonBandwidth [cm]'
            [idx] = subFilter(range,STATS.cm.band_hz);
        case 'NonBandwidth [cm_max]'
            [idx] = subFilter(range,STATS.cm_max.band_hz);
        case 'Slope [cm]'
            [idx] = subFilter(range,STATS.cm.slope);
        case 'Slope [cm_max]'
            [idx] = subFilter(range,STATS.cm_max.slope);
        case 'Start Frequency [cm]'
            [idx] = subFilter(range,STATS.cm.start_hz);
        case 'Start Frequency [cm_max]'
            [idx] = subFilter(range,STATS.cm_max.start_hz);
        case 'End Frequency [cm]'
            [idx] = subFilter(range,STATS.cm.end_hz);
        case 'End Frequency [cm_max]'
            [idx] = subFilter(range,STATS.cm_max.end_hz);
        case 'Min Frequency [cm]'
            [idx] = subFilter(range,STATS.cm.min_hz);
        case 'Min Frequency [cm_max]'
            [idx] = subFilter(range,STATS.cm_max.min_hz);
        case 'Max Frequency [cm]'
            [idx] = subFilter(range,STATS.cm.max_hz);
        case 'Max Frequency [cm_max]'
            [idx] = subFilter(range,STATS.cm_max.max_hz);
        case 'Peak Frequency'
            [idx] = subFilter(range,STATS.cm.peak_freq);
        case 'P2P RL'
            [idx] = subFilter(range,STATS.cm.p2p_rl);
        case 'RMS RL [-3dB] [cm]'
            [idx] = subFilter(range,STATS.cm.rms_rl_m3);
        case 'RMS RL [-3dB] [cm_max]'
            [idx] = subFilter(range,STATS.cm_max.rms_rl_m3);
        case 'RMS RL [-10dB] [cm]'
            [idx] = subFilter(range,STATS.cm.rms_rl_m10);
        case 'RMS RL [-10dB] [cm_max]'
            [idx] = subFilter(range,STATS.cm_max.rms_rl_m10);
        case 'RMS RL [90%] [cm]'
            [idx] = subFilter(range,STATS.cm.rms_rl_e90);
        case 'RMS RL [90%] [cm_max]'
            [idx] = subFilter(range,STATS.cm_max.rms_rl_e90);
        case 'RMS RL [97%] [cm]'
            [idx] = subFilter(range,STATS.cm.rms_rl_e97);
        case 'RMS RL [97%] [cm_max]'
            [idx] = subFilter(range,STATS.cm_max.rms_rl_e97);
        case 'RMS RL [GPL] [cm]'
            [idx] = subFilter(range,STATS.cm.rms_rl_gplSE);
        case 'RMS RL [GPL] [cm_max]'
            [idx] = subFilter(range,STATS.cm_max.rms_rl_gplSE);
        case 'SEL [-3dB] [cm]'
            [idx] = subFilter(range,STATS.cm.SEL_m3);
        case 'SEL [-3dB] [cm_max]'
            [idx] = subFilter(range,STATS.cm_max.SEL_m3);
        case 'SEL [-10dB] [cm]'
            [idx] = subFilter(range,STATS.cm.SEL_m10);
        case 'SEL [-10dB] [cm_max]'
            [idx] = subFilter(range,STATS.cm_max.SEL_m10);
        case 'SEL [90%] [cm]'
            [idx] = subFilter(range,STATS.cm.SEL_e90);
        case 'SEL [90%] [cm_max]'
            [idx] = subFilter(range,STATS.cm_max.SEL_e90);
        case 'SEL [97%] [cm]'
            [idx] = subFilter(range,STATS.cm.SEL_e97);
        case 'SEL [97%] [cm_max]'
            [idx] = subFilter(range,STATS.cm_max.SEL_e97);
        case 'SEL [GPL] [cm]'
            [idx] = subFilter(range,STATS.cm.SEL_gplSE);
        case 'SEL [GPL] [cm_max]'
            [idx] = subFilter(range,STATS.cm_max.SEL_gplSE);
        otherwise
            disp('Invalid metric input, see readme for options')
            return
    end

    index_dump{n} = idx;

end

% Combine filtered calls
idx_final = index_dump{1};
for j = 1:length(index_dump)-1
    idx_final = intersect(idx_final,index_dump{j+1});
end

figCalls = plot_calls(idx_final);

end



