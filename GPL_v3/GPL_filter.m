function [GPL_struct] = GPL_filter(GPL_struct,parm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will run a binary filter on contour measurements produced
% in a prior function. If the contour does not pass the minimum/maximum
% requirements, a zero will be placed into the output structure. Detections
% given zeros can be removed with eliminate_filter.m
% Note that start/end frequency and minimum/maximum frequency are tailored
% in this filter for downsweeps only.

% Based on code written by Tyler Helble
% Expanded, modified, and documented by Ian
% 04/03/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Loop over all current detections 
for k = 1:length(GPL_struct) 

    % Filters for cm: all islands
    if isfield(GPL_struct,'cm')
        
        % Contour Duration:
        if GPL_struct(k).cm.duration_sec < parm.filter_parm.values.cm_min_duration_s
            GPL_struct(k).filter_v3.cm.duration_sec = 0;
        else
            GPL_struct(k).filter_v3.cm.duration_sec = 1;
        end
        % Frequency Bandwidth:
        if GPL_struct(k).cm.freq_bandwidth_hz < parm.filter_parm.values.cm_min_freq_bandwidth_hz
            GPL_struct(k).filter_v3.cm.freq_bandwidth_hz = 0;
        else
            GPL_struct(k).filter_v3.cm.freq_bandwidth_hz = 1;
        end
        % Start Frequency:
        if GPL_struct(k).cm.start_freq_hz < parm.filter_parm.values.cm_start_freq_hz
            GPL_struct(k).filter_v3.cm.start_freq_hz = 0;
        else
            GPL_struct(k).filter_v3.cm.start_freq_hz = 1;
        end
        % End Frequency
        if GPL_struct(k).cm.end_freq_hz > parm.filter_parm.values.cm_end_freq_hz
            GPL_struct(k).filter_v3.cm.end_freq_hz = 0;
        else
            GPL_struct(k).filter_v3.cm.end_freq_hz = 1;
        end
        % Minimum Frequency:
        if GPL_struct(k).cm.min_freq_hz > parm.filter_parm.values.cm_min_freq_hz
            GPL_struct(k).filter_v3.cm.min_freq_hz = 0;
        else
            GPL_struct(k).filter_v3.cm.min_freq_hz = 1;
        end
        % Maximum Frequency:
        if GPL_struct(k).cm.max_freq_hz < parm.filter_parm.values.cm_max_freq_hz
            GPL_struct(k).filter_v3.cm.max_freq_hz = 0;
        else
            GPL_struct(k).filter_v3.cm.max_freq_hz = 1;
        end
        % Slope:
        if (GPL_struct(k).cm.slope > parm.filter_parm.values.cm_slope_lower) && ...
                (GPL_struct(1).cm.slope < parm.filter_parm.values.cm_slope_upper)
            GPL_struct(k).filter_v3.cm.slope = 1;
        else
            GPL_struct(k).filter_v3.cm.slope = 0;
        end
        
    end % If: cm exists
    
    
    
    % Filters for cm_max: Single Strongest Contour
    if isfield(GPL_struct,'cm_max')
        
        % Contour Duration:
        if GPL_struct(k).cm_max.duration_sec < parm.filter_parm.values.cm_max_min_duration_s
            GPL_struct(k).filter_v3.cm_max.duration_sec = 0;
        else
            GPL_struct(k).filter_v3.cm_max.duration_sec = 1;
        end
        % Frequency Bandwidth:
        if GPL_struct(k).cm_max.freq_bandwidth_hz < parm.filter_parm.values.cm_max_min_freq_bandwidth_hz
            GPL_struct(k).filter_v3.cm_max.freq_bandwidth_hz = 0;
        else
            GPL_struct(k).filter_v3.cm_max.freq_bandwidth_hz = 1;
        end
        % Start Frequency:
        if GPL_struct(k).cm_max.start_freq_hz < parm.filter_parm.values.cm_max_start_freq_hz
            GPL_struct(k).filter_v3.cm_max.start_freq_hz = 0;
        else
            GPL_struct(k).filter_v3.cm_max.start_freq_hz = 1;
        end
        % End Frequency
        if GPL_struct(k).cm_max.end_freq_hz > parm.filter_parm.values.cm_max_end_freq_hz
            GPL_struct(k).filter_v3.cm_max.end_freq_hz = 0;
        else
            GPL_struct(k).filter_v3.cm_max.end_freq_hz = 1;
        end
        % Minimum Frequency:
        if GPL_struct(k).cm_max.min_freq_hz > parm.filter_parm.values.cm_max_min_freq_hz
            GPL_struct(k).filter_v3.cm_max.min_freq_hz = 0;
        else
            GPL_struct(k).filter_v3.cm_max.min_freq_hz = 1;
        end
        % Maximum Frequency:
        if GPL_struct(k).cm_max.max_freq_hz < parm.filter_parm.values.cm_max_max_freq_hz
            GPL_struct(k).filter_v3.cm_max.max_freq_hz = 0;
        else
            GPL_struct(k).filter_v3.cm_max.max_freq_hz = 1;
        end
        % Slope:
        if (GPL_struct(k).cm_max.slope > parm.filter_parm.values.cm_max_slope_lower) && ...
                (GPL_struct(1).cm_max.slope < parm.filter_parm.values.cm_max_slope_upper)
            GPL_struct(k).filter_v3.cm_max.slope = 1;
        else
            GPL_struct(k).filter_v3.cm_max.slope = 0;
        end
        
    end % If: cm_max exists
        
    
    
    % Filters for cm_max2: Single Second Strongest Contour
    if isfield(GPL_struct,'cm_max2')
        
        % Contour Duration:
        if GPL_struct(k).cm_max2.duration_sec < parm.filter_parm.values.cm_max2_min_duration_s
            GPL_struct(k).filter_v3.cm_max2.duration_sec = 0;
        else
            GPL_struct(k).filter_v3.cm_max2.duration_sec = 1;
        end
        % Frequency Bandwidth:
        if GPL_struct(k).cm_max2.freq_bandwidth_hz < parm.filter_parm.values.cm_max2_min_freq_bandwidth_hz
            GPL_struct(k).filter_v3.cm_max2.freq_bandwidth_hz = 0;
        else
            GPL_struct(k).filter_v3.cm_max2.freq_bandwidth_hz = 1;
        end
        % Start Frequency:
        if GPL_struct(k).cm_max2.start_freq_hz < parm.filter_parm.values.cm_max2_start_freq_hz
            GPL_struct(k).filter_v3.cm_max2.start_freq_hz = 0;
        else
            GPL_struct(k).filter_v3.cm_max2.start_freq_hz = 1;
        end
        % End Frequency
        if GPL_struct(k).cm_max2.end_freq_hz > parm.filter_parm.values.cm_max2_end_freq_hz
            GPL_struct(k).filter_v3.cm_max2.end_freq_hz = 0;
        else
            GPL_struct(k).filter_v3.cm_max2.end_freq_hz = 1;
        end
        % Minimum Frequency:
        if GPL_struct(k).cm_max2.min_freq_hz > parm.filter_parm.values.cm_max2_min_freq_hz
            GPL_struct(k).filter_v3.cm_max2.min_freq_hz = 0;
        else
            GPL_struct(k).filter_v3.cm_max2.min_freq_hz = 1;
        end
        % Maximum Frequency:
        if GPL_struct(k).cm_max2.max_freq_hz < parm.filter_parm.values.cm_max2_max_freq_hz
            GPL_struct(k).filter_v3.cm_max2.max_freq_hz = 0;
        else
            GPL_struct(k).filter_v3.cm_max2.max_freq_hz = 1;
        end
        % Slope:
        if (GPL_struct(k).cm_max2.slope > parm.filter_parm.values.cm_max2_slope_lower) && ...
                (GPL_struct(1).cm_max2.slope < parm.filter_parm.values.cm_max2_slope_upper)
            GPL_struct(k).filter_v3.cm_max2.slope = 1;
        else
            GPL_struct(k).filter_v3.cm_max2.slope = 0;
        end
        
    end % If: cm_max2 exists
        

    
end % For: loop over all detections from current window.