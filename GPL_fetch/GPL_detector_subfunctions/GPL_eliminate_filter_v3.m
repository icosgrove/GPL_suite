%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPL_eliminate_filter_v3 will remove detections from an output .mat file
% created by GPL_v3 in process_HARP_v3.m. Detections are removed depending
% on if they pass or fail certain contour measurement criteria. 
%
% Potential measurements:
%   - Contour Duration
%   - Contour Frequency Bandwidth
%   - Contour Frequency Start/End
%   - Contour Frequency Minimum/Maximum
%   - Contour Slope
%
% These measurements are taken for the contour 'cm' which is all islands of
% the detection, 'cm_max' which is the single strongest energy island, and
% 'cm_max2' which is the second strongest energy island. Measurements are
% already marked as a 1/0 in the GPL process, filter parameters should be
% modified in GPL_parameter_input_v3.m. This script simply removes
% detections that fail whichever filter(s) the user specifies.
%
% Tyler Helble originally had a certain process for removing false
% detections. That process can be used here by setting filter_v2 = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GPL_v2 eliminate filter
% Switch on/off GPL_v2 filter elimination
filter_v2 = 1;

if filter_v2 == 1
    
    % Loop over all detections
    for n = 1:length(hyd.detection.calls)
        
        % Check if cm_max filter data exists
        if isfield(hyd.detection.calls(n).filter_v3.cm_max)
        if isfield(hyd.detection.calls(n).filter_v3.cm_max.duration_sec)
        if isfield(hyd.detection.calls(n).filter_v3.cm_max.freq_bandwidth_hz)
        if isfield(hyd.detection.calls(n).filter_v3.cm_max.slope)
            
            % Check if cm_max passes duration, bandwidth, and slope
            if (hyd.detection.calls(n).filter_v3.cm_max.duration_sec + ...
                    hyd.detection.calls(n).filter_v3.cm_max.freq_bandwidth_hz + ...
                    hyd.detection.calls(n).filter_v3.cm_max.slope) == 3
                continue % Keep detection if it passes
                
            else % If it doesn't pass cm_max, check is cm_max2 passes
                
                % Check if cm_max2 filter data exists
                if isfield(hyd.detection.calls(n).filter_v3.cm_max2)
                if isfield(hyd.detection.calls(n).filter_v3.cm_max2.duration_sec)
                if isfield(hyd.detection.calls(n).filter_v3.cm_nax2.freq_bandwidth_hz)
                if isfield(hyd.detection.calls(n).filter_v3.cm_max2.slope)
                    
                    % Check if cm_max2 passes duration, bandwidth, and slope
                    if (hyd.detection.calls(n).filter_v3.cm_max2.duration_sec + ...
                            hyd.detection.calls(n).filter_v3.cm_max2.freq_bandwidth_hz + ...
                            hyd.detection.calls(n).filter_v3.cm_max2.slope) == 3
                        continue % Keep detection if it passes
                    else
                        hyd.detection.calls(n) = []; % Remove detection if it fails cm_max and cm_max2
                    end
                    
                end % End Field checks cm_max2
                end % ""
                end % ""
                end % ""
                
                hyd.detection.calls(n) = []; % Remove detection if it fails cm_max and cm_max2 doesn't exist
                
            end % If: cm_max filter checks
            
        end % End field checks cm_max
        end % ""
        end % ""
        end % ""
        
    end % For: Loop over all calls
    
    % Save new output
    save('filtered_detections_GPL_v2_35_90_NUNAT_SB_01_210811_190000_df100.mat','hyd')
    
end % If: Filter_v2 is on






%% GPL_v3 Eliminate Filter
% Remove detections by the filtering criteria of your choice

if filter_v2 == 0
    
    
    
end






















            
            
    