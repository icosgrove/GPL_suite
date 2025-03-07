function [adhoc_append] = appendAdhocDetection(adhoc_append,adhoc_detection)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will append selected Adhox detections to the current
% running total. A function is used because multiple Adhoc need to be
% handled carefully
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Case: First set of adhocs
if isempty(adhoc_append)
    adhoc_append = adhoc_detection;
    return
end

% Append adhocs one by one
if ~isempty(adhoc_detection)
    for n = 1:length(adhoc_detection)
        adhoc_append = [adhoc_append; adhoc_detection(n)];
    end
end