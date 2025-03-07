function [single_save] = removeEmptySingleContour(single_calls)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will go through manual calls matched to a single GPL
% detection and remove those that have non-empty contours.
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through and save detections that have nonempty indices
k = 1;
for n = 1:length(single_calls)  
    if ~isempty(single_calls(n).gpl_match.cm.index)
        single_save.calls(k) = single_calls(n);
        k=k+1;
    end
end