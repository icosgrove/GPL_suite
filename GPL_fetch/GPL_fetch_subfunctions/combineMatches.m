function [combined] = combineMatches(multiple_calls,parm,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will combine together two GPL detection calls, used as a
% subfunction for when three or more calls must be combined. 
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sample start/end
combined.start_time = multiple_calls(1).start_time;
combined.end_time = multiple_calls(2).end_time;

% Remove empty contours and combine
combined.cm = removeEmptyContours(multiple_calls(1).cm,multiple_calls(2).cm,parm,n);
combined.cm_max = removeEmptyContours(multiple_calls(1).cm_max,multiple_calls(2).cm_max,parm,n);
combined.cm_max2 = removeEmptyContours(multiple_calls(1).cm_max2,multiple_calls(2).cm_max2,parm,n);  

% Waveform/spectral noise/recieved level
if isnan(multiple_calls(1).cm_max_waveform) 
    if isnan(multiple_calls(2).cm_max_waveform) 
        combined.cm_max_waveform = nan; % Both nan
    else
        combined.cm_max_waveform = multiple_calls(2).cm_max_waveform; % 2 is valid, 1 fails
    end
else % 1 passes
    if isnan(multiple_calls(2).cm_max_waveform)  % 1 is valid, 2 fails
        combined.cm_max_waveform = multiple_calls(1).cm_max_waveform;
    else
        combined.cm_max_waveform = [multiple_calls(1).cm_max_waveform; multiple_calls(2).cm_max_waveform]; % Both pass
    end
end          
if isnan(multiple_calls(1).spec_noise) && isnan(multiple_calls(2).spec_noise)
    combined.spec_noise = nan;
else
    combined.spec_noise = mean([multiple_calls(~isnan([multiple_calls(:).spec_noise])).spec_noise]);
end
if isnan(multiple_calls(1).spec_rl) && isnan(multiple_calls(2).spec_rl)
    combined.spec_rl = nan;
else
    combined.spec_rl = mean([multiple_calls(~isnan([multiple_calls(:).spec_rl])).spec_rl]);
end

% Julian start/end
combined.julian_start_time = multiple_calls(1).julian_start_time;
combined.julian_end_time = multiple_calls(2).julian_end_time;