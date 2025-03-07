function [dec] = xwavSelectWarning(xwav_start,manual)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function wil check if the next detection pick lies outside, or close
% to outside of the current range of time covered by the selected xwavs.
% Written: Ian 08/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%

dec = [];
proceed = false;

% Sort manual time
loc = find(sort([manual xwav_start]) == manual);
if loc == 1 % Manual time is before xwav range 
    dispErrorFlag(15)
    while ~proceed
        dec = input('Enter "c" to proceed, or "b" to break out: ','s');
        if strcmp(dec,'c') || strcmp(dec,'b')
            proceed = true;
        else
            dispErrorFlag(16)
        end
    end
    return
elseif loc == length(xwav_start) + 1 % Manual time is after xwav range
    dispErrorFlag(15)
    while ~proceed
        dec = input('Enter "c" to proceed, or "b" to break out: ','s');
        if strcmp(dec,'c') || strcmp(dec,'b')
            proceed = true;
        else
            dispErrorFlag(16)
        end
    end
    return
else
    return
end
    
