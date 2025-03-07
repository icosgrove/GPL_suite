function [fetched_call] = clearSavedWindows(fetched_call,parm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script erases saved windows from output of GPL Fetch
% Written: Ian 08/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear out saved windows in output if necessary
if parm.fp1.spectrogram == 0
    fetched_call.spectrogram = [];
end
if parm.fp1.gpl_contour_window == 0
    fetched_call.gpl_contour_window = [];
end
if parm.fp1.cm_maxsubplot == 1
    if parm.fp1.savecm_max == 0
        fetched_call.cm_max_window = [];
    end
end