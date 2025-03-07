function [strings] = IDstate(state)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will control which text should be displayed when breaking
% out of one pairing/window state and reverting back to previous one.
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Identify previous pairing state
if ~isempty(state{1})
    switch state{1}
        case 'r'
            strings{1} = 'Special case pairing: Regular entry';
        case 'm'
            strings{1} = 'Special Case Mode: Multiple GPL Calls';
        case 'e'
            strings{1} = 'Special Case Mode: GPL Miss';
        case 'a'
            strings{1} = 'Special Case Mode: Adhoc';
        case 's'
            strings{1} = 'Special Case Mode: Detection Swap';
        case 'f'
            strings{1} = 'Special Case Mode: False Log Entry';
        case 'und'
            strings{1} = 'Special Case Mode: Undetermined';
        otherwise
            strings{1} = 'Special Case Mode: Idle';
    end
else
    strings{1} = 'Special Case Mode: Idle';
end

if ~isempty(state{2})
    switch state{2}
        case 'rp1'
            strings{1} = strcat(strings{1},' + Reprocess: Bad Box');
        case 'rp2'
            strings{1} = strcat(strings{1},' + Reprocess: Weak/Empty');
    end
end
