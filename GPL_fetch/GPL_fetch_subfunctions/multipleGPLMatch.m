function [fetched_call,breakoutFlag] = multipleGPLMatch(GPL_struct,fetched_call)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to mathc multiple GPL detections to one manual
% detection.
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[numdet,breakoutFlag] = inputCallID('Enter number of GPL detections that span the single call: ',GPL_struct,0);
if breakoutFlag == 1 % Terminate if breakout requested
    return
end
if isempty(numdet)
    dispErrorFlag(6)
    numdet = input('Enter number of GPL detections that span the single call: ');
end

% Loop over number of matches and record each GPL detection. 
for n1 = 1:numdet
    [dup_id,breakoutFlag] = inputCallID('Enter GPL Detection #s: ',GPL_struct,0);
    if breakoutFlag == 1 % Breakout option
        fetched_call.gpl_match = [];
        fetched_call.note = cell(1,1);
        return
    end
    if isempty(dup_id)
        dispErrorFlag(6)
        [dup_id,breakoutFlag] = inputCallID('Enter GPL Detection #s: ',GPL_struct,0);
        if breakoutFlag == 1 % Breakout option
            fetched_call.gpl_match = [];
            fetched_call.note = cell(1,1);
            return
        end
    end
    structInd(n1) = dup_id;
    if dup_id == 0 % GPL Miss for re-do case
        fetched_call.gpl_match = nan;
        fetched_call.note = {'GPL Miss'};
    end
    if dup_id > length(GPL_struct) % Over indexed struct
        error_flag = 1;
        break
    end
    if isempty(fetched_call.gpl_match)
        fetched_call.gpl_match = GPL_struct(dup_id);
    else
        fetched_call.gpl_match(n1) = GPL_struct(dup_id);
    end
    if numdet ~= 1
        fetched_call.note = {'Multiple GPL'};
    else
        fetched_call.note = {''};
    end
end
