function [error_flag] = checkCallID(call_id,GPL_struct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checks if what the user entered is actually an integer, and allows them
% to proceed if so.
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(call_id)
    error_flag = 6;
    return
end

% Check if conversion was successful and the number is an integer
if isnumeric(call_id) && isreal(call_id) && isfinite(call_id) && (floor(call_id) == call_id) 
    if call_id >= 0
        if call_id > length(GPL_struct)
            error_flag = 12;
        else
            error_flag = 0;
        end
    else
        error_flag = 7;
    end
else
    error_flag = 7;
end
    
