function [call_id,breakoutFlag] = inputCallID(string,GPL_struct,allowZero)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will requets user input for different parts of the fetch
% process and check if the input is an integer > (&sometimes=) 0. It will
% catch errors and force re-dos until the request is valid. 
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

valid = false; error_flag = 1;
while error_flag ~= 0 % Catch Bad inputs
    while ~valid % Catch input exceptions
        try
            call_id = input(string, 's'); % Pairing input
            if strcmp(call_id,'b') % Breakout of current special input mode
                call_id = [];
                breakoutFlag = 1;
                return
            else
                breakoutFlag = 0;
            end
            call_id = str2double(call_id);
            valid = true;
        catch
            dispErrorFlag(11);
            valid = false;
        end
    end
    
    % Error check for inputs
    [error_flag] = checkCallID(call_id,GPL_struct);
    if error_flag == 0 % Deal with 0
        if allowZero == 0
            if call_id == 0
                error_flag = 13;
            end
        end
    end
    dispErrorFlag(error_flag);
    
    % Re-try input if error check doesn't pass
    if error_flag ~= 0
        valid = false;
    end

end