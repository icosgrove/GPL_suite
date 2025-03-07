function breakoutClear(mode,state)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function clears out text indicators when the user breakout option
% 'b' is used.
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Breaking out of current special case mode.\n')
position1 = [0.3, -0.22, 0]; % pairing
position2 = [0.3, -0.3, 0]; % window

% Remove breakout state text
boxes1 = findobj(gca,'Type','text','Position',position1);
if ~isempty(boxes1) % clear previous pairing indicator
    if mode == 1
        delete(boxes1)
    end
end
boxes2 = findobj(gca,'Type','text','Position',position2);
if ~isempty(boxes2) % clear previous window indicator
    if mode == 2
        delete(boxes2)
    end
end

% Identify previous states
[strings] = IDstate(state);
specialModeIndicator(strings,1)