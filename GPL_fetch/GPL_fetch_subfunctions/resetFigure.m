function fig = resetFigure()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function checks in the GPL fetch figure has been created or not and
% removes it and resets it if needed. 
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tag = 'tag';
fig_handle = findobj('Type','figure','Tag',tag);
if isempty(fig_handle) % No figure, make one
    fig = figure(7);
    set(fig,'Tag',tag)
else % Clear figure and make new one
    clf(fig_handle)
    fig = figure(7);
    set(fig,'Tag',tag)
end