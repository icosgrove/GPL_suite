% setColorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup the colorbar and y2 axis stuff
% Adapted from Triton code.
% Written: Ian 10/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if figPARAMS.colorbar.on == 1 % Colorbar
    axes1.FontSize = figPARAMS.y2FontSize;
    if figPARAMS.y2labelBold == 1
        axes1.FontWeight = 'bold';
    end

    figPARAMS.colorbar.cb = colorbar('vert');
    yl = get(figPARAMS.colorbar.cb,'YLabel');
    set(yl,'String',figPARAMS.colorbar.YLabelString);
    minc = figPARAMS.colorbar.sp_min;
    maxc = figPARAMS.colorbar.sp_max;
    set(figPARAMS.colorbar.cb,'YLim',[minc maxc])
    cbYtick = get(figPARAMS.colorbar.cb,'Ytick');  % current colorbar ticks
    newYtick = (cbYtick - figPARAMS.brightness) / (figPARAMS.contrast/100);  % reverse of bright/constrast
    set(figPARAMS.colorbar.cb,'YtickLabel',num2str(newYtick'));  % new labels
end