function chgContour(fetched_call,parm,b,c,zero_idx,fig,time_range,freq_range)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function applies changes to the GPL contour window
% Written: Ian 08/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set brightness and contrast
parm.fp1.gpl_brightness = b;
parm.fp1.gpl_contrast = c;
fetched_call.gpl_contour_window = (c/100).*fetched_call.gpl_contour_window + b;
fetched_call.gpl_contour_window(zero_idx) = parm.fp1.noise_floor;

% Change figure
figure(fig); subplot(3,1,3); hold on
image(time_range,freq_range,fetched_call.gpl_contour_window)
axs = gca;
look4text = findall(axs,'Type','text','Units','normalized','Position',[1 -0.18 0]);
delete(look4text)
BC_text(parm,2);
hold off

end