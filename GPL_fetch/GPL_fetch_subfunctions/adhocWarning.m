function adhocWarning(fig,t_handle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will display a warning and visually indicate to the user
% that they shouldn't do adhocs yet because more calls are coming in the
% same window.
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Adhoc Status: WAIT, more detections coming in this window.\n')
figure(fig)
for n = 1:2
    set(t_handle,'Visible','off')
    pause(0.3)
    set(t_handle,'Visible','on')
    pause(0.3)
end