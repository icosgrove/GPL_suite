function confirmBox(mstart_sec, mend_sec, parm, fetched_call, cross_flag, confirmFlag, fig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add yellow boxes to both subplots indicating chosen detection.
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(fig);

% Identify the number of subplots
figgcf = gcf;
axH = findall(figgcf,'Type','axes'); % Axis handles
subplotH = axH(arrayfun(@(h) isempty(get(h, 'Tag')), axH)); % Locate handles that may only be plots
nSubplots = length(subplotH);

if nSubplots == 3
    subplot(3,1,3); hold on
    [rectHandle] = det_boxes(nan, nan, mstart_sec, mend_sec, parm, fetched_call, cross_flag, confirmFlag, fig);
    hold off
elseif nSubplots == 2
    subplot(3,1,2); hold on
    [rectHandle] = det_boxes(nan, nan, mstart_sec, mend_sec, parm, fetched_call, cross_flag, confirmFlag, fig);
    hold off
end

% subplot(2,1,1); hold on
% det_boxes(nan, nan, mstart_sec, mend_sec, parm, fetched_call, cross_flag, confirmFlag)
% hold off

pause(parm.fp1.box_wt) % Pause for user to see box