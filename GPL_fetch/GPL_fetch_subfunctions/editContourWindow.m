% editContourWindow.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script will edit the color values on the GPL contour window for
% cases where the normal imagesc route is not sufficient.
% Written: Ian 08/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Entering GPL Contour Settings')
subplot(3,1,3); hold on
disp('Switching contour to be based off of the colorbar')

hold off