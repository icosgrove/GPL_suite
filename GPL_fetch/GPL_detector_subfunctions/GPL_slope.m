function [slope_ci, slope, slope_y_int] = GPL_slope(contour, parm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will compute a weighted slope of the given contour using
% functions from the Curve Fitting Toolbox. The outputs include a
% confidence interval and a y-intercept of the slope to be used for
% plotting in a seperate function.

% Using code written by Tyler Helble
% Created, documented, and tested by Ian
% 04/02/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find number of time bins of the contour
nonz = find(sum(contour)); % nonzero columns 
num_Tbins = length(nonz); 

if num_Tbins > 2 % Length check

    % Find contour indices
    [contour_y,contour_x,contour_nonz] = find(contour);

    % Compute linear fit slope for the contour
    [fitresult, ~] = weighted_slope(contour_x, contour_y, contour_nonz.^2);

    % Compute confidence of the fit
    ci = confint(fitresult,0.95);
    confidence = abs(ci(1)-ci(2)); 

    % Save relevant otuputs
    slope_ci = confidence;
    slope = fitresult.p1; 
    slope_y_int = fitresult.p2 + parm.FreqBinLo - 1;
    
else
    
    slope_ci = nan;
    slope = nan;
    slope_y_int = nan;

end