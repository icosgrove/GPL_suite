function [r,c] = getSPIndex(figPARAMS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create row, column index for each subplot
% Written: Ian 10/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:figPARAMS.subsPerFig
    r(k) = ceil(k/figPARAMS.subplotSize(2));
    c(k) = k-(figPARAMS.subplotSize(2)*(r(k)-1));
end