function [cm_max] = reprocCmMax(cm_max)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will check if cm_max still has multiple islands and return
% only one if it is the case.
% Using code written by Tyler Helble.
% Written: Ian 08/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Locate indices, size
maxI = find(cm_max);
[sz1,sz2] = size(cm_max);

% Redo the chain of islands and find max island
[chain] = island1x(maxI,sz1);
[~,island_energy_sum,island_start] = island21(cm_max,chain);
msk = zeros(sz1,sz2);  
[~,ks0] = sort(island_energy_sum);

% Set max island as cm_max
msk(chain(island_start(ks0(end))+1:island_start(ks0(end)+1)-1)) = 1;
cm_max = msk.*cm_max;