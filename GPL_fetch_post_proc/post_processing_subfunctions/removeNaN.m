function [STATS] = removeNaN(STATS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove any NaN entries from each measurement stat
% Written: Ian 09/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cm
fieldnames_cm = fields(STATS.cm);
for n = 1:length(fieldnames_cm)
    vec = STATS.cm.(fieldnames_cm{n});
    STATS.cm.(fieldnames_cm{n}) = vec(~isnan(vec));
end

% Cm_max
fieldnames_cm_max = fields(STATS.cm_max);
for n = 1:length(fieldnames_cm_max)
    vec = STATS.cm_max.(fieldnames_cm_max{n});
    STATS.cm_max.(fieldnames_cm_max{n}) = vec(~isnan(vec));
end