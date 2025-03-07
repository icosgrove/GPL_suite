function [idx] = subFilter(range,vec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove calls based on metric
% Written: Ian 09/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vec_copy = vec(~isnan(vec)); % remove nans
SD = std(vec_copy);
if range(1) == 0
    min = min(vec_copy)-1;
else
    min = mean(vec_copy) + range(1)*SD;
end
if range(2) == inf
    max = max(vec_copy)+1;
else
    max = mean(vec_copy) + range(2)*SD;
end
if range(1) == 0
    if range(2) == inf
        idx = find(vec); % why would anyone do this
    else % < x
        idx = find(vec < max);
    end
else
    if range(2) == inf % > x
        idx = find(vec > min);
    else % within range
        idx = intersect(find(vec>min),find(vec<max));
    end
end