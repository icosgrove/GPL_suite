function areCompat = structcompatable(s1, s2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will compare the length of fields of two structures and
% return logical true if they may be combined (have the same set of fields)
% or logical false if they do not. 
% Written by Ian Cosgrove 06/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = fieldnames(s1);
f2 = fieldnames(s2);

if length(f1) == length(f2) 
    areCompat = all(strcmp(sort(f1),sort(f2)));
else
    areCompat = 0;
end