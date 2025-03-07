function [matrix] = GPL_full(label,k,GPL_struct)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPL_full will re-construct call contours that are stored in a 16-bit,
% compressed format from GPL_sparse. The values are re-scaled and placed
% back in the correct place to represent the contour of the call. 

% Written by Tyler Helble
% Tested, documented, and cleaned by Ian, but not fundamentally modified
% 04/01/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Recall contour data back into the workspace 
eval(strcat('values=double(GPL_struct(k).',label,'.values);')); % sp values of contour
eval(strcat('index=double(GPL_struct(k).',label,'.index);')); % contour indices
eval(strcat('peak=GPL_struct(k).',label,'.scale;')); % peak energy (for re-scaling)
eval(strcat('siz=GPL_struct(k).',label,'.size;')); % contour size 

% Recreate matrix of the contour
matrix = zeros(siz(1),siz(2)); % Size of the contour
matrix(index) = values*peak/(2^16-1); % Values re-scaled and placed into index location