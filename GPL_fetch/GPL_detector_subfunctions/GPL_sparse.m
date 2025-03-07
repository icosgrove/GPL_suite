function [GPL_struct] = GPL_sparse(matrix,label,k,GPL_struct)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPL_sparse will take the contours 'cm', 'cm_max', 'cm_max2', convert the
% spectrogram values of the contours into 16-bit integers, indentifies the
% size of the spectrogram, its maximum energy value (for re-scaling back to
% original values in contour recreation, done in GPL_full), and saves them
% to GPL_struct. After this step the detection timestamps (in samples) and
% the corresponding contours are completed. 

% Written by Tyler Helble
% Tested, documented, and cleaned by Ian, but not fundamentally modified
% 04/01/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Locate contour indices/values
index=find(matrix); % all nonzero indices of the matrix, or the indices of the contour
values=matrix(index); % the spectrogram values of the contour

% Normalize values by the maximum energy value, scale to be 16-bit integers
peak = max(values); % Max Energy
ivalues = uint16(round(values/peak*(2^16-1)));

% Size of the contour (Frequency by time)
siz = size(matrix);

% Convert indices to form uint16
index = uint16(index);

% Load everything just created into GPL_struct
eval(strcat('GPL_struct(k).',label,'.values=ivalues;')); 
eval(strcat('GPL_struct(k).',label,'.index=index;'));
eval(strcat('GPL_struct(k).',label,'.scale=peak;'));
eval(strcat('GPL_struct(k).',label,'.size=siz;'));