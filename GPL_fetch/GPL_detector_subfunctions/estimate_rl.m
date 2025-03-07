function [spec_rl] = estimate_rl(cm,sp_orig)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will estimate the recieved level of a contour vs. the
% full window that the contour is found in.

% Written by Tyler Helble
% Documented, tested, and cleaned by Ian, but not fundamentally modified.
% 04/02/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Create a new slate containing the contour only, with values from original
% spectrogram, then convolve (smooth)
[contour_ind] = find(cm); 
new_slate = zeros(size(cm)); 
new_slate(contour_ind) = sp_orig(contour_ind); % Set contour vals using original spectrogram
pad = conv2(new_slate,ones(3,3)/9,'same'); % Convolution

% Create a noise slate of the region NOT part of the contour only.
noise_slate = sp_orig; 
noise_slate(pad ~= 0) = 0;

% Find indices of noise ONLY, vector of noise values, indices of contour
% ONLY
[x1,y1] = find(noise_slate);
z1 = noise_slate(noise_slate ~= 0);
[x2,y2] = find(cm);

% Use griddata to remove bins that fall around the boundary of the contour
z2 = griddata(x1,y1,z1,x2,y2,'nearest');

% Use that to reduce the values of the original spectrogram.
if ~isempty(z2) % Rare case z2 comes back empty
    new_slate(contour_ind) = new_slate(contour_ind) - z2;
end

% Remove a few bins <0 on boundary.
new_slate(new_slate < 0) = 0;

% Set RL as RMS of summed values (summed over time)
spec_rl = sqrt(mean(sum(new_slate.^2,1))); 
