function [island_area,island_energy_sum,island_start] = island21(bas,chain)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Island21 will compute energy sums for each individual island. This is
% used to determine which island is the strongest and second strongest
% contour, which are the primary modes for filtration via contour size
% measurements later.

% Written by Tyler Helble
% Tested, documented, and cleaned by Ian, but not fundamentally modified
% 03/28/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Pre-processing 
island_area(1) = chain(1); % Length (# of indices) of first (or only) island.
loc = 1; 
chain_length = length(chain); % total number of indices for all islands
     

% Loop through the chain of islands (Format: [length(island #1), indices of the
% spectrogram containing that island length(island #2) indices ...]) And 
% save size of each island only to 'island_area'
i = 2;
while (loc + island_area(i-1) + 1) < (chain_length + 1) % Break when indexed beyond number of islands

    loc = loc + island_area(i-1) + 1; % Locate index of the next island size value
    island_area(i) = chain(loc); % Save size of the island

    i = i + 1;
end
    

% Find the location in 'chain' pertaining to the actual first index of an
% island (not the size of the island)
island_start(1) = 1; 
for i = 1:length(island_area) % Loop over all islands
    island_start(i+1) = island_start(i) + island_area(i) + 1; % Start point is previous start + island size
end

% Sum all the energy pertaining to each island from the whitened
% spectrogram normalized over mean noise level of that sp.
for i = 1:length(island_area) % loop over all islands
    island_energy_sum(i) = sum(bas(chain(island_start(i)+1:island_start(i+1)-1))); 
end
