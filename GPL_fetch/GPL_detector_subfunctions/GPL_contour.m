function [cm,cm_max,cm_max2] = GPL_contour(bas0,parm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPL_contour will take each detection individually, perform manipulations
% and identify where in time and frequency its energy spikes occurred. The
% contour 'cm' is all of these high-energy pixels. Pixels are separated from
% each other by regions of low energy, creating islands. 'cm_max' is the
% island with the most energy, 'cm_max2' is the island with the second most
% amount of energy. Measurements are later taken of these contours to
% determine if the detection is a signal of interest or not.

% Written by Tyler Helble
% Tested, documented, and cleaned by Ian but not fundamentally modified.
% 02/10/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Allocate 
% Add a one bin buffer around the call
[sz1,sz2] = size(bas0); 
bas = zeros(sz1+2,sz2+2); 
bas(2:end - 1,2:end - 1) = bas0; 

% Allocate Contour matrices 
cm = 0*bas;
cm_max = cm;
cm_max2 = cm;
ks0 = 0;


%% Apply whitening and normalization to the call sp data
% Reshape the call spectrogram into a single vector 
[sz1,sz2] = size(bas); 
basgot = reshape(bas,sz1*sz2,1); 

% Whiten the vector and normalize by it's mean noise level
[~,mu] = whiten_vec(basgot); 
qre = bas/mu - 1; 

% If not provided in input args, define the pizel energy cutoff
if nargin == 1 
    cutoff = 10; 
end


%% Identify Islands

% Define minimum number of pixels for a detection to be considered a contour.
a_min = parm.MinIslandSize;

% Locate indices of whitened/normalized pixel levels above the energy
% cutoff. These are pixels considered as part of the contour
k = find(qre > parm.ContourCutoff); 


% Identify islands of the contour if it passes minimum length filter
if length(k) > a_min 
    
    % Set passed pixel indices to '1' in a window the size of the call
    msk0 = zeros(sz1,sz2); 
    msk0(k) = 1; 
    
    % Perform a 2-D convolution using the central part on the passed pixels
    dm = msk0; 
    dm = conv2(dm,ones(3,3)/9,'same');
    
    % get rid of low-lying pixels after convolution
    dm(dm > 2/9) = 1; 
    dm(dm < 1) = 0; 
    
    % Apply the convolution (with some pixels on the edges removed in
    % previous step) to the original pixels. This has the effect of adding
    % some extra pixels above cutoff around the call.
    dm = dm + msk0; 

    % Locate new set of pixel indices
    kk = find(dm); 

    
    % Locate islands in the call. Choice of expanded radius available. An
    % island is a region of the contour with enough energy to be considered
    % relevant to the overall energy sum. Each island is separated from
    % others on all sides by regions of lower energy (presumed not to
    % contribute to the contour)
    if parm.island_radius == 1
        chain = island1x(kk,sz1); % Each pixel's direct neighbors are checked for high energy
    else
        chain = radial_island(kk,sz1,sz2,parm); % A customizable-sized region is checked around each pixel
    end

    
    % Process each island: Extract it's size, where it begins in
    % 'chain', and the each individual energy sum. (Energy sum is taken
    % from the whitened spectrogram normalized over the sp's mean noise)
    [island_area,island_energy_sum,island_start] = island21(bas,chain);

    
    %% Identify Contours: cm, cm_max, cm_max2
    
    % Proceed with contour creation for the strongest islands
    if max(island_area) > a_min % Check for at least one viable island
        
        % Consider all islands above a small size cutoff percentage of the
        % largest island.
        ks = find(island_area > parm.smaller_island_cutoff*max(island_area)); 

        % Sort the islands by energy sum value
        msk = zeros(sz1,sz2);  
        [~,ks0] = sort(island_energy_sum); % Sorted weakest-strongest by island number
        
        % Locate the indices of the original detection for strongest island
        msk(chain(island_start(ks0(end))+1:island_start(ks0(end)+1)-1)) = 1;
        % Single contour for island with most energy
        cm_max = msk.*bas;

        % Reprocess to fix cases where cm_max returns two islands
        [cm_max] = reprocCmMax(cm_max);
     
        % Process a second strongest contour if >1 islands exist
        if length(ks0) > 1 
            
            msk = zeros(sz1,sz2);
            
            % Indices of second strongest island only
            msk(chain(island_start(ks0(end-1))+1:island_start(ks0(end-1)+1)-1)) = 1;
            % Single contour for island with second most amount of energy
            cm_max2 = msk.*bas; 
            
        end
     
        % Process ALL islands found in the detection in a single contour (size cutoff still applies) 
        msk = zeros(sz1,sz2); 
        for cm_loop = 1:length(ks) % Loop over all islands
            msk(chain(island_start(ks(cm_loop))+1:island_start(ks(cm_loop)+1)-1)) = 1; % combine the contours together 
        end
        % Single contour with all islands in it.
        cm = msk.*bas; 

        
    end % End if: at least one contour possible

end % End if: detection is large enough to be considered for contours

% Remove 1 bin padding around contour spectrograms. (1 bin radius only)
if parm.island_radius == 1
    cm = cm(2:end-1,2:end-1); 
    cm_max = cm_max(2:end-1,2:end-1); 
    if length(ks0) > 1
        cm_max2 = cm_max2(2:end-1,2:end-1); 
    end
end
