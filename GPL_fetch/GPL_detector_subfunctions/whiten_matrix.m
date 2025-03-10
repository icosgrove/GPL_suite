function [sp_whiten,qs,mu] = whiten_matrix(sp,fac)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function whiten_matrix will produce a vector containing the mean
% noise level for each frequency bin, the original window spectrogram with
% the mean noise subtracted from each value, and the signal-absent window
% with mean noise subttratced. 

% The function calls a  subfunction 'base3x' which extracts indices of
% signal absense from the window. Those indices are used to find the mean
% noise level for each frequency bin.

% Written by Tyler Helble
% Tested, documented, and cleaned by Ian, but not fundamentally modified.
% 02/02/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% If the enhanced mean noise level (fac) is not defined, then set it to 1
if nargin==1 
    fac=1;
end

[sz1,sz2] = size(sp); 
% Allocate array with same # frequency bins and half the # of time bins
qs = zeros(sz1,ceil(sz2/2)+2) + nan; 



%% Extract signal-absent time bins using 'base3x' for each frequency bin

% Loop through each frequency bin
for i = 1:sz1 
    
    % Apply base3x to current frequency bin over ALL time bins
    [ks] = base3x(sp(i,:));
    
    % Use 'middle energy' indices from base3x to create 'qs'
    % Qs will be a spectrogram with signal-present time bins removed.
    qs(i,1:length(ks)) = sp(i,ks); 

end

% Remove NaN from 'qs'
sqs = sum(qs);
qs = qs(:,isfinite(sqs)); % eliminate NaN columns




%% Find mean noise level for each frequency bin, mean-reduced spectrogram, and mean-reduced signal-absent space.

[~,sz4] = size(qs); % number of time bins after base3x removed high/low energy

% Find the mean of each leftover frequency bins after base3x 
mu = mean(qs,2);

% The frequency means are expanded to cover all time bins, and then
% subtracted from the original spectrogram. Each time bin has the mean
% energy over it's frequency bin subtracted from it.
sp_whiten = sp - fac*(mu*ones(1,sz2));

% Frequency means are subtratced from the columns that they were found in,
% which leaves a very flat image of signal-absent space
qs = qs - fac*mu*ones(1,sz4);
