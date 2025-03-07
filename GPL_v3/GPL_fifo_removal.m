function [sp] = GPL_fifo_removal(sp_raw, parm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will remove 50 Hz harmonics from the spectrogram. Those
% farmonics are typically unwanted noise produced by onboard electronics on
% a HARP. The function will locate the number of integer multiples of 50Hz
% occur in the requested frequency range, remove them, and relace them with
% the average of the bins directly above and/or below the harmonic. 

% Written, documented, and tested by Ian Cosgrove 
% 04/08/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract structure elements
fund_fifo = parm.fund_fifo_freq;
freq_lo = parm.freq_lo;
freq_hi = parm.freq_hi;

% Allocate output
sp = sp_raw;

%% Identify the if/what FIFO frequencies lie in the frequency range

% Identify maximum realized frequency from FFT
max_freq = ceil(parm.SampleFreq / 2);

% Terminate if fundamental FIFO is out of possible frequency range
if fund_fifo > max_freq
    return
end

% Identify the number of FIFO harmonics in the frequenyt range
num_fifo_harmonics = floor(max_freq / fund_fifo);

% Find which harmonics fall in frequency range
fifo_range = fund_fifo:fund_fifo:num_fifo_harmonics*fund_fifo; % Vector of possible harmonics
fifo_range = sort([freq_lo freq_hi fifo_range]); % sorted freq limits into possibilities
fifo_lo = find(fifo_range == freq_lo); % locate low freq position
fifo_hi = find(fifo_range == freq_hi); % locate high freq position

%% Process different cases - Single bin removal

% Case: No FIFO frequencies lie in the requested frequency range
if fifo_hi - fifo_lo == 1
    return
end

% Case: Min frequency is a FIFO band
if numel(fifo_lo) > 1
    if parm.fifo3bin == 0 % 1 bin removal
        avg = (sp_raw(fifo_range(fifo_lo(1))+1,:) + sp_raw(fifo_range(fifo_lo(1))-1,:)) ./ 2;
        sp(fifo_range(fifo_lo(1)),:) = avg;
    else % 3 bin removal
       avg = (sp_raw(fifo_range(fifo_lo(1))+2,:) + sp_raw(fifo_range(fifo_lo(1))-2,:)) ./ 2;
       avg_rep = [avg; avg; avg]; 
       sp(fifo_range(fifo_lo(1))-1:fifo_range(fifo_lo(1))+1,:) = avg_rep; 
    end
end

% Case: Max frequency is FIFO band, and possibly ceiling
if numel(fifo_hi) > 1
    if parm.fifo3bin == 0 % 1 bin removal
        if fifo_range(fifo_hi(1)) == max_freq % Max freq. limit is a FIFO band and the resolution ceiling
            sp(max_freq,:) = sp_raw(max_freq - 1,:);
        else % Max freq. limit is just a FIFO band
            avg = (sp_raw(fifo_range(fifo_hi(1))+1,:) + sp_raw(fifo_range(fifo_hi(1))-1,:)) ./ 2;
            sp(fifo_range(fifo_hi(1)),:) = avg;
        end
    else % 3 bin removal
        if fifo_range(fifo_hi(1)) == max_freq % Max freq. limit is a FIFO band and the resolution ceiling
            rep = [sp_raw(max_freq - 2,:); sp_raw(max_freq - 2,:)];
            sp(max_freq-1:max_freq,:) = rep;
        else % Max freq. limit is just a FIFO band
            avg = (sp_raw(fifo_range(fifo_hi(1))+2,:) + sp_raw(fifo_range(fifo_hi(1))-2,:)) ./ 2;
            avg_rep = [avg; avg; avg];
            sp(fifo_range(fifo_hi(1))-1:fifo_range(fifo_hi(1))+1,:) = avg_rep;
        end
    end
end

% Case: FIFO frequencies lie strictly within the requested limits
leftover = fifo_hi(1) - fifo_lo(end) - 1;
if leftover > 0
    if parm.fifo3bin == 0 % 1 bin removal
        for n = 1:leftover % loop over FIFO bands within limits
            replace_freq = fifo_range(fifo_lo(end) + n);
            avg = (sp_raw(replace_freq+1,:) + sp_raw(replace_freq-1,:)) ./2; % Average from above and below
            sp(replace_freq,:) = avg; % Replace
        end
    else % 3 bin removal
        for n = 1:leftover % loop over FIFO bands within limits
            replace_freq = fifo_range(fifo_lo(end) + n);
            avg = (sp_raw(replace_freq+2,:) + sp_raw(replace_freq-2,:)) ./2; % Average from above and below
            avg_rep = [avg; avg; avg];
            sp(replace_freq-1:replace_freq+1,:) = avg_rep; % Replace
        end
    end
end

end % Function 
