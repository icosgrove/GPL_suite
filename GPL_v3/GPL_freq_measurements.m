function [start_freq, end_freq, min_freq, max_freq, abs_bandwidth] = GPL_freq_measurements(contour,parm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will measure frequency characteristics for a given contour.
% The start and end frequency is the mean frequency value corresponding the
% leftmost and rightmost side of the contour (beginning and end in time).
% The maximum and minimum frequency correspond to the absolute max and min
% frequency values of the contour, regardless of where they occur in time.

% Written and documented by Ian Cosgrove
% Based on call contours created using code written by Tyler Helble
% 04/03/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Frequency Start/End
% Sum energy over time and locate nonzero bins
nonz_time = find(sum(contour));

% Locate start and end time bins of the contour
contour_start_bin = min(nonz_time);
contour_end_bin = max(nonz_time);

% Calculate mean frequency of these start/end bins
start_freq = mean(find(contour(:,contour_start_bin))) + parm.FreqBinLo;
start_freq = start_freq * (parm.SampleFreq/parm.fftl) - 2; % Convert to Hz

end_freq = mean(find(contour(:,contour_end_bin))) + parm.FreqBinLo;
end_freq = end_freq * (parm.SampleFreq/parm.fftl) - 2; % Convert to Hz


%% Frequency Min/Max
% Sum energy over frequecy and find nonzero bins
nonz_freq = find(sum(contour,2));

% Find largest/smallest frequency value
max_freq = max(nonz_freq) + parm.FreqBinLo - 2; % Account for bin 0 and indexing starting at 1
max_freq = max_freq * (parm.SampleFreq/parm.fftl); % Convert to Hz

min_freq = min(nonz_freq) + parm.FreqBinLo - 2;
min_freq = min_freq * (parm.SampleFreq/parm.fftl); % Convert to Hz

abs_bandwidth = max_freq - min_freq + 1;
