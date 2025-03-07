function [peak_freq] = computePeakFreq(data,parm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function uses the methods from Triton v1.95.20230315 to create a
% frequency spectrum of a subset of data. All credit goes to the creators
% of Triton, as this is an adaptation of their code.
% Written: Ian 09/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove DC offset
noDC_data = detrend(data,'constant');

% Use Welch's method to find PSD, find max in BP range
if length(data) >= parm.fftl
    win = hanning(parm.fftl);
    fftl = parm.fftl;
else
    win = hanning(length(data));
    fftl = length(data);
end
nOL = parm.fftl - round(parm.fftOverlap);
if nOL >= fftl
    nOL = fftl-1;
end
[Pxx,~] = pwelch(noDC_data,win,nOL,fftl,parm.SampleFreq);
Pxx = 10*log10(Pxx); % dBs
Pxx = Pxx(parm.bp_min_f:parm.bp_max_f); % within BP corners only
peak_freq = find(Pxx==max(Pxx)) + parm.bp_min_f - 1;