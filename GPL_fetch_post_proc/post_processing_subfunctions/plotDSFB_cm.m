function [v] = plotDSFB_cm(STATS,v,histPARAMS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the Duration, slope, and frequency bandwidth for cm contours
% Written: Ian 09/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Duration (s)
figure(v); 
sgtitle('Cm: Full Contour')
subplot(2,2,1)
% [uniq_vals,~,index] = unique(dur_s);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
histogram(STATS.cm.dur_s,'BinWidth',histPARAMS.cm.dur.binWidth)
title('Duration'); xlabel('Time [s]'); ylabel('Counts')
box on; grid on

% Slope (Hz/s)
subplot(2,2,2)
histogram(STATS.cm.slope,'BinWidth',histPARAMS.cm.slope.binWidth)
title('Slope'); xlabel('Slope [Hz/s]'); ylabel('Counts')
box on; grid on

% Nonzero Bandwidth (Hz)
subplot(2,2,3)
% [uniq_vals,~,index] = unique(band_hz);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
histogram(STATS.cm.band_hz,'BinWidth',histPARAMS.cm.nonzFB.binWidth)
title('Nonzero Frequency Bandwidth'); xlabel('Bandwidth [Hz]'); ylabel('Counts')
box on; grid on

% Full Bandwidth (Hz)
subplot(2,2,4)
% [uniq_vals,~,index] = unique(abs_band_hz);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
histogram(STATS.cm.abs_band_hz,'BinWidth',histPARAMS.cm.fullFB.binWidth)
title('Absolute Frequency Bandwidth'); xlabel('Bandwidth [Hz]'); ylabel('Counts')
box on; grid on; v=v+1;