function [v] = plotDSFB_cm_max(STATS,v,histPARAMS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the Duration, slope, and frequency bandwidth for cm_max contours
% Written: Ian 09/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Duration (s)
figure(v); 
sgtitle('Cm_max: Single Strongest Contour','Interpreter','none')
subplot(2,2,1)
% [uniq_vals,~,index] = unique(dur_s);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
histogram(STATS.cm_max.dur_s,'BinWidth',histPARAMS.cm_max.dur.binWidth)
title('Duration'); xlabel('Time [s]'); ylabel('Counts')
box on; grid on

% Slope (Hz/s)
subplot(2,2,2)
histogram(STATS.cm_max.slope,'BinWidth',histPARAMS.cm_max.slope.binWidth)
title('Slope'); xlabel('Slope [Hz/s]'); ylabel('Counts')
box on; grid on

% Nonzero Bandwidth (Hz)
subplot(2,2,3)
% [uniq_vals,~,index] = unique(band_hz);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
histogram(STATS.cm_max.band_hz,'BinWidth',histPARAMS.cm_max.nonzFB.binWidth)
title('Nonzero Frequency Bandwidth'); xlabel('Bandwidth [Hz]'); ylabel('Counts')
box on; grid on

% Full Bandwidth (Hz)
subplot(2,2,4)
% [uniq_vals,~,index] = unique(abs_band_hz);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
histogram(STATS.cm_max.abs_band_hz,'BinWidth',histPARAMS.cm_max.fullFB.binWidth)
title('Absolute Frequency Bandwidth'); xlabel('Bandwidth [Hz]'); ylabel('Counts')
box on; grid on; v=v+1;