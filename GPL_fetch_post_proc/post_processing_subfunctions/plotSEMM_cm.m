function [v] = plotSEMM_cm(STATS,v,histPARAMS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the start/end/min/mas frequency for cm contours
% Written: Ian 09/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start Frequency (Hz)
figure(v); 
sgtitle('Cm: Full Contour')
subplot(2,2,1) 
histogram(STATS.cm.start_hz,'BinWidth',histPARAMS.cm.startf.binWidth)
title('Start Frequency'); xlabel('Frequency [Hz]'); ylabel('Counts')
box on; grid on; 

% End Frequency (Hz)
subplot(2,2,2) 
histogram(STATS.cm.end_hz,'BinWidth',histPARAMS.cm.endf.binWidth)
title('End Frequency'); xlabel('Frequency [Hz]'); ylabel('Counts')
box on; grid on; 

% Min. Frequency (Hz)
subplot(2,2,3)
% [uniq_vals,~,index] = unique(min_hz);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
histogram(STATS.cm.min_hz,'BinWidth',histPARAMS.cm.minf.binWidth)
title('Minimum Frequency'); xlabel('Frequency [Hz]'); ylabel('Counts')
box on; grid on; 

% Max. Frequency (Hz)
subplot(2,2,4)
% [uniq_vals,~,index] = unique(max_hz);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
histogram(STATS.cm.max_hz,'BinWidth',histPARAMS.cm.maxf.binWidth)
title('Maximum Frequency'); xlabel('Frequency [Hz]'); ylabel('Counts')
box on; grid on; v=v+1;