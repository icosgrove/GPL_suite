function [v] = plotSEMM_cm_max(STATS,v,histPARAMS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the start/end/min/mas frequency for cm_max contours
% Written: Ian 09/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start Frequency (Hz)
figure(v); 
sgtitle('Cm_max: Single Strongest Contour','Interpreter','none')
subplot(2,2,1) 
histogram(STATS.cm_max.start_hz,'BinWidth',histPARAMS.cm_max.startf.binWidth)
title('Start Frequency'); xlabel('Frequency [Hz]'); ylabel('Counts')
box on; grid on; 

% End Frequency (Hz)
subplot(2,2,2) 
histogram(STATS.cm_max.end_hz,'BinWidth',histPARAMS.cm_max.endf.binWidth)
title('End Frequency'); xlabel('Frequency [Hz]'); ylabel('Counts')
box on; grid on; 

% Min. Frequency (Hz)
subplot(2,2,3)
% [uniq_vals,~,index] = unique(min_hz);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
histogram(STATS.cm_max.min_hz,'BinWidth',histPARAMS.cm_max.minf.binWidth)
title('Minimum Frequency'); xlabel('Frequency [Hz]'); ylabel('Counts')
box on; grid on; 

% Max. Frequency (Hz)
subplot(2,2,4)
% [uniq_vals,~,index] = unique(max_hz);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
histogram(STATS.cm_max.max_hz,'BinWidth',histPARAMS.cm_max.maxf.binWidth)
title('Maximum Frequency'); xlabel('Frequency [Hz]'); ylabel('Counts')
box on; grid on; v=v+1;