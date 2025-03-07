function [v] = plotSEL_cm_max(STATS,v,histPARAMS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the sound exposure level for cm_max contours
% Written: Ian 09/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SEL (cm_max)
figure(v); 
sgtitle('Sound Exposure Level: Cm_max','Interpreter','none')
subplot(2,3,1)
histogram([STATS.cm_max.SEL_m3],'BinWidth',histPARAMS.cm_max.sel3.binWidth)
title('-3dB Endpoints'); xlabel('RL [dB re 1 \muPa^2/Hz]'); ylabel('Counts')
box on; grid on
subplot(2,3,4)
histogram([STATS.cm_max.SEL_m10],'BinWidth',histPARAMS.cm_max.sel10.binWidth)
title('-10dB Endpoints'); xlabel('RL [dB re 1 \muPa^2/Hz]'); ylabel('Counts')
box on; grid on
subplot(2,3,2)
histogram([STATS.cm_max.SEL_e90],'BinWidth',histPARAMS.cm_max.sel90.binWidth)
title('90% Total Energy Endpoints'); xlabel('RL [dB re 1 \muPa^2/Hz]'); ylabel('Counts')
box on; grid on
subplot(2,3,5)
histogram([STATS.cm_max.SEL_e97],'BinWidth',histPARAMS.cm_max.sel97.binWidth)
title('97% Total Energy Endpoints'); xlabel('RL [dB re 1 \muPa^2/Hz]'); ylabel('Counts')
box on; grid on
subplot(2,3,3)
histogram([STATS.cm_max.SEL_gplSE],'BinWidth',histPARAMS.cm_max.selGPL.binWidth)
title('GPL-marked Endpoints'); xlabel('RL [dB re 1 \muPa^2/Hz]'); ylabel('Counts')
box on; grid on; v=v+1;