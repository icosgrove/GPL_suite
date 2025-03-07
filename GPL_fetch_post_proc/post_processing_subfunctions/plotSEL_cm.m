function [v] = plotSEL_cm(STATS,v,histPARAMS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%_cm%%%%%%%%%%%%%%%
% Plot the sound exposure level for cm contours
% Written: Ian 09/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SEL (cm)
figure(v); 
sgtitle('Sound Exposure Level: Cm')
subplot(2,3,1)
histogram([STATS.cm.SEL_m3],'BinWidth',histPARAMS.cm.sel3.binWidth)
title('-3dB Endpoints'); xlabel('RL [dB re 1 \muPa^2/Hz]'); ylabel('Counts')
box on; grid on
subplot(2,3,4)
histogram([STATS.cm.SEL_m10],'BinWidth',histPARAMS.cm.sel10.binWidth)
title('-10dB Endpoints'); xlabel('RL [dB re 1 \muPa^2/Hz]'); ylabel('Counts')
box on; grid on
subplot(2,3,2)
histogram([STATS.cm.SEL_e90],'BinWidth',histPARAMS.cm.sel90.binWidth)
title('90% Total Energy Endpoints'); xlabel('RL [dB re 1 \muPa^2/Hz]'); ylabel('Counts')
box on; grid on
subplot(2,3,5)
histogram([STATS.cm.SEL_e97],'BinWidth',histPARAMS.cm.sel97.binWidth)
title('97% Total Energy Endpoints'); xlabel('RL [dB re 1 \muPa^2/Hz]'); ylabel('Counts')
box on; grid on
subplot(2,3,3)
histogram([STATS.cm.SEL_gplSE],'BinWidth',histPARAMS.cm.selGPL.binWidth)
title('GPL-marked Endpoints'); xlabel('RL [dB re 1 \muPa^2/Hz]'); ylabel('Counts')
box on; grid on; v=v+1;
