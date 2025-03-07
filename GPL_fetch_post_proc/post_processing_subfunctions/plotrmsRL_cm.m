function [v] = plotrmsRL_cm(STATS,v,histPARAMS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the rms received level for the cm contours
% Written: Ian 09/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RMS RL (cm)
figure(v); 
sgtitle('RMS Received Level: Cm')
subplot(2,3,1)
histogram([STATS.cm.rms_rl_m3],'BinWidth',histPARAMS.cm.rms3.binWidth)
title('-3dB Endpoints'); xlabel('RL [dB re 1 \muPa]'); ylabel('Counts')
box on; grid on; 
subplot(2,3,4)
histogram([STATS.cm.rms_rl_m10],'BinWidth',histPARAMS.cm.rms10.binWidth)
title('-10dB Endpoints'); xlabel('RL [dB re 1 \muPa]'); ylabel('Counts')
box on; grid on; 
subplot(2,3,2)
histogram([STATS.cm.rms_rl_e90],'BinWidth',histPARAMS.cm.rms90.binWidth)
title('90% Total Energy Endpoints'); xlabel('RL [dB re 1 \muPa]'); ylabel('Counts')
box on; grid on; 
subplot(2,3,5)
histogram([STATS.cm.rms_rl_e97],'BinWidth',histPARAMS.cm.rms97.binWidth)
title('97% Total Energy Endpoints'); xlabel('RL [dB re 1 \muPa]'); ylabel('Counts')
box on; grid on; 
subplot(2,3,3)
histogram([STATS.cm.rms_rl_gplSE],'BinWidth',histPARAMS.cm.rmsGPL.binWidth)
title('GPL-marked Endpoints'); xlabel('RL [dB re 1 \muPa]'); ylabel('Counts')
box on; grid on; v=v+1;
