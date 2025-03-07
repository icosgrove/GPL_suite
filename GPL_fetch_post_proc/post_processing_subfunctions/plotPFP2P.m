function [v] = plotPFP2P(STATS,v,histPARAMS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the peak frequency and P2P RL for cm & cm_max contours
% Written: Ian 09/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Peak frequency (cm & cm_max)
figure(v); 
sgtitle('Peak Frequency & P2P RL (Cm & Cm_max)','Interpreter','none')
subplot(2,2,1)
histogram(STATS.cm.peak_freq,'BinWidth',histPARAMS.cm.pf.binWidth)
% [uniq_vals,~,index] = unique(peak_freq);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
title('Peak Frequency: Cm'); xlabel('Frequency [Hz]'); ylabel('Counts')
box on; grid on; 

subplot(2,2,2)
histogram(STATS.cm_max.peak_freq,'BinWidth',histPARAMS.cm_max.pf.binWidth)
% [uniq_vals,~,index] = unique(peak_freq_max);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
title('Peak Frequency: Cm_max','Interpreter','none'); xlabel('Frequency [Hz]'); ylabel('Counts')
box on; grid on;

% P2P RL (cm & cm_max)
subplot(2,2,3)
histogram(STATS.cm.p2p_rl,'BinWidth',histPARAMS.cm.p2pRL.binWidth)
% [uniq_vals,~,index] = unique(p2p_rl);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
title('Peak-to-peak Received Level: Cm'); xlabel('Received Level [dB re 1 \muPa]'); ylabel('Counts')
box on; grid on;

subplot(2,2,4)
histogram(STATS.cm_max.p2p_rl,'BinWidth',histPARAMS.cm_max.p2pRL.binWidth)
% [uniq_vals,~,index] = unique(p2p_rl_max);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
title('Peak-to-peak Received Level: Cm_max','Interpreter','none'); xlabel('Received Level [dB re 1 \muPa]'); ylabel('Counts')
box on; grid on; v=v+1;
