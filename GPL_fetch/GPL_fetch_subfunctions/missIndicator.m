function missIndicator(fig,mstart_sec,mend_sec,fetched_call,parm,mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will create and indicator for a GPL miss by putting a slash
% through the manual call
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(fig)
subplot(3,1,3); hold on
if mode == 1 % Miss orange x
    plot3([mstart_sec mend_sec],[parm.freq_hi parm.freq_lo],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-','Color',[245 138 56]./255,'LineWidth',4)
    plot3([mstart_sec mend_sec],[parm.freq_lo parm.freq_hi],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-','Color',[245 138 56]./255,'LineWidth',4)
elseif mode == 2 % False entry: Red x
    plot3([mstart_sec mend_sec],[parm.freq_hi parm.freq_lo],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-','Color',[255 0 0]./255,'LineWidth',4)
    plot3([mstart_sec mend_sec],[parm.freq_lo parm.freq_hi],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-','Color',[255 0 0]./255,'LineWidth',4)
elseif mode == 3 % Swap bluish x
    plot3([mstart_sec mend_sec],[parm.freq_hi parm.freq_lo],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-','Color',[4 204 197]./255,'LineWidth',4)
    plot3([mstart_sec mend_sec],[parm.freq_lo parm.freq_hi],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-','Color',[4 204 197]./255,'LineWidth',4)
end
hold off