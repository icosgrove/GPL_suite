% Old scatter plots for displaying Fetch measurements
% Used by Ian 07-09/2024

% % Manual-GPL Match Scatters
% 
% % Duration vs. RL Scatter (Cm & Cm_max)
% figure(v); hold on
% scatter(gpl_rl,dur_s,10,'r','MarkerFaceColor','r')
% scatter(gpl_rl,dur_s_cm_max,10,'r','MarkerFaceColor','b','MarkerEdgeColor','b')
% title('Matches: RL vs. Duration (cm & cm-max)'); ylabel('Duration [s]'); xlabel('Received Level')
% box on; grid on; v=v+1;
% legend('Cm','Cm-max')
% hold off
% 
% % Slope vs. RL Scatter (cm & cm_max)
% figure(v); hold on
% scatter(gpl_rl,slope,8,'r','MarkerFaceColor','r')
% scatter(gpl_rl,slope_cm_max,8,'r','MarkerFaceColor','b','MarkerEdgeColor','b')
% title('Matches: RL vs. Slope (cm & cm-max)'); ylabel('Slope [Hz/s]'); xlabel('Received Level')
% box on; grid on; v=v+1;
% legend('Cm','Cm-max')
% hold off
% 
% % Frequency Nonzero Bandwidth vs. RL Scatter (cm & cm_max)
% figure(v); hold on
% scatter(gpl_rl,band_hz,10,'r','MarkerFaceColor','r')
% scatter(gpl_rl,band_hz_cm_max,10,'r','MarkerFaceColor','b','MarkerEdgeColor','b')
% title('Matches: RL vs. Nonzero Bandwidth (cm & cm-max)'); ylabel('Bandwidth [Hz]'); xlabel('Received Level')
% box on; grid on; v=v+1;
% legend('Cm','Cm-max')
% hold off
% 
% % Frequency Absolute Bandwidth vs. RL Scatter (Cm & Cm_max)
% figure(v); hold on
% scatter(gpl_rl,abs_band_hz,10,'r','MarkerFaceColor','r')
% scatter(gpl_rl,abs_band_hz_cm_max,10,'r','MarkerFaceColor','b','MarkerEdgeColor','b')
% title('Matches: RL vs. Absolute Bandwidth (cm & cm-max)'); ylabel('Bandwidth [Hz]'); xlabel('Received Level')
% box on; grid on; v=v+1;
% legend('Cm','Cm-max')
% hold off
% 
% % Frequency Start/End/Min/Max vs. RL Scatter (Cm & Cm_max)
% figure(v); 
% subplot(2,2,1); hold on % Start 
% scatter(gpl_rl,start_hz,10,'r','MarkerFaceColor','r')
% scatter(gpl_rl,start_hz_cm_max,10,'r','MarkerFaceColor','b','MarkerEdgeColor','b')
% title('Matches: RL vs. Start Frequency (cm & cm-max)'); xlabel('Received Level'); ylabel('Start Frequency [Hz]')
% box on; grid on; legend('Cm','Cm-max'); hold off
% 
% subplot(2,2,2); hold on % End
% scatter(gpl_rl,end_hz,10,'r','MarkerFaceColor','r')
% scatter(gpl_rl,end_hz_cm_max,10,'r','MarkerFaceColor','b','MarkerEdgeColor','b')
% title('Matches: RL vs. End Frequency (cm & cm-max)'); xlabel('Received Level'); ylabel('End Frequency [Hz]')
% box on; grid on; legend('Cm','Cm-max'); hold off
% 
% subplot(2,2,3); hold on % Min
% scatter(gpl_rl,min_hz,10,'r','MarkerFaceColor','r')
% scatter(gpl_rl,min_hz_cm_max,10,'r','MarkerFaceColor','b','MarkerEdgeColor','b')
% title('Matches: RL vs. Minimum Frequency (cm & cm-max)'); xlabel('Received Level'); ylabel('Minimum Frequency [Hz]')
% box on; grid on; legend('Cm','Cm-max'); hold off
% 
% subplot(2,2,4); hold on % Max
% scatter(gpl_rl,max_hz,10,'r','MarkerFaceColor','r')
% scatter(gpl_rl,max_hz_cm_max,10,'b','MarkerFaceColor','b','MarkerEdgeColor','b')
% title('Matches: RL vs. Maximum Frequency (cm & cm-max)'); xlabel('Received Level'); ylabel('Maximum Frequency [Hz]')
% box on; grid on; legend('Cm','Cm-max'); hold off














