%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% OUTDATED, use easyFetchCallAnalysis.m
%
%
% GPL Fetch Process 1: Manual-GPL Detection Matching Analysis. Load in
% matched calls and select what statistics or measurements you want. Other
% PR information may be obtained if a full GPL run of the data with the
% same params has been done. Some reprocessing is done. Use
% easyFetchCallAnalysis.m for more in-deptch statistical analysis.
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load FP1 output file 
load('hyd_fetch.mat')
parm = hyd_fetch.parm;
calls_append = hyd_fetch.calls;

%% User Input

% If known: Enter total detections for this set
total_gpl_detections = 364176;
total_manual_detections = 547;


%% Analysis

% Process FP1 match cases
note_vals = {calls_append.note};
miss_calls = calls_append(strcmp(note_vals, 'GPL Miss'));
overlap_calls = calls_append(strcmp(note_vals, 'Window Overlap'));
multiple_calls = calls_append(strcmp(note_vals, 'Multiple GPL'));
single_calls = calls_append(strcmp(note_vals, ''));

% Process singles for empty contours
nonempty = removeEmptySingleContour(single_calls);
single_save = nonempty.calls;

% Process Overlaps
overlapOut = windowOverlapCombine(overlap_calls,parm);
overlap_save = overlapOut(strcmp({overlapOut.note}, ''));

% Process Multiple GPL detections
multipleOut = combineMultipleGPL(multiple_calls,parm);
multiple_save = multipleOut(cellfun('isempty', {multipleOut.note}));

% PR Statistics
TP = length(single_calls) + length(overlap_save) + length(multipleOut);
FN = length(miss_calls);
FP = total_gpl_detections - TP;  
R = TP / (TP + FN);
P = TP / (TP + FP);

% Calls with measurements that may be plotted
plot_calls = [single_save overlap_save multiple_save];

%%%%%%%%% Plots

%% Manual-GPL Matches: Cm
% Extract Match Measurements
v=1;
for k = 1:length(plot_calls)
    dur_s(k) = plot_calls(k).gpl_match.cm.duration_sec;
    slope(k) = plot_calls(k).gpl_match.cm.slope;
    band_hz(k) = plot_calls(k).gpl_match.cm.freq_bandwidth_hz;
    start_hz(k) = plot_calls(k).gpl_match.cm.start_freq_hz;
    end_hz(k) = plot_calls(k).gpl_match.cm.end_freq_hz;
    min_hz(k) = plot_calls(k).gpl_match.cm.min_freq_hz;
    max_hz(k) = plot_calls(k).gpl_match.cm.max_freq_hz - 2 ;
    abs_band_hz(k) = plot_calls(k).gpl_match.cm.abs_bandwidth_hz;
    gpl_rl(k) = plot_calls(k).gpl_match.spec_rl;
end

% % Duration (s): Matches only
% figure(v)
% [uniq_vals,~,index] = unique(dur_s);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
% title('Match: Duration'); xlabel('Time [s]'); ylabel('Counts')
% box on; grid on; v=v+1;
% 
% % Slope (Hz/s): Matches only
% figure(v) 
% histogram(slope,'BinWidth',1)
% title('Match: Slope'); xlabel('Slope [Hz/s]'); ylabel('Counts')
% box on; grid on; v=v+1;
% 
% % Nonzero Bandwidth (Hz): Matches only
% figure(v)
% [uniq_vals,~,index] = unique(band_hz);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
% title('Match: Nonzero Frequency Bandwidth'); xlabel('Bandwidth [Hz]'); ylabel('Counts')
% box on; grid on; v=v+1;
% 
% % Full Bandwidth (Hz): Matches only
% figure(v)
% [uniq_vals,~,index] = unique(abs_band_hz);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
% title('Match: Absolute Frequency Bandwidth'); xlabel('Bandwidth [Hz]'); ylabel('Counts')
% box on; grid on; v=v+1;

% % Start Frequency (Hz): Matches only
% figure(v); subplot(2,2,1) 
% histogram(start_hz,'BinWidth',1)
% title('Match: Start Frequency (Cm)'); xlabel('Frequency [Hz]'); ylabel('Counts')
% box on; grid on; v=v+1;
% 
% % Start Frequency (Hz): Matches only
% subplot(2,2,2) 
% histogram(end_hz,'BinWidth',1)
% title('Match: End Frequency (Cm)'); xlabel('Frequency [Hz]'); ylabel('Counts')
% box on; grid on; v=v+1;
% 
% % Min. Frequency (Hz): Matches only
% subplot(2,2,3)
% [uniq_vals,~,index] = unique(min_hz);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
% title('Match: Minimum Frequency (Cm)'); xlabel('Frequency [Hz]'); ylabel('Counts')
% box on; grid on; v=v+1;
% 
% % Max. Frequency (Hz): Matches only
% subplot(2,2,4)
% [uniq_vals,~,index] = unique(max_hz);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
% title('Match: Maximum Frequency (Cm)'); xlabel('Frequency [Hz]'); ylabel('Counts')
% box on; grid on; v=v+1;

% % GPL RMS Recieved Level
% figure(v)
% histogram(gpl_rl,'BinWidth',0.01)
% title('Match: GPL RMS Recieved Level (Cm)'); xlabel('Recieved Level'); ylabel('Counts')
% box on; grid on; v=v+1;



%% Manual-GPL Matches: Cm_max
% Extract Match Measurements
for k = 1:length(plot_calls)
    dur_s_cm_max(k) = plot_calls(k).gpl_match.cm_max.duration_sec;
    slope_cm_max(k) = plot_calls(k).gpl_match.cm_max.slope;
    band_hz_cm_max(k) = plot_calls(k).gpl_match.cm_max.freq_bandwidth_hz;
    start_hz_cm_max(k) = plot_calls(k).gpl_match.cm_max.start_freq_hz;
    end_hz_cm_max(k) = plot_calls(k).gpl_match.cm_max.end_freq_hz;
    min_hz_cm_max(k) = plot_calls(k).gpl_match.cm_max.min_freq_hz;
    max_hz_cm_max(k) = plot_calls(k).gpl_match.cm_max.max_freq_hz;
    abs_band_hz_cm_max(k) = plot_calls(k).gpl_match.cm_max.abs_bandwidth_hz;
end

% % Duration (s): Matches only
% figure(v) 
% [uniq_vals,~,index] = unique(dur_s_cm_max);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
% title('Match: Duration (Cm-max)'); xlabel('Time [s]'); ylabel('Counts')
% box on; grid on; v=v+1;
% 
% % Slope (Hz/s): Matches only
% figure(v) 
% histogram(slope_cm_max,'BinWidth',1)
% title('Match: Slope (Cm-max)'); xlabel('Slope [Hz/s]'); ylabel('Counts')
% box on; grid on; v=v+1;
% 
% % Nonzero Bandwidth (Hz): Matches only
% figure(v)
% [uniq_vals,~,index] = unique(band_hz_cm_max);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
% title('Match: Nonzero Frequency Bandwidth (Cm-max'); xlabel('Bandwidth [Hz]'); ylabel('Counts')
% box on; grid on; v=v+1;
% 
% % Full Bandwidth (Hz): Matches on
% figure(v)
% [uniq_vals,~,index] = unique(abs_band_hz_cm_max);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
% title('Match: Absolute Frequency Bandwidth (Cm-max)'); xlabel('Bandwidth [Hz]'); ylabel('Counts')
% box on; grid on; v=v+1;
% 
% % Start Frequency (Hz): Matches only
% figure(v); subplot(2,2,1)
% histogram(start_hz_cm_max,'BinWidth',1)
% title('Match: Start Frequency (Cm-max)'); xlabel('Frequency [Hz]'); ylabel('Counts')
% box on; grid on; v=v+1;
% 
% % End Frequency (Hz): Matches only
% subplot(2,2,2)
% histogram(end_hz_cm_max,'BinWidth',1)
% title('Match: End Frequency (Cm-max)'); xlabel('Frequency [Hz]'); ylabel('Counts')
% box on; grid on; v=v+1;
% 
% % Min. Frequency (Hz): Matches only
% subplot(2,2,3)
% [uniq_vals,~,index] = unique(min_hz_cm_max);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
% title('Match: Minimum Frequency (Cm-max)'); xlabel('Frequency [Hz]'); ylabel('Counts')
% box on; grid on; v=v+1;
% 
% % Max. Frequency (Hz): Matches only
% subplot(2,2,4)
% [uniq_vals,~,index] = unique(max_hz_cm_max);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
% title('Match: Maximum Frequency (Cm-max)'); xlabel('Frequency [Hz]'); ylabel('Counts')
% box on; grid on; v=v+1;


%% Manual-GPL Match Scatters

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


%% Aggregate GPL detections
% Load detection files for mass measurement plotting
[filename, pathname] = uigetfile('*.mat'); % Select File
cwd = pwd; 
cd(pathname) % Set current directory to path containing xwav file
file_dir = pwd; 
addpath(pwd); 
detfiles = dir('*.mat'); % Variable files contains xwav data
cd(cwd); % Set current directory back to current working directory
s=1; call_total = 0;
for n = 1:length(detfiles)
    fprintf('\nLoading file %i\n',n)
    load(sprintf(detfiles(n).name))
    fprintf('Processing file %i\n',n)
    call_total = call_total + length(hyd.detection.calls);
    for k = 1:length(hyd.detection.calls)
        if ~isempty(hyd.detection.calls(k).cm.index)
            dur_s_total(s) = hyd.detection.calls(k).cm.duration_sec;
            slope_total(s) = hyd.detection.calls(k).cm.slope;
            band_hz_total(s) = hyd.detection.calls(k).cm.freq_bandwidth_hz;
            start_hz_total(s) = hyd.detection.calls(k).cm.start_freq_hz;
            end_hz_total(s) = hyd.detection.calls(k).cm.end_freq_hz;
            min_hz_total(s) = hyd.detection.calls(k).cm.min_freq_hz;
            max_hz_total(s) = hyd.detection.calls(k).cm.max_freq_hz;
            abs_band_hz_total(s) = hyd.detection.calls(k).cm.abs_bandwidth_hz;
            gpl_rl_total(s) = hyd.detection.calls(k).spec_rl;
            s=s+1;
        end
    end
end

% % Duration (s): All GPL detections
% figure(v)
% [uniq_vals,~,index] = unique(dur_s_total);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
% title('All GPL: Duration'); xlabel('Time [s]'); ylabel('Counts')
% box on; grid on; v=v+1;
% 
% % Slope (Hz/s): All GPL detections
% figure(v) 
% histogram(slope_total,'BinWidth',1)
% title('All GPL: Slope'); xlabel('Slope [Hz/s]'); ylabel('Counts')
% box on; grid on; v=v+1;
% 
% % Nonzero Bandwidth (Hz): All GPL detections
% figure(v)
% [uniq_vals,~,index] = unique(band_hz_total);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
% title('All GPL: Nonzero Frequency Bandwidth'); xlabel('Bandwidth [Hz]'); ylabel('Counts')
% box on; grid on; v=v+1;
% 
% % Full Bandwidth (Hz): All GPL detections
% figure(v)
% [uniq_vals,~,index] = unique(abs_band_hz_total);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
% title('All GPL: Absolute Frequency Bandwidth'); xlabel('Bandwidth [Hz]'); ylabel('Counts')
% box on; grid on; v=v+1;
% 
% % Start Frequency (Hz): All GPL detections
% figure(v)
% subplot(2,2,1)
% histogram(start_hz_total,'BinWidth',1)
% title('All GPL: Start Frequency'); xlabel('Frequency [Hz]'); ylabel('Counts')
% box on; grid on; v=v+1;
% 
% % Start Frequency (Hz): All GPL detections
% subplot(2,2,2) 
% histogram(end_hz_total,'BinWidth',1)
% title('All GPL: End Frequency'); xlabel('Frequency [Hz]'); ylabel('Counts')
% box on; grid on; v=v+1;
% 
% % Min. Frequency (Hz): All GPL detections
% subplot(2,2,3) 
% [uniq_vals,~,index] = unique(min_hz_total);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
% title('All GPL: Minimum Frequency'); xlabel('Frequency [Hz]'); ylabel('Counts')
% box on; grid on; v=v+1;
% 
% % Max. Frequency (Hz): All GPL detections
% subplot(2,2,4)
% [uniq_vals,~,index] = unique(max_hz_total);
% counts = accumarray(index,1);
% bar(uniq_vals,counts);
% title('All GPL: Maximum Frequency'); xlabel('Frequency [Hz]'); ylabel('Counts')
% box on; grid on; v=v+1;

% GPL RMS Recieved Level
figure(v)
histogram(gpl_rl_total,'BinWidth',0.04)
title('All GPL RMS Recieved Level'); xlabel('Recieved Level'); ylabel('Counts')
box on; grid on; v=v+1;


%% Scatters: Matches placed within all GPL detections

% % Duration vs. RL Scatter (Cm & Cm_max)
% figure(v); hold on
% scatter(gpl_rl_total,dur_s_total,10,'b','MarkerFaceColor','b','MarkerEdgeColor','b')
% scatter(gpl_rl,dur_s,10,'r','MarkerFaceColor','r')
% title('All GPL & Cm: RL vs. Duration'); ylabel('Duration [s]'); xlabel('Received Level')
% box on; grid on; v=v+1;
% legend('All GPL','Match: Cm')
% hold off
% 
% Slope vs. RL Scatter (cm & cm_max)
% figure(v); hold on
% scatter(gpl_rl_total,slope_total,8,'b','MarkerFaceColor','b','MarkerEdgeColor','b')
% scatter(gpl_rl,slope,8,'r','MarkerFaceColor','r')
% title('All GPL & Cm: RL vs. Slope'); ylabel('Slope [Hz/s]'); xlabel('Received Level')
% box on; grid on; v=v+1;
% legend('All GPL','Match: Cm')
% hold off
% 
% Frequency Nonzero Bandwidth vs. RL Scatter (cm & cm_max)
% figure(v); hold on
% scatter(gpl_rl_total,band_hz_total,10,'b','MarkerFaceColor','b','MarkerEdgeColor','b')
% scatter(gpl_rl,band_hz,10,'r','MarkerFaceColor','r')
% title('All GPL & Cm: RL vs. Nonzero Bandwidth'); ylabel('Bandwidth [Hz]'); xlabel('Received Level')
% box on; grid on; v=v+1;
% legend('All GPL','Cm')
% hold off
% 
% % Frequency Absolute Bandwidth vs. RL Scatter (Cm & Cm_max)
% figure(v); hold on
% scatter(gpl_rl_total,abs_band_hz_total,10,'b','MarkerFaceColor','b','MarkerEdgeColor','b')
% scatter(gpl_rl,abs_band_hz,10,'r','MarkerFaceColor','r')
% title('All GPL & Cm: RL vs. Absolute Bandwidth'); ylabel('Bandwidth [Hz]'); xlabel('Received Level')
% box on; grid on; v=v+1;
% legend('All GPL','Cm')
% hold off
% 
% Frequency Start/End/Min/Max vs. RL Scatter (Cm & Cm_max)
figure(v); 
subplot(2,2,1); hold on % Start 
scatter(gpl_rl_total,start_hz_total,10,'b','MarkerFaceColor','b','MarkerEdgeColor','b')
scatter(gpl_rl,start_hz,10,'r','MarkerFaceColor','r')
title('All GPL & Cm: RL vs. Start Frequency'); xlabel('Received Level'); ylabel('Start Frequency [Hz]')
box on; grid on; legend('GPL All','Cm'); hold off

subplot(2,2,2); hold on % End
scatter(gpl_rl_total,end_hz_total,10,'b','MarkerFaceColor','b','MarkerEdgeColor','b')
scatter(gpl_rl,end_hz,10,'r','MarkerFaceColor','r')
title('All GPL & Cm: RL vs. End Frequency'); xlabel('Received Level'); ylabel('End Frequency [Hz]')
box on; grid on; legend('GPL All','Cm'); hold off

subplot(2,2,3); hold on % Min
scatter(gpl_rl_total,min_hz_total,10,'b','MarkerFaceColor','b','MarkerEdgeColor','b')
scatter(gpl_rl,min_hz,10,'r','MarkerFaceColor','r')
title('All GPL & Cm: RL vs. Minimum Frequency'); xlabel('Received Level'); ylabel('Minimum Frequency [Hz]')
box on; grid on; legend('GPL All','Cm'); hold off

subplot(2,2,4); hold on % Max
scatter(gpl_rl_total,max_hz_total,10,'b','MarkerFaceColor','b','MarkerEdgeColor','b')
scatter(gpl_rl,max_hz,10,'r','MarkerFaceColor','r')
title('All GPL & Cm: RL vs. Maximum Frequency'); xlabel('Received Level'); ylabel('Maximum Frequency [Hz]')
box on; grid on; legend('GPL All','Cm'); hold off












