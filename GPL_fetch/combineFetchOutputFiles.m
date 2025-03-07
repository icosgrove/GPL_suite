% Combine Fetch output files together 
% Written: Ian 09/2024

clear; clc

% Load Files
f1 = load('bb_90_fetch_file#2.mat');
f2 = load('bb_90_fetch_file#3.mat');
% f3 = load('bp_30_fetch_file#3.mat');
% f4 = load('bp_30_fetch_file#4.mat');

% Output name
output_name = 'all_bb_90_fetch.mat';

% Remove empty fields
files = [f1 f2];% f3 f4];
for n = 1:length(files)
    if isfield(files(n).hyd_fetch.calls,'spectrogram')
        files(n).hyd_fetch.calls = rmfield(files(n).hyd_fetch.calls,'spectrogram');
    end
    if isfield(files(n).hyd_fetch.calls,'gpl_contour_window')
        files(n).hyd_fetch.calls = rmfield(files(n).hyd_fetch.calls,'gpl_contour_window');
    end
    if isfield(files(n).hyd_fetch.calls,'cm_max_window')
        files(n).hyd_fetch.calls = rmfield(files(n).hyd_fetch.calls,'cm_max_window');
    end
end
f1 = files(1); f2 = files(2);% f3 = files(3); f4 = files(4);

% Add padding for missed win start time
% for n = 1:length(f1.hyd_fetch.window_flag)
%     f1.hyd_fetch.window_flag(n).win_start_time = [];
% end
        
% concatenate
hyd_fetch.calls = [f1.hyd_fetch.calls, f2.hyd_fetch.calls];%, f3.hyd_fetch.calls, f4.hyd_fetch.calls];
hyd_fetch.adhoc_detections = [f1.hyd_fetch.adhoc_detections; f2.hyd_fetch.adhoc_detections];% f3.hyd_fetch.adhoc_detections; f4.hyd_fetch.adhoc_detections];
hyd_fetch.window_flag = [f1.hyd_fetch.window_flag; f2.hyd_fetch.window_flag];% f3.hyd_fetch.window_flag; f4.hyd_fetch.window_flag];
hyd_fetch.parm = f1.hyd_fetch.parm;
hyd_fetch.manual_log = f1.hyd_fetch.manual_log;

save(output_name,'hyd_fetch')