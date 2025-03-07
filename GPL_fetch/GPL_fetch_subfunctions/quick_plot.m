% Plot
n=34; n1 = 1;
mat = zeros(hyd_fetch.calls(n).gpl_match(n1).cm.size(1), hyd_fetch.calls(n).gpl_match(n1).cm.size(2));
mat(hyd_fetch.calls(n).gpl_match(n1).cm.index) = hyd_fetch.calls(n).gpl_match(n1).cm.values;
figure(8); imagesc(mat); set(gca,'YDir','normal')
% Timestamps
manual_start_time = hyd_fetch.calls(n).manual_start_time % Manual pick
paired_gpl_start_time = datestr(hyd_fetch.calls(n).gpl_match(n1).julian_start_time,'yyyy-mm-dd HH:MM:SS.FFF') % GPL match

% Adhoc
k = 4;
mat = zeros(hyd_fetch.adhoc_detections(k).cm.size(1),hyd_fetch.adhoc_detections(k).cm.size(2));
mat(hyd_fetch.adhoc_detections(k).cm.index) = hyd_fetch.adhoc_detections(k).cm.values;
figure(9); imagesc(mat); set(gca,'YDir','normal')
adhoc_start_time = datestr(hyd_fetch.adhoc_detections(k).julian_start_time,'yyyy-mm-dd HH:MM:SS.FFF')