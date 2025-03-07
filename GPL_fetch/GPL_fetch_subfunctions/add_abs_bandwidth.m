% % Add absolute bandwidth
% for m = 1:length(hyd_fetch.calls)
%     for n = 1:length(hyd_fetch.calls(m).gpl_match)
%         if ~isempty([hyd_fetch.calls(m).gpl_match(n)])
%             if ~isempty([hyd_fetch.calls(m).gpl_match(n).cm.index])
%                 hyd_fetch.calls(m).gpl_match(n).cm.abs_bandwidth_hz = hyd_fetch.calls(m).gpl_match(n).cm.max_freq_hz - hyd_fetch.calls(m).gpl_match(n).cm.min_freq_hz + 1;
%             end
%             if ~isempty([hyd_fetch.calls(m).gpl_match(n).cm_max.index])
%                 hyd_fetch.calls(m).gpl_match(n).cm_max.abs_bandwidth_hz = hyd_fetch.calls(m).gpl_match(n).cm_max.max_freq_hz - hyd_fetch.calls(m).gpl_match(n).cm_max.min_freq_hz + 1;
%             end
%             if ~isempty([hyd_fetch.calls(m).gpl_match(n).cm_max2.index])
%                 hyd_fetch.calls(m).gpl_match(n).cm_max2.abs_bandwidth_hz = hyd_fetch.calls(m).gpl_match(n).cm_max2.max_freq_hz - hyd_fetch.calls(m).gpl_match(n).cm_max2.min_freq_hz + 1;
%             end
%         end
%     end
% end

% Subtract 2 Hz from start/end freq
for m = 1:length(h yd_fetch.calls)
    for n = 1:length(hyd_fetch.calls(m).gpl_match)
        if ~isempty([hyd_fetch.calls(m).gpl_match(n)])
            if ~isempty([hyd_fetch.calls(m).gpl_match(n).cm.index])
                hyd_fetch.calls(m).gpl_match(n).cm.start_freq_hz = hyd_fetch.calls(m).gpl_match(n).cm.start_freq_hz-2;
                hyd_fetch.calls(m).gpl_match(n).cm.end_freq_hz = hyd_fetch.calls(m).gpl_match(n).cm.end_freq_hz-2;
            end
            if ~isempty([hyd_fetch.calls(m).gpl_match(n).cm_max.index])
                hyd_fetch.calls(m).gpl_match(n).cm_max.start_freq_hz = hyd_fetch.calls(m).gpl_match(n).cm_max.start_freq_hz-2;
                hyd_fetch.calls(m).gpl_match(n).cm_max.end_freq_hz = hyd_fetch.calls(m).gpl_match(n).cm_max.end_freq_hz-2; 
            end
            if ~isempty([hyd_fetch.calls(m).gpl_match(n).cm_max2.index])
                hyd_fetch.calls(m).gpl_match(n).cm_max2.start_freq_hz = hyd_fetch.calls(m).gpl_match(n).cm_max2.start_freq_hz-2; 
                hyd_fetch.calls(m).gpl_match(n).cm_max2.end_freq_hz = hyd_fetch.calls(m).gpl_match(n).cm_max2.end_freq_hz-2;
            end
        end
    end
end