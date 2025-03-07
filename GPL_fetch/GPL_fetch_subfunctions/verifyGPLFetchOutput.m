function verifyGPLFetchOutput(hyd_fetch,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is intended to provide the user with information to verify
% that the outputs of GPL fetch make sense and line up with what is
% expected. 
% Written: Ian 08/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table = hyd_fetch.manual_log;
parm = hyd_fetch.parm;

if ~isempty(hyd_fetch.calls(n).gpl_match)
    out_mst = hyd_fetch.calls(n).manual_start_time;
    out_win = hyd_fetch.calls.window_start_time;
    if size(hyd_fetch.calls(n).gpl_match,1) > 1
        return
    end
    out_cm = zeros(hyd_fetch.calls(n).gpl_match.cm.size(1),hyd_fetch.calls(n).gpl_match.cm.size(2));
    out_cm(hyd_fetch.calls(n).gpl_match.cm.index) = hyd_fetch.calls(n).gpl_match.cm.values;
    
    freq_range = (parm.freq_lo:parm.SampleFreq/parm.fftl:parm.freq_hi);
    time_range = (1:hyd_fetch.calls(n).gpl_match.cm.size(2))*(parm.fftOverlap/parm.SampleFreq);
    slope = hyd_fetch.calls(n).gpl_match.cm.slope.*time_range + hyd_fetch.calls(n).gpl_match.cm.slope_y_int - 1;
    
    figure(9); hold on
    ax = gca;
    ax.OuterPosition = [0 0.3 1 0.7];
    imagesc(time_range,freq_range,out_cm)
    plot(time_range,slope,'w--','LineWidth',1.5)
    xlabel('Time [s]'); ylabel('Frequency [Hz]')
    % colormap('jet')
    title(sprintf('Pairing %i',n))
    axis tight
    text('Units','normalized','Position',[0 -.15 0],'String',sprintf('Window: %s',out_win))
    text('Units','normalized','Position',[0 -.25 0],'String',sprintf('Table ST: %s',datestr(table.StartTime(n),'yyyy-mm-dd HH:MM:SS.FFF')))
    text('Units','normalized','Position',[0 -.35 0],'String',sprintf('Manual ST: %s',out_mst))
    text('Units','normalized','Position',[0,-.45,0],'String',sprintf('GPL ST: %s',datestr(hyd_fetch.calls(n).gpl_match.julian_start_time,'yyyy-mm-dd HH:MM:SS.FFF')))
    text('Units','normalized','Position',[1 -.15 0],'String',sprintf('Duration: %.1fs',hyd_fetch.calls(n).gpl_match.cm.duration_sec),'HorizontalAlignment','right')
    text('Units','normalized','Position',[1 -.25 0],'String',sprintf('Slope: %.4fHz/s',hyd_fetch.calls(n).gpl_match.cm.slope),'HorizontalAlignment','right')
    text('Units','normalized','Position',[1 -.35 0],'String',sprintf('Nonz. Freq. Band.: %.1fHz',hyd_fetch.calls(n).gpl_match.cm.freq_bandwidth_hz),'HorizontalAlignment','right')
    text('Units','normalized','Position',[1 -.45 0],'String',sprintf('Abs. Freq. Band.: %.1fHz',hyd_fetch.calls(n).gpl_match.cm.abs_bandwidth_hz),'HorizontalAlignment','right')
    text('Units','normalized','Position',[1 -.55 0],'String',sprintf('Start/End Freq.: %.1f,%.1fHz',hyd_fetch.calls(n).gpl_match.cm.start_freq_hz,hyd_fetch.calls(n).gpl_match.cm.end_freq_hz),'HorizontalAlignment','right')
    text('Units','normalized','Position',[1 -.65 0],'String',sprintf('Min/Max Freq.: %.1f,%.1fHz',hyd_fetch.calls(n).gpl_match.cm.min_freq_hz,hyd_fetch.calls(n).gpl_match.cm.max_freq_hz),'HorizontalAlignment','right')
    disp('Waiting for keystroke..')
    pause
    clf
else
    disp('No detection paired to display.')
end