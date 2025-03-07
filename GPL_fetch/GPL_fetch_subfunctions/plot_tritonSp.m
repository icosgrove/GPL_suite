% plot_tritonSp.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script contains the plotting code for the spectrogram repurposed
% from Triton.
% Code is sourced from plot_specgram.m and mkspecgram.m from Triton version
% 1.95.20230315. All credit goes to the Scripps Whale Acoustics Laboratory.
% I am repurposing the code with minor changes.
% Written: Ian 08/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

handles.sp = image(time,freq,fetched_call.spectrogram); % Spectrogram
set(handles.sp,'YData',[parm.fp1.freqDelimitLo,parm.fp1.freqDelimitHi])
set(handles.sp,'XData',[0,(length(baseline0)*parm.fftOverlap/parm.SampleFreq)])
set(gca,'TickDir','in'); set(gca,'Layer','Top')
colormap(parm.fp1.colormap) % Colormap 
if parm.fp1.colorbar.switch == 1 % Colorbar
    parm.fp1.colorbar.cb = colorbar('vert');
    yl = get(parm.fp1.colorbar.cb,'YLabel');
    set(yl,'String',parm.fp1.colorbar.YLabelString);
    minc = parm.fp1.colorbar.sp_context_min;
    maxc = parm.fp1.colorbar.sp_context_max;
    set(parm.fp1.colorbar.cb,'YLim',[minc maxc])
    cbYtick = get(parm.fp1.colorbar.cb,'Ytick');  % current colorbar ticks
    newYtick = (cbYtick - parm.fp1.brightness) / (parm.fp1.contrast/100);  % reverse of bright/constrast
    set(parm.fp1.colorbar.cb,'YtickLabel',num2str(newYtick'));  % new labels
end
if parm.fp1.sp_det_boxes == 1
   if cross_flag == 0 
       [rectHandle] = det_boxes(start_sec, finish_sec, mstart_sec, mend_sec, parm, fetched_call, 0, 0, fig);
   else
       [rectHandle] = det_boxes(start_sec, finish_sec, mstart_sec, mend_sec, parm, fetched_call, cross_flag, 0, fig);
   end 
else
    [rectHandle] = det_boxes([], [], mstart_sec, mend_sec, parm, fetched_call, 0, 0, fig);
end
title('Spectrogram [Triton]') 
xlabel('Time [s]'); ylabel('Frequency [Hz]')
text('Units', 'normalized', 'Position', [3/parm.window_length, -0.18, 0], ... % Window start time text
    'String', sprintf(datestr(window_jst+datenum([0 0 0 0 0 cross_flag*parm.window_length]))), 'HorizontalAlignment', 'right','FontSize',10);
box on; set(gca,'Layer','top')
BC_text(parm,1);
if parm.fp1.tf.switch == 1
    text('Units', 'normalized', 'Position', [10/parm.window_length, -0.18, 0], ... 
        'String',strcat('TF File: ',parm.fp1.tf.file),'HorizontalAlignment','left','Interpreter','none')
end
axis([0 (time_range(end)) parm.fp1.freqDelimitLo parm.fp1.freqDelimitHi])