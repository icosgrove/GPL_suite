function [call_ready, window_ready] = prepareWindowOutput(calls_TEMP,window_TEMP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will combine together relevant fields for window crossing
% cases so that it is just one element in the output structure, instead of
% two. 
% Written: Ian 08/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

window_ready = [];

call_ready.manual_start_time = calls_TEMP(1).manual_start_time;
call_ready.manual_sample_start = calls_TEMP(1).manual_sample_start;
call_ready.manual_sample_end = calls_TEMP(1).manual_sample_end;
if isfield(calls_TEMP(1),'spectrogram')
    call_ready.spectrogram = {calls_TEMP(1).spectrogram calls_TEMP(2).spectrogram};
end
if isfield(calls_TEMP(2),'gpl_contour_window')
    call_ready.gpl_contour_window = {calls_TEMP(1).gpl_contour_window calls_TEMP(2).gpl_contour_window};
end
if isfield(calls_TEMP(2),'cm_max_window')
    call_ready.cm_max_window = {calls_TEMP(1).cm_max_window calls_TEMP(2).cm_max_window};
end
call_ready.window_start_time = {calls_TEMP(1).window_start_time calls_TEMP(2).window_start_time};
call_ready.gpl_match = [calls_TEMP(1).gpl_match, calls_TEMP(2).gpl_match];
call_ready.note = {calls_TEMP(1).note calls_TEMP(2).note 'Window Overlap'};
if ~isempty(window_TEMP)
    for n = 1:length(window_TEMP)
        window_ready = [window_ready; window_TEMP(n)];
    end
end