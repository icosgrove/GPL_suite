function multipleOut = combineMultipleGPL(multiple_calls,parm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will take calls that were split into multiple GPL
% detections and combine them into one contour with new associated
% measurements.
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Append detections together and do measurements
for n = 1:length(multiple_calls) % Loop over multiple cases
    select_calls = multiple_calls(n);
    combined_calls.manual_start_time = select_calls.manual_start_time;
    combined_calls.manual_sample_start = select_calls.manual_sample_start;
    combined_calls.manual_sample_end = select_calls.manual_sample_end;
    if parm.fp1.spectrogram == 1
        combined_calls.spectrogram = select_calls.spectrogram;
    end
    if parm.fp1.whitened_sp == 1
        combined_calls.whitened_sp = select_calls.whitened_sp;
    end
    if parm.fp1.gpl_contour_window == 1
        combined_calls.gpl_contour_window = select_calls.gpl_contour_window;
    end
    combined_calls.window_start_time = multiple_calls.window_start_time;
    num_span = length(select_calls.gpl_match);
    for m = 1:num_span - 1 % Loop over span
        temp_calls = select_calls.gpl_match(1:2);
        [combined] = combineMatches(temp_calls,parm,n);
        if m == num_span - 1
            clear temp_calls
            temp_calls = combined;
            break
        end
        if length(select_calls.gpl_match) > 2 % Replace 2 solo with combined and move on to combine with 3+
            select_calls.gpl_match(m+1) = combined;
            select_calls.gpl_match(1:m) = [];
        end
    end
    combined_calls.gpl_match = temp_calls;
    combined_calls.note = '';
    multipleOut(n) = combined_calls;
    if isempty(multipleOut(n).gpl_match.cm.index)
        multipleOut(n).note = 'Empty Contour';
    end
    clear select_calls temp_calls
end
