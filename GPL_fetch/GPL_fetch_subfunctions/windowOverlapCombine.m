function overlapOut = windowOverlapCombine(overlap_calls,parm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will process window overlap calls from GPL Fetch process 1.
% Overlap cases will be checked if they are missed partially or totally. If
% the two matches are valid, they will be combined together.
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = 1;
for n = 1:2:length(overlap_calls)  
    if isempty(overlap_calls(n).gpl_match) % Front Miss
        if isempty(overlap_calls(n+1).gpl_match) % Back Miss
            overlapOut(k) = overlap_calls(n); 
            overlapOut(k).note = 'GPL Miss'; 
            if isempty(overlapOut(k).gpl_match.cm.index)
                overlapOut(k).note = 'Empty Contour';
            else
                overlapOut(k).note = '';
            end
            k=k+1;
        else % Front Miss, Back Pass
            overlapOut(k) = overlap_calls(n+1); 
            overlapOut(k).note = ''; 
            if isempty(overlapOut(k).gpl_match.cm.index)
                overlapOut(k).note = 'Empty Contour';
            else
                overlapOut(k).note = '';
            end
            k=k+1;
        end
    else % Front Pass
        if isempty(overlap_calls(n+1).gpl_match) % Back Miss
            overlapOut(k) = overlap_calls(n); 
            overlapOut(k).note = ''; 
            if isempty(overlapOut(k).gpl_match.cm.index)
                overlapOut(k).note = 'Empty Contour';
            else
                overlapOut(k).note = '';
            end
            k=k+1;
        else % Front Pass, Back Pass: Combine contours
            overlapOut(k).manual_start_time = overlap_calls(n).manual_start_time;
            overlapOut(k).manual_sample_start = overlap_calls(n).manual_sample_start;
            overlapOut(k).manual_sample_end = overlap_calls(n+1).manual_sample_end;
            if parm.fp1.spectrogram == 1
                overlapOut(k).spectrogram = [overlap_calls(n).spectrogram overlap_calls(n+1).spectrogram];
            end
            if parm.fp1.whitened_sp == 1
                overlapOut(k).whitened_sp = [overlap_calls(n).whitened_sp overlap_calls(n+1).whitened_sp];
            end
            if parm.fp1.gpl_contour_window == 1
                overlapOut(k).gpl_contour_window = [overlap_calls(n).gpl_contour_window overlap_calls(n+1).gpl_contour_window];
            end
            overlapOut(k).window_start_time = {overlap_calls(n).window_start_time overlap_calls(n+1).window_start_time};
            overlapOut(k).gpl_match.start_time = overlap_calls(n).gpl_match.start_time;
            overlapOut(k).gpl_match.end_time = overlap_calls(n).gpl_match.start_time + overlap_calls(n+1).gpl_match.end_time;
            overlapOut(k).gpl_match.cm = combineContours(overlap_calls(n).gpl_match.cm,overlap_calls(n+1).gpl_match.cm,parm);
            overlapOut(k).gpl_match.cm_max = combineContours(overlap_calls(n).gpl_match.cm_max,overlap_calls(n+1).gpl_match.cm_max,parm);
            overlapOut(k).gpl_match.cm_max2 = combineContours(overlap_calls(n).gpl_match.cm_max2,overlap_calls(n+1).gpl_match.cm_max2,parm);
            overlapOut(k).gpl_match.spec_noise = (overlap_calls(n).gpl_match.spec_noise + overlap_calls(n+1).gpl_match.spec_noise) / 2;
            overlapOut(k).gpl_match.spec_rl = (overlap_calls(n).gpl_match.spec_rl + overlap_calls(n+1).gpl_match.spec_rl) / 2;
            overlapOut(k).gpl_match.julian_start_time = overlap_calls(n).gpl_match.julian_start_time;
            overlapOut(k).gpl_match.julian_end_time = overlap_calls(n+1).gpl_match.julian_end_time;
            if isempty(overlapOut(k).gpl_match.cm.index)
                overlapOut(k).note = 'Empty Contour';
            else
                overlapOut(k).note = '';
            end
            k=k+1;
        end
    end
end
            