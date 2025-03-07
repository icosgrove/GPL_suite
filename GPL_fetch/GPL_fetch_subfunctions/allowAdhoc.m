function [allowAdhocFlag] = allowAdhoc(mst_julian,met_julian,sample_range,SampleRangeStart,SampleRangeEnd,xwav_start,parm,n,file_ident)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will identify how many calls are flagged in a specific
% window. This is done so that later, when plotting, the user will be
% notified that other detections are in the window and to wait for Adhoc
% detections.
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Don't do deck test adjustment beyond disk 1
if file_ident > 1
    parm.pre_offset = 1;
end

% Sort the next detection into a window and if it matches the current
% window set flag to have adhoc wait
if n ~= length(mst_julian)
    
    % Convert next call into samples
    next_call_sampStart = round(seconds(datetime(mst_julian(n+1),'ConvertFrom','datenum')-datetime(xwav_start(file_ident),'ConvertFrom','datenum'))*parm.SampleFreq)+(parm.pre_offset-1)*parm.nrec;
    next_call_sampEnd = round(seconds(datetime(met_julian(n+1),'ConvertFrom','datenum')-datetime(xwav_start(file_ident),'ConvertFrom','datenum'))*parm.SampleFreq)+(parm.pre_offset-1)*parm.nrec;

    % Sort and check if next call is in the same window
    startSort = find(sort([sample_range next_call_sampStart]) == next_call_sampStart);
    endSort = find(sort([sample_range next_call_sampEnd]) == next_call_sampEnd);
    if startSort ~= 1 % Adjust back to actual window
        startSort = startSort - 1;
        endSort = endSort - 1;
    end
    if startSort == SampleRangeStart % Same start window
        if endSort == SampleRangeEnd % Same end window
            allowAdhocFlag = 0; % Another detection is in the same window
        else % Implies the end is in the next window, so next call is window overlap
            allowAdhocFlag = 1; % Allow adhoc since window overlap will auto-move forward
        end
    else % Different Start window
        allowAdhocFlag = 1;
    end
else % End of manual calls
    allowAdhocFlag = 1;
end
