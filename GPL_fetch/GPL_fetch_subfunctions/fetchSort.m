function [SampleRangeStart, SampleRangeEnd, error_flag, file_string, sample_range, detStart, detEnd, detStart_rel, detEnd_rel, file_ident] = fetchSort(mst_julian, met_julian, n, xwav_start, parm, xwav_struct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fetchsort will sort detections into their respective xwav and identify
% the sample range start/end that should be processed to view the correct
% GPL detections. 
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sort detection into an xwav
StartSort = sort([mst_julian(n) xwav_start(1,:)]); % Sort
EndSort = sort([met_julian(n) xwav_start(1,:)]);
StartLoc = find(StartSort == mst_julian(n)); % Locate index
EndLoc = find(EndSort == met_julian(n));
if StartLoc == EndLoc
    if StartLoc == 1
        file_ident = StartLoc;
    else
        file_ident = StartLoc - 1; % File selected
    end
    error_flag = 0;
else % Detection crosses file
    error_flag = 4;
    return
end

% Create sample range possibilities for the file
file_string = sprintf('xwav%i',file_ident);
sample_range = 1:parm.nrec:xwav_struct.(file_string).TotalSamples;
if sample_range(end) < xwav_struct.(file_string).TotalSamples
    sample_range(end+1) = xwav_struct.(file_string).TotalSamples; % Assuming total samples is int. multiple of nrec
end

% Don't do deck test adjustment beyond disk 1
if file_ident > 1
    parm.pre_offset = 1;
end

% Calculate manual start and end time in samples since start of xwav  
detStart = round(seconds(datetime(mst_julian(n),'ConvertFrom','datenum')-datetime(xwav_start(file_ident),'ConvertFrom','datenum'))*parm.SampleFreq)+(parm.pre_offset-1)*parm.nrec; 
detEnd = round(seconds(datetime(met_julian(n),'ConvertFrom','datenum')-datetime(xwav_start(file_ident),'ConvertFrom','datenum'))*parm.SampleFreq)+(parm.pre_offset-1)*parm.nrec; 

% Sort start/end samples into sample range and determine relevant windows
if file_ident == StartLoc
    SampleRangeStart = find(sort([sample_range detStart]) == detStart);
    SampleRangeEnd = find(sort([sample_range detEnd]) == detEnd);
else
    SampleRangeStart = find(sort([sample_range detStart]) == detStart) - 1;
    SampleRangeEnd = find(sort([sample_range detEnd]) == detEnd) - 1;
end

% Find the detection start/end samples relative to start of window
detStart_rel = detStart - (SampleRangeStart-1)*parm.nrec;
detEnd_rel = detEnd - (SampleRangeEnd-1)*parm.nrec;