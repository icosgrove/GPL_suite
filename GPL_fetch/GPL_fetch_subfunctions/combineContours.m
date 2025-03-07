function combCM = combineContours(cm1,cm2,parm,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will combine adjacent contours together as one, for the
% cases where a call becomes split due to it crossing a window boundary
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create both slates
mat1 = zeros(cm1.size(1), cm1.size(2));
mat2 = zeros(cm2.size(1), cm2.size(2));

% Fill in slates with contours
mat1(cm1.index) = double(cm1.values)*cm1.scale./(2^16-1);
mat2(cm2.index) = double(cm2.values)*cm2.scale./(2^16-1);

if size(mat2,1) == (size(mat1,1) + 2) % cm_max2 extra bins issue
    mat2 = mat2(2:end-1,:);
elseif size(mat1,1) == (size(mat2,1) + 2)
    mat1 = mat1(2:end-1,:);
end

% Combine slates 
newSlate = [mat1 mat2];

% Create new contour measurements
index = find(newSlate);
values = newSlate(index);
scale = max(values);
combCM.values = newSlate(index)/scale*(2^16-1);
combCM.index = index;
combCM.scale = scale;
combCM.size = size(newSlate);

if isfield(cm1,'duration_sec') && isfield(cm2,'duration_sec') % Duration
    nonz = find(sum(newSlate)); % Locate nonzero columns
    num_Tbins = length(nonz); 
    combCM.duration_sec = num_Tbins * (parm.fftOverlap/parm.SampleFreq); % seconds
    combCM.duration_bin = num_Tbins; % bins
end
if isfield(cm1,'slope') && isfield(cm2,'slope') % Slope
    [slope_ci, slope, slope_y_int] = GPL_slope(newSlate,parm);
    combCM.slope_ci = slope_ci;
    combCM.slope = slope*(parm.SampleFreq/parm.fftl)/(parm.fftOverlap/parm.SampleFreq);
    combCM.slope_y_int = slope_y_int;
end
if isfield(cm1,'freq_bandwidth_hz') && isfield(cm2,'freq_bandwidth_hz') % Freq. Bandwidth
    num_Fbins = length(find(sum(newSlate,2))); % Find # of nonzero freq bins
    combCM.freq_bandwidth_hz = num_Fbins*(parm.SampleFreq/parm.fftl); % Hz
    combCM.freq_bandwidth_bin = num_Fbins; % bins
end
if isfield(cm1,'start_freq_hz') && isfield(cm2,'start_freq_hz')
    [start_freq, end_freq, min_freq, max_freq, abs_bandwidth_hz] = GPL_freq_measurements(newSlate,parm);
    combCM.start_freq_hz = start_freq;
    combCM.end_freq_hz = end_freq;
    combCM.min_freq_hz = min_freq;
    combCM.max_freq_hz = max_freq;
    combCM.abs_bandwidth_hz = abs_bandwidth_hz;
end

% figure(1)
% subplot(1,3,1); imagesc(mat1)
% subplot(1,3,2); imagesc(mat2)
% subplot(1,3,3); imagesc(newSlate)