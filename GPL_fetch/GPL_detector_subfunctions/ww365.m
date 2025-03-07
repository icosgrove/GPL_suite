function [waveform,length_waveform] = ww365(x_ref,cm,parm,sp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will recreate the time-domain waveform of the strongest
% contour only. It will remove all low-noise regions around the contour,
% apply modulation, and produce a vector of the waveform with some
% zero-padding. 

% Written by Tyler Helble in 2014 version of GPL detector.
% Tested, documented, cleaned, established in GPL_v2 and v3 by Ian, but not
% fundamentally modified.
% 04/01/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialize local variables:
 
cm_copy = cm; % Strongest contour 
x = x_ref; % Raw audio data of cm
win = hamming(parm.fftl); % Hamming widndow
sp_squared = sp.^2; % element-wise squared spectrogram


%% Process Strongest Contour 

% Replace contour pixel locations with '1'
cm_max_binary = cm_copy;
cm_max_binary(cm_max_binary ~= 0) = 1; 

% Convolve and replace pixel values with '1' again. (2-bin
% paddding around the original contour)
low_noise_pixels = conv2(cm_max_binary,ones(3,3)/9,'same'); 
low_noise_pixels(cm_max_binary ~= 0) = 1; 

% Locate pixels NOT part of the contour
low_noise_pixels = 1 - low_noise_pixels; 
low_noise_index = find(low_noise_pixels); 

% Make number of low-energy pixels even if necessary
num_low_noise_indices = length(low_noise_index); % # of non-contour pixels
if mod(num_low_noise_indices,2) == 1 
 num_low_noise_indices = num_low_noise_indices - 1; % make length even
end

% Rare case if there are NO low-energy pixels
if num_low_noise_indices == 0 
    waveform = nan; 
    length_waveform = 0;
    return
end

% Use 2-D convolution to smooth the contour many times, each time replacing
% the areas of low-noise back to their original values. 
sp_copy = sp;
low_noise_copy = sp_copy(low_noise_index); 
for k = 1:500 
    sp_copy = conv2(sp_copy,ones(3,3)/9,'same'); % Convolve
    sp_copy(low_noise_index) = low_noise_copy; % Restore low-energy regions
end

 
% Locate the values of cm_max from the convolved spectrogram from the prev.
% step.
[cm_y,cm_x] = find(cm_copy); 
cm_max_ind = (parm.FreqBinHi - parm.FreqBinLo+1)*(cm_x-1) + cm_y; 
cm_max_smooth_vals = sp_copy(cm_max_ind); 

% Subtract smoothed values from absolute value of original contour sp
sp_copy2 = sp;
sp_copy2(cm_max_ind) = sp_squared(cm_max_ind).^.5 - cm_max_smooth_vals; 

% Normalize previous results, only save indices that make up cm_max
sp_copy2(sp_copy2 < 0) = 0; 
ratio = (sp_copy2./sp_squared.^.5).*cm_max_binary;

% Calculate number of time bins to loop over
offset = (parm.fftl-2*parm.fftOverlap)/2; 
sample_sec = ((length(x) - parm.fftl) / parm.fftOverlap) + 1;

% Contour is shorter than FFTL, too short so return. 
if sample_sec <= 0 
    waveform = nan; 
    length_waveform = 0;
    return
end

% Loop over time bins of the contour and compute each section of the
% time-domain waveform.
for j = 1:sample_sec 
    
    % Take FFT of the current time bin
    start = (j-1)*parm.fftOverlap + 1; 
    finish = start + parm.fftl - 1; 
    fft_col = fft(x(start:finish).*win); 

    % Create realized frequency range from FFTL, fill in the relevant bins
    % with enhanced sp produced in prior steps
    profile_filter = zeros(parm.fftl/2+1, 1); 
    profile_filter(parm.FreqBinLo:parm.FreqBinHi) = ratio(:,j); 
    profile_filter = profile_filter'; 
    
    % Append enhanced values to each other and scale using FFT values.
    profile_filter = [profile_filter,fliplr(profile_filter(2:end-1))]'; % Append
    fft_col = fft_col.*profile_filter; % Scale
    
    % Reproduce the time-domain signal using modulated/enhanced contour
    fft_col = real(ifft(fft_col))./win;
    timeD_signal(j,:) = fft_col(offset+1:end-offset); 
    timeD_signal_copy(j,:) = timeD_signal(j,:); 

end % Time bin loop

% Modulate the first row of the time-domain signal. 
z = linspace(-pi/2,pi/2,2*parm.fftOverlap);
left = cos(z).^2; 
waveform = left.*timeD_signal(1,:); 

% Loop through the time bibs and apply the modulation to the time-domain
% signal.
for j = 2:sample_sec 
    
    % Append space and add modulated signal
    waveform = [waveform,zeros(1,floor(parm.fftOverlap))];
    waveform(end-2*parm.fftOverlap+1:end) = waveform(end-2*parm.fftOverlap+1:end) + left.*timeD_signal_copy(j,:);

end
 
% Normalize and scale the waveform using original audio sample values
v = waveform/norm(waveform); % Normalize
length_waveform = length(waveform); 
remove_offset = x(1+offset:length_waveform+offset); 
waveform = v'*(v*remove_offset); 

% Add zero padding to the ends of the waveform, over the offset
pad = length(x) - length_waveform - offset; 
waveform = [zeros(1,floor(offset)),waveform',zeros(1,floor(pad))]'; 
