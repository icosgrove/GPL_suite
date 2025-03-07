function [sp,parm] = GPL_fft_fetch(data,parm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPL_fft take the Fast Fourier Transform of the current window of audio
% data with a given FFTL and overlap set by the user in the input
% parameters.

% Written by Tyler Helble
% Documented, tested, and slightly modified by Ian
% 04/08/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Pre-processing

win = hamming(parm.fftl); % Hamming Window
sp = zeros(parm.fftl , parm.NumTimeBins); % Spectrogram Window
x = data;

% Transpose to column vector if necessary
[x1,x2] = size(x);    
if x2 > x1
    x = x'; 
end 

%% Take FFT slices and create spectrogram

for j = 1:(parm.NumTimeBins) % parm.fftl/parm.fftOverlap)
    
    % Start/End Samples 
    start = (j-1)*parm.fftOverlap + 1; % Account for FFT overlap
    finish = start + parm.fftl - 1; 
    
    % FFT and Spectrogram 
    q = fft(x(start:finish).*win); % FFT
    sp(:,j) = abs(q); % Create Spectrogram

end

% Remove pre-padding
% sp = sp(:,parm.fftl/parm.fftOverlap:end);

