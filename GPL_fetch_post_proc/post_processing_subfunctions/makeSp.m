function [freq,time,sp,figPARAMS] = makeSp(data,figPARAMS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create spectrograms using the spectrogram.m function and following the
% methods used in Triton
% Some code repurposed from Triton: credit goes to the Scripps Whale
% Acoustics Laboratory.
% Written: Ian 09/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if figPARAMS.fftl > length(data)
    window = hanning(length(data));
else
    window = hanning(figPARAMS.fftl);
end
nOverlapSamp = (figPARAMS.overlapPerc/100)*figPARAMS.fftl;
if nOverlapSamp > length(data)
    disp('Warning: Overlap % is too large for call size')
    nOverlapSamp = length(data)-1;
end
[~,freq,time,sp] = spectrogram(data,window,nOverlapSamp,figPARAMS.fftl,figPARAMS.sampleFreq);

% Find indices of frequency start/end
figPARAMS.freq_loc_min = find(freq==max(freq(freq <= figPARAMS.freqLo)));
figPARAMS.freq_loc_max = find(freq==min(freq(freq >= figPARAMS.freqHi)));
sp = sp(figPARAMS.freq_loc_min:figPARAMS.freq_loc_max,:);
c_freqRange = freq(figPARAMS.freq_loc_min:figPARAMS.freq_loc_max);

% Apply main TF to obtain units of dB re 1 uPa^2/Hz
if figPARAMS.tf_on == 1
    tfreq = figPARAMS.tf(:,1);
    vals = figPARAMS.tf(:,2);
    [~,ia_f,ic_f] = unique(tfreq); % check for repeated values and discard
    if length(ia_f) == length(ic_f)
        freq_f = tfreq;
        uppc_f = vals;
    else
        freq_f = tfreq(ia_f);
        uppc_f = vals(ia_f);
    end
    extrap_tf = interp1(freq_f,uppc_f,c_freqRange,'linear','extrap');
    sz_f = size(sp);
    bwdb_f = 10*log10(figPARAMS.fftl/figPARAMS.sampleFreq);
    sp = extrap_tf*ones(1,sz_f(2)) + sp + bwdb_f*ones(sz_f(1),sz_f(2));
end

% Apply Brightness/Contrast
sp = (figPARAMS.contrast/100).*sp + figPARAMS.brightness;

% Configure Colorbar
if figPARAMS.colorbar.on == 1
    figPARAMS.colorbar.Cmin = (figPARAMS.contrast / 100) * -40 + figPARAMS.brightness;
    figPARAMS.colorbar.YLabelString = '';
    if figPARAMS.tf_on == 1
        figPARAMS.colorbar.YLabelString = 'Spectrum Level [dB re 1\muPa^2/Hz]';
    else
        figPARAMS.colorbar.YLabelString = 'Spectrum Level [dB re counts^2/Hz]';
    end
    figPARAMS.colorbar.sp_min = min(min(sp));
    figPARAMS.colorbar.sp_max = max(max(sp));    
    if isinf(figPARAMS.colorbar.sp_min)
        figPARAMS.colorbar.sp_min = figPARAMS.colorbar.Cmin;
    end
    if isinf(figPARAMS.colorbar.sp_max)
        figPARAMS.colorbar.sp_max = 100;
    end
    if figPARAMS.colorbar.sp_min < figPARAMS.colorbar.Cmin
        figPARAMS.colorbar.sp_min = figPARAMS.colorbar.Cmin;
    end
end