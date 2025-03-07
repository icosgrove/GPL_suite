function [freq_f,time,sp_context,parm] = makeTritonSp(sub_data,parm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makeTritonSp.m will create a spectrogram of the given subset of data
% using the same methods that Triton uses. The option to apply a transfer
% function and control brightness and contrast is available. 
% Code is sourced from plot_specgram.m and mkspecgram.m from Triton version
% 1.95.20230315. All credit goes to the Scripps Whale Acoustics Laboratory.
% This is my recreation of their code. 
% Written: Ian 08/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spectrogram creation
window = hanning(parm.fftl);
if length(window) > length(sub_data)
    window = hanning(length(sub_data));
end
nOverlapSamp = (parm.fftOverlapPercent/100)*parm.fftl;
if nOverlapSamp > length(sub_data)
    nOverlapSamp = length(sub_data) -1;
end
[~,freq_f,time,sp_context] = spectrogram(sub_data,window,nOverlapSamp,parm.fftl,parm.SampleFreq);
sp_context = sp_context(parm.fp1.ContextFreqBinLo:parm.fp1.ContextFreqBinHi,:);
c_freqRange = freq_f(parm.fp1.ContextFreqBinLo:parm.fp1.ContextFreqBinHi);

% Apply first TF to obtain units of dB re 1 counts^2/Hz
sp_context = 10*log10(abs(sp_context));

% Apply main TF to obtain units of dB re 1 uPa^2/Hz
if parm.fp1.tf.switch == 1
    if exist(parm.fp1.tf.file)
        tfunc = load(parm.fp1.tf.file);
        parm.fp1.tf.freq = tfunc(:,1);
        parm.fp1.tf.vals = tfunc(:,2);
        [~,ia_f,ic_f] = unique(parm.fp1.tf.freq); % check for repeated values and discard
        if length(ia_f) == length(ic_f)
            freq_f = parm.fp1.tf.freq;
            uppc_f = parm.fp1.tf.vals;
        else
            freq_f = parm.fp1.tf.freq(ia_f);
            uppc_f = parm.fp1.tf.vals(ia_f);
        end
        extrap_tf = interp1(freq_f,uppc_f,c_freqRange,'linear','extrap');
        sz_f = size(sp_context);
        bwdb_f = 10*log10(parm.fftl/parm.SampleFreq);
        sp_context = extrap_tf*ones(1,sz_f(2)) + sp_context + bwdb_f*ones(sz_f(1),sz_f(2));
    else
        dispErrorFlag(24)
    end
end

% Apply Brightness/Contrast
sp_context = (parm.fp1.contrast/100).*sp_context + parm.fp1.brightness;

% Configure Colorbar
if parm.fp1.colorbar.switch == 1
    parm.fp1.colorbar.Cmin = (parm.fp1.contrast / 100) * parm.fp1.noise_floor + parm.fp1.brightness;
    parm.fp1.colorbar.YLabelString = '';
    if parm.fp1.tf.switch == 1
        parm.fp1.colorbar.YLabelString = 'Spectrum Level [dB re 1\muPa^2/Hz]';
    else
        parm.fp1.colorbar.YLabelString = 'Spectrum Level [dB re counts^2/Hz]';
    end
    parm.fp1.colorbar.sp_context_min = min(min(sp_context));
    parm.fp1.colorbar.sp_context_max = max(max(sp_context));    
    if isinf(parm.fp1.colorbar.sp_context_min)
        parm.fp1.colorbar.sp_context_min = parm.fp1.colorbar.Cmin;
    end
    if isinf(parm.fp1.colorbar.sp_context_max)
        parm.fp1.colorbar.sp_context_max = 100;
    end
    if parm.fp1.colorbar.sp_context_min < parm.fp1.colorbar.Cmin
        parm.fp1.colorbar.sp_context_min = parm.fp1.colorbar.Cmin;
    end
end
