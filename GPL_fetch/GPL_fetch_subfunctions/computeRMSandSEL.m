function [rms_rl,sel] = computeRMSandSEL(data,parm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the RMS RL and finds the sound exposure level of
% the given snippet of a waveform. Units are dB re 1 count
% Reference: Madsen 2005
% Written: Ian 09/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = length(data);
rms_rl = 10*log10((1/T) * sum(data.^2));
sel = rms_rl + 10*log10(T);
