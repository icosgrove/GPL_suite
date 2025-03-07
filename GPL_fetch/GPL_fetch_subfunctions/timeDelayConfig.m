function parm = timeDelayConfig(parm,SampleRangeStart)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will set the time delay as needed to handle the correct
% sample range for overlapping FFT windows. The delay is ignored if it is
% the first set of samples in the data or turned off entirely.
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (SampleRangeStart == parm.pre_offset) || (SampleRangeStart == 1)
    parm.time_delay = 0;
else
    if parm.fetch.time_delay_switch == 1
        parm.time_delay = parm.fftOverlapPercent*1e-2*parm.fftl/parm.SampleFreq/2;
    else
        parm.time_delay = 0;
    end
end