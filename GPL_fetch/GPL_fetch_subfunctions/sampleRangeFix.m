function [start_fix, end_fix] = sampleRangeFix(sample_range, SampleRangeStart, parm, xwav_struct, file_string)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will adjust the sample range in order to account for
% beginning and end cases of the audio data
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% No time delay for the first window 
parm = timeDelayConfig(parm,SampleRangeStart);

% Adjust sample range for overlap delay issue
if sample_range(SampleRangeStart)-(parm.time_delay*parm.SampleFreq) < 0
    start_fix = 1;
else
    if parm.fetch.time_delay_switch == 0
        start_fix = sample_range(SampleRangeStart);
    else
        start_fix = sample_range(SampleRangeStart)-(parm.time_delay*parm.SampleFreq);
    end
end
if (sample_range(SampleRangeStart)+parm.nrec+parm.fftOverlapPercent*0.01*parm.fftl-1-(parm.time_delay*parm.SampleFreq)) > xwav_struct.(file_string).TotalSamples
    end_fix = xwav_struct.(file_string).TotalSamples;
else
    if parm.fetch.time_delay_switch == 0
        end_fix = sample_range(SampleRangeStart)+parm.nrec-1;
    else
        end_fix = sample_range(SampleRangeStart)+parm.nrec+parm.fftOverlapPercent*0.01*parm.fftl-1-(parm.time_delay*parm.SampleFreq);
    end
end
