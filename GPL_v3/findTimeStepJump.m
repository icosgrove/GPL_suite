function [j_flag, new_startTime,ideal_idx] = findTimeStepJump(param_s,parm,num_loops,step_size)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will locate a jump in time due to a deck test. This occurs
% when the hydrophone is turned on on a ship for some time, and then the
% effort jumps to when it is in the water.
% Written: Ian 08/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Expected time jump between windows
expTimeStep = datenum([0 0 0 0 0 parm.window_length]);

% Find jumps in time larger than expected+3% margin
df = diff(param_s.ltsahd.dnumStart);
jump = find(df > (expTimeStep+0.03*expTimeStep)); 
if isempty(jump) % no deck test
    j_flag = [];
    new_startTime = [];
    ideal_idx = [];
    disp('No deck test detected, please turn off parm.deck_test_adjustment.')
    return
end

% Find the first sample start time after the deck test.
if parm.pad == 1
    ideal_sample_range = 1:parm.nrec:(length(df)*parm.nrec);
    sample_range = 1:step_size:(num_loops*step_size+step_size);
    ideal_samp = parm.nrec*jump;
    samp_loc = sort([sample_range ideal_samp])==ideal_samp;
    j_flag = sample_range(samp_loc);
    new_startTime = find(j_flag==sample_range);
    ideal_idx = find(sort([j_flag ideal_sample_range])==j_flag);
else
    sample_range = 1:step_size:(num_loops*step_size+step_size);
    j_flag = sample_range(jump+1);
    new_startTime = [];
    ideal_idx = [];
end

if (length(j_flag) > 1) || (length(new_startTime) > 1) || (length(ideal_idx) > 1)
    disp('Issue with deck test timestamp adjustment.')
    pause
end