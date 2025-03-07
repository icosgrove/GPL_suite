function [calls_append, adhoc_append] = estFetchRL(calls_append,adhoc_append,parm,rl_data,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will estimate the peak-to-peak and rms RL of each call in 
% the GPL pairing and extra adhoc detections. 
% Written: Ian 08/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if parm.fp1.p2p_rms_RL == 1

    % Paired calls
    if ~isempty(calls_append.gpl_match)

            % Loop through pairings and process waveforms for RLs
            for k = 1:size(calls_append.gpl_match,2)

                % Pre-processing on data chunk
                stime = calls_append.gpl_match(k).start_time;
                endtime = calls_append.gpl_match(k).end_time;
                if (flag.pre == 0) && (flag.post == 0)
                    stime = stime + parm.nrec; 
                    endtime = endtime + parm.nrec;
                end
                if parm.fp1.rl_time_fix == 1 % Apply time delay adjustment
                    stime = stime + parm.fp1.tdelay;
                    endtime = endtime + parm.fp1.tdelay;
                end
                % % Apply large window that centers the call
                dur = endtime - stime + 1; 
                stime = stime - parm.fp1.rms_super_window_length*dur; 
                endtime = endtime + parm.fp1.rms_super_window_length*dur;
                if stime < 1 % ensure no 0/negative samples
                    stime = 1;
                end
                if endtime > length(rl_data) % ..and no over-indexing
                    endtime = length(rl_data);
                end
                data = rl_data(stime:endtime); % waveform + large padding

                % Continue processing on each contour type, if they exist
                if ~isempty(calls_append.gpl_match(k).cm.values) % cm: all contours
                    [p2p_RL,rms_RL,sel,peak_freq,RL_note] = computeFetchRL(calls_append.gpl_match(k).cm,parm,data,dur);
                    calls_append.gpl_match(k).cm.peak_freq = peak_freq;
                    calls_append.gpl_match(k).cm.p2p_RL = p2p_RL;
                    calls_append.gpl_match(k).cm.rms_RL = rms_RL;
                    calls_append.gpl_match(k).cm.SEL = sel;
                    calls_append.gpl_match(k).cm.RL_note = RL_note;
                end
                if ~isempty(calls_append.gpl_match(k).cm_max.values) % cm_max: single strongest contour
                    [p2p_RL,rms_RL,sel,peak_freq,RL_note] = computeFetchRL(calls_append.gpl_match(k).cm_max,parm,data,dur);
                    calls_append.gpl_match(k).cm_max.peak_freq = peak_freq;
                    calls_append.gpl_match(k).cm_max.p2p_RL = p2p_RL;
                    calls_append.gpl_match(k).cm_max.rms_RL = rms_RL;
                    calls_append.gpl_match(k).cm_max.SEL = sel;
                    calls_append.gpl_match(k).cm_max.RL_note = RL_note;
                end
                if ~isempty(calls_append.gpl_match(k).cm_max2.values) % cm_max2: 2nd strongest contour
                    [p2p_RL,rms_RL,sel,peak_freq,RL_note] = computeFetchRL(calls_append.gpl_match(k).cm_max2,parm,data,dur);
                    calls_append.gpl_match(k).cm_max2.peak_freq = peak_freq;
                    calls_append.gpl_match(k).cm_max2.p2p_RL = p2p_RL;
                    calls_append.gpl_match(k).cm_max2.rms_RL = rms_RL;
                    calls_append.gpl_match(k).cm_max2.SEL = sel;
                    calls_append.gpl_match(k).cm_max2.RL_note = RL_note;
                end

            end % For: loop over paired call(s)
    end % GPL matches

    % Adhoc detections
    if ~isempty(adhoc_append)
        for k = 1:size(adhoc_append,1)

            % Pre-processing on data chunk
            stime = adhoc_append(k).start_time;
            endtime = adhoc_append(k).end_time;
            if (flag.pre == 0) && (flag.post == 0)
                stime = stime + parm.nrec; 
                endtime = endtime + parm.nrec;
            end
            if parm.fp1.rl_time_fix == 1 % Apply time delay adjustment
                stime = stime + parm.fp1.tdelay;
                endtime = endtime + parm.fp1.tdelay;
            end
            % % Apply large window that centers the call
            dur = endtime - stime + 1; 
            stime = stime - parm.fp1.rms_super_window_length*dur; 
            endtime = endtime + parm.fp1.rms_super_window_length*dur;
            if stime < 1 % ensure no 0/negative samples
                stime = 1;
            end
            if endtime > length(rl_data) % ..and no over-indexing
                endtime = length(rl_data);
            end
            data = rl_data(stime:endtime); % waveform + large padding

            % Continue processing on each contour type, if they exist
            if ~isempty(adhoc_append(k).cm.values) % cm: all contours
                [p2p_RL,rms_RL,sel,peak_freq,RL_note] = computeFetchRL(adhoc_append(k).cm,parm,data,dur);
                adhoc_append(k).cm.peak_freq = peak_freq;
                adhoc_append(k).cm.p2p_RL = p2p_RL;
                adhoc_append(k).cm.rms_RL = rms_RL;
                adhoc_append(k).cm.SEL = sel;
                adhoc_append(k).cm.RL_note = RL_note;
            end
            if ~isempty(adhoc_append(k).cm_max.values) % cm: all contours
                [p2p_RL,rms_RL,sel,peak_freq,RL_note] = computeFetchRL(adhoc_append(k).cm_max,parm,data,dur);
                adhoc_append(k).cm_max.peak_freq = peak_freq;
                adhoc_append(k).cm_max.p2p_RL = p2p_RL;
                adhoc_append(k).cm_max.rms_RL = rms_RL;
                adhoc_append(k).cm_max.SEL = sel;
                adhoc_append(k).cm_max.RL_note = RL_note;
            end
            if ~isempty(adhoc_append(k).cm_max2.values) % cm: all contours
                [p2p_RL,rms_RL,sel,peak_freq,RL_note] = computeFetchRL(adhoc_append(k).cm_max2,parm,data,dur);
                adhoc_append(k).cm_max2.peak_freq = peak_freq;
                adhoc_append(k).cm_max2.p2p_RL = p2p_RL;
                adhoc_append(k).cm_max2.rms_RL = rms_RL;
                adhoc_append(k).cm_max2.SEL = sel;
                adhoc_append(k).cm_max2.RL_note = RL_note;
            end

        end
    end
end