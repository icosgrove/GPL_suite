function [GPL_struct] = GPL_v3(sub_data,parm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function GPL_v3 will process a window of xwav audio data through the
% GPL algorithm. That algorithm will identify signsls of interest as
% indicated by the input parameters. The selected calls will have their
% contours extracted, those contours will be measured, and preliminary
% filters will be applied to the detections. 

% Theory behind the GPL algorithm was created and code was written by Tyler
% Helble. Testing, documentation, cleaning, some modifications, and some
% new (as of v3) were written by Ian Cosgrove with Joshua Jones.
% 02/03/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Take the FFT and produce a spectrogram of the data - Tyler Helble
[sp] = GPL_fft(sub_data,parm);

% In anticipation of FIFO removal, GPL_fft does not restrict the frequency
% range of the produced spectrogram to the requested limits, this is done 
% after FIFO removal (If FIFO removal is not wanted, it will do it next).


%% FIFO Removal - Ian Cosgrove 
% Optional removal of 50 Hz harmonics from the spectrogram before further
% processing. 
if parm.fifo_removal == 1
    
    [sp] = GPL_fifo_removal(sp,parm);
    
end


% Restrict spectrogram to requested frequency range
sp = sp(parm.FreqBinLo:parm.FreqBinHi,:);


%% Whiten and spectrogram - Tyler Helble
[sp_whiten] = GPL_whiten(sp,parm);

% Helble et al. 2012 eq 7

% sp_whiten: The whitened spectrogram. Mean noise values for each frequency
% (and optionally time) bin are found and used to reduce and normalize the
% values of the spectrogram. A 2-D convolution is also performed.


% Detections will be found with whitened spectrogram later
sp_loop = sp_whiten;



%% Extract and replictae quiet (non-signal) portions of the slate - Tyler Helble
[quiet_whiten, quiet_fft, quiet_base, noise_floor, blocked, baseline0] = GPL_quiet(sp,sp_whiten,parm);

% quiet_whiten: Signal absent portions of the whitened spectrogram 
% quiet_fft: Signal absent portions of original spectrogram
% quiet_base: Mean noise level of entire slate
% noise_floor: Minimum noise value for the current window
% blocked: A vector where '0' is signal presence and '1' is signal absense,
% one value for each time bin.
% baseline0: Energy sums for each time bin, used for identifying signal
% presence/absence.






%% Apply GPL Algorithm to find Detections - Tyler Helble

% Loop variable allocation
mcalls = []; 
iteration_count = 0; 
num_calls = 2; 
 
% Identify bin ranges for summation later
low_off = parm.SumFreqBinLo - parm.FreqBinLo + 1; 
high_off = parm.SumFreqBinHi - parm.FreqBinLo + 1;

% Loop to find units above detection threshold: Find calls, remove them,
% search the slate for more calls, and repeat for the number of loops
% defined by NumLoops parameter OR terminate if no calls are found in the 
% first loop, or any successive loop.
while iteration_count < parm.NumLoops && num_calls > 1 
       
    % Normalize looping (whitened) spectrogram by RSS values in frequency and time   
    norm_v = sp_loop./(ones(parm.NumFreqBins,1)*sum(sp_loop.^2).^(1/2)); % Time
    norm_h = sp_loop./(sum(sp_loop'.^2).^(1/2)'*ones(1,parm.NumTimeBins)); % Frequency
    
    % Whiten both matrices
    norm_v = whiten_matrix(norm_v')';
    norm_h = whiten_matrix(norm_h);
 
    % Use exponents to find signal presence
    bas = abs(norm_v).^parm.xp1.*abs(norm_h).^parm.xp2;
    
    % Find the energy sum for each solumn, whiten, and normalize by mea
    % noise value of the entire slate
    base_in = sum(bas(low_off:high_off,:).^2); 
    [b0] = whiten_vec(base_in'); 
    base_in = b0'/quiet_base; 
  
    
    
    %%% Locate signal present and signal absent time bins
    [base_out,calls] = GPL_cropping(base_in,parm.noise_ceiling*noise_floor,...
                             parm.thresh*noise_floor);   
    % Outputs of GPL_cropping:
    % base_out: Vector with detections set to '0'
    % calls: Time bin start/end indices of each detection
         
    
             
        % If no calls are found, lc = 0, and loop will be terminated.                  
        [num_calls,~] = size(calls); 
   
        
    
    %%% Finalize detgection start/end times found in GPL_cropping
    unblocked = zeros(parm.NumTimeBins,1); 
    
    % Set the indices of detections to be a '1' for vector of the entire window.
    for j = 1:num_calls 
        unblocked(calls(j,1):calls(j,2)) = 1; 
    end

    % Remove signal absent bins AND the straddle bins, since they are below
    % noise ceiling (they bordered a detection on the quiet side)
    tst = unblocked - blocked;
    tst(tst < 1) = 0; 

    
    % Locate index where 0 (signal absense) turns to a 1 (signal presence).
    % Those are the start bins of each detection
    st = find(diff(tst) > 0) + 1; 
    
        % Diff won't work if the first bin is part of a detection, so that case
        % is handled by appending a 1 to leading element.
        if tst(1) == 1 
            st = [1,st']'; 
        end
    
    % Locate index where a 1 (signal presence) changes to a 0 (signal
    % absence). There are the finish bins of each detection.
    fn = find(diff(tst) < 0); 
    
        % Account for the case where diff won't work 
        if tst(end) == 1
            fn = [fn',parm.NumTimeBins]';
        end
    
    % Start and end times of signal presence (detection) are saved. 
    calls = [st,fn]; 

    %%% Pre-merging detection is complete 

    
    
    
    
    %%% Merge adjacent calls together than fall within a bin range. This
    %%% accounts for calls where there is an energy drop within the call.
    
    if parm.filter_parm.switches.AdjacentCallMerger == 0
    
        % Locate detections where the start and end times are within the
        % maximum bin cutoff for two adjacent detections to be merged into one.
        [~,k2] = sort(calls(:,1)); 
        st = calls(k2,1);
        fn = calls(k2,2);  
        k = find(st(2:end) - fn(1:end-1) < parm.filter_parm.values.AdjacentCallBinNum);

        % For adjacent calls within cutoff, remove the end/start time that
        % falsely split up the overall call.
        omit = length(k); 
        if omit > 0
           n = 1:length(st);
           st = st(setdiff(n,k + 1)); 
           fn = fn(setdiff(n,k));
        end
        calls = [st,fn]; % Redefine calls with falsely split detections combined
    
    end % End call merger
    
    
        % Locate zero length calls and remove them
        dur = diff(calls'); 
        calls = calls(dur ~= 0,:); 
    
    
    % Append calls from future loops onto list already created
    mcalls = [mcalls',calls']'; 
    
    % Step loop iteration counter forward
    iteration_count = iteration_count + 1; 
    
    % Find number of new calls found, if no new ones are found the loop
    % terminates.
    [num_calls,~] = size(calls); 
   
     
    % Block out time bins of new calls found in current loop so they are
    % not repeated in future loops.
    blocked = blocked + unblocked; 
    blocked(blocked > 1) = 1; % Some old calls may add to 2, adjust that to 1 

    
    
    % Iteration change: Locate time bins that are correlated with a
    % detection, remove then from the quiet whitened spectrogram and
    % continue the loop with these quieter bins. 
    ksel = find(base_out == 0);
    
    % Replace the signal present windows with signal absent whitened
    % values. This ensures they will not be processed again
    sp_loop(:,ksel) = quiet_whiten(:,ksel); 

end % End while: Unit selection is complete


%%%%%%%%%%%%%%%%%% Detection identification is complete %%%%%%%%%%%%%%%%%%




%% Set a flag if a detection hits the window boundary (only for non-padding process)
% Check if the final detections end time is the same as the length of the
% window.

boundary_cross_flag_start = 0; % By default, stays zero
boundary_cross_flag_end = 0;

if parm.pad == 0
    if numel(mcalls) > 0 
        end_time = mcalls(end,2);
        start_time = mcalls(1,1);
        if end_time == length(base_in) % Detection reaches window end
            boundary_cross_flag_end = 1;
        elseif start_time == 1 % Detection starts immediately (assumed rollover)
            boundary_cross_flag_start = 1;
        end
    end
end



%% Call Duration Filter - Tyler Helble

% First round of filters is based only on call duration
[start,finish,times] = GPL_CallDurationFilter(mcalls,parm);
% start: call start times in bins relative to window start
% finish: call finish times in bins relative to window start
% times: start/finish times in samples relative to window start


%% If Tyler's padding is being used, combine calls that cross slate boundaries
% Detection start times within parm.pad*1 seconds of window end will be
% ignored and processed in the next window as it can be assumed that they
% cross the boundary. 

if parm.pad ~= 0 
    siz = size(times);
    if siz(1) > 0

        t_pad = parm.pad*parm.SampleFreq; % Calculate padding amount 
        k1 = find(times(:,1) > t_pad); % Save detection start times after left padding    
        k2 = find(times(:,1) <= parm.nrec - t_pad); % " " before right padding 
        k1 = intersect(k1,k2); 
        
        % Save detections only within the padded times
        times = times(k1,:);
        start = start(k1,:);
        finish = finish(k1,:);
        
    end

end

%% Perform call contour creation and measurements, load data in GPL_struct 

% Create GPL_struct/clear it from previous window
GPL_struct = [];

ncalls = size(times);

if ncalls(1) > 0 % Skip this step if there are no detections
    
    for k = 1:ncalls(1) % Loop through all calls

        % Load call start/end times into GPL_struct
        GPL_struct(k).start_time = times(k,1); 
        GPL_struct(k).end_time = times(k,2); 
        
        % Load Boundary Flag into GPL_struct
        if (boundary_cross_flag_end == 1) || (boundary_cross_flag_start == 1) % One/both flags on
            if k == ncalls(1) % Potential end case 
                if boundary_cross_flag_end == 1 % Confirmed end case
                    GPL_struct(k).boundary_cross_flag = 1; 
                else
                    GPL_struct(k).boundary_cross_flag = 0;
                end
            elseif k == 1 % Potential start case
                if boundary_cross_flag_start == 1 % Confirmed start case
                    GPL_struct(k).boundary_cross_flag = 1;
                else
                    GPL_struct(k).boundary_cross_flag = 0;
                end
            else % Middle case (no possibility of crossing)
                GPL_struct(k).boundary_cross_flag = 0;
            end
        else % No flag
            GPL_struct(k).boundary_cross_flag = 0;
        end

    end % For: Loop over calls

    %% Call Contour Creation - Tyler Helble
    % Process each detection one by one, extract three contours of the
    % detection, one showing all time-frequency locations of high
    % energy that triggered the detection in the first place, one 
    % showing the single strongest 'island' that makes up the
    % detection, and the single second strongest island if applicable.
    % A time-domain waveform may also be extracted for the strongest
    % contour. 
    if parm.template_on == 1 

        [GPL_struct] = GPL_template(GPL_struct,sp,sp_whiten,start,finish,sub_data,parm);

    end


    %% Call Contour Measurements - Tyler Helble
    % Measure the duration, slope, and frequency characteristics of the
    % contours.
    if parm.measurements_on == 1

        [GPL_struct] = GPL_measurements(GPL_struct,sp,quiet_fft,start,finish,parm);

    end
    
    
    %% Filter based on contour measurements
    % Filter structure is creating indicating 0/1 if a contour passes or
    % not. Removing calls based on these can be done in an external func.
    if parm.filter_on == 1
        
         [GPL_struct] = GPL_filter(GPL_struct,parm);

    end

    
end % If: At least one detection 


%% Optional Plotting: Will display for each window and pause the algorithm - Tyler Helble
% Display a plot of the whitened spectrogram, energy summation and
% detection thresholds, and signal-absent whitened spectrogram.

if parm.plot == 1
    
    time_range = (1:length(baseline0))*parm.fftOverlap/parm.SampleFreq;
  
    figure(7); 
    
    % Plot the whitened spectrogram 
    subplot(3,1,1)
    freq_range = linspace(parm.FreqBinLo,parm.FreqBinHi,parm.FreqBinHi - parm.FreqBinLo + 1);
    imagesc(time_range, freq_range, 20*log10(abs(sp_whiten)/max(max(sp_whiten))),[-40,0]);
    axis xy; title('Whitened Spectrogram');
    xlabel('Time [s]'); ylabel('Frequency [Hz]'); 

    % Plot energy sum, detection thresholds, and detection regions
    subplot(3,1,2); hold on
    plot(time_range, log10(abs(baseline0)),'k'); % Energy Sum as black line
    for j = 1:length(start) % Plot detection regions as red line
        plot(time_range(start(j):finish(j)),log10(abs(baseline0(start(j):finish(j)))),'r');
    end
    noise_ceil_line = log10(parm.noise_ceiling*noise_floor); % Find noise_floor
    det_thresh_line = log10(parm.thresh*noise_floor); % find the noise ceiling
    plot(time_range,noise_ceil_line+0*time_range,'k--'); % plot noise floor as black dashed
    plot(time_range,det_thresh_line+0*time_range,'g--'); % plot noise ceiling as green dashed
    xlabel('Time [s]'); ylabel('Energy')
    axis tight
    hold off
    
    % Plot whitened signal ABSENT whitened spectrogram
    subplot(3,1,3)
    imagesc(time_range,freq_range,20*log10(abs(sp_loop)/max(max(sp_whiten))),[-40,0]); 
    axis xy; title('Signal Absent Whitened Spectrogram');
    xlabel('Time [s]'); ylabel('Frequency [Hz]');
    
    drawnow % overwrite previous 75s plot
    fprintf('\nPlot displayed: press any key to continue.\n')
    pause % wait for any key strike
    
    clf
    
end % If: plotting is on

