function [GPL_struct,fetched_call,error_flag,adhoc_detection,window_flag,endflag,userErrorFlag,parm] = GPL_v3_fetch(sub_data,samp_s,samp_e,parm,window_jst,mst_j,met_j,mst_num,cross_flag,process_num,manual_log,allowAdhocFlag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function GPL_v3 will process a window of xwav audio data through the
% GPL algorithm. That algorithm will identify signsls of interest as
% indicated by the input parameters. The selected calls will have their
% contours extracted, those contours will be measured, and preliminary
% filters will be applied to the detections. 

% Theory behind the GPL algorithm was created and code was written by Tyler
% Helble. Testing, documentation, cleaning, some modifications, and some
% new functions (as of v3) were written by Ian Cosgrove with Joshua Jones.
% 02/03/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

redo_flag = 0;
error_flag = 0;
adhoc_detection = [];
window_flag = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Take the FFT and produce a spectrogram of the data - Tyler Helble
[sp,parm] = GPL_fft_fetch(sub_data,parm);

% In anticipation of FIFO removal, GPL_fft does not restrict the frequency
% range of the produced spectrogram to the requested limits, this is done 
% after FIFO removal (If FIFO removal is not wanted, it will do it next).

% Load sample start into fetch struct
if process_num == 1
    if isnumeric(manual_log{mst_num,5}) % Switch out of datenum if necessary
        fetched_call.manual_start_time = datestr(manual_log{mst_num,5},'yyyy-mm-dd HH:MM:SS.FFF');
    else
        fetched_call.manual_start_time = string(manual_log{mst_num,5});
    end
end
fetched_call.manual_sample_start = samp_s;
fetched_call.manual_sample_end = samp_e;


%% FIFO Removal - Ian Cosgrove 
% Optional removal of 50 Hz (or your choice) harmonics from the spectrogram before further
% processing. 
if parm.fifo_removal == 1

    [sp] = GPL_fifo_removal(sp,parm);

end

% Save spectrogram of the full context frequency range 
if parm.fp1.triton_sp == 0 % Use GPL spectrogram for top subplot
    sp_context = sp(parm.fp1.ContextFreqBinLo:parm.fp1.ContextFreqBinHi,:);
    fetched_call.spectrogram = sp_context;
elseif parm.fp1.triton_sp == 1 % Use Triton methods for spectrogram creation, apply TF, B/C, colorbar
    [freq,time,sp_context,parm] = makeTritonSp(sub_data,parm);
    fetched_call.spectrogram = sp_context;
end


% Restrict GPL spectrogram to requested frequency range
sp = sp(parm.FreqBinLo:parm.FreqBinHi,:);


%% Whiten and spectrogram - Tyler Helble
[sp_whiten] = GPL_whiten(sp,parm);

% sp_whiten: The whitened spectrogram. Mean noise values for each frequency
% (and optionally time) bin are found and used to reduce and normalize the
% values of the spectrogram. A 2-D convolution is also performed.

if parm.fp1.whitened_sp == 1
    fetched_call.whitened_sp = sp_whiten;
end

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


%%%%%%%%%%%%%%%%%% Detection identification is complete %%%%%%%%%%%%%%%%%%%



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
    
end % If: At least one detection 



%% Place contours back into full length Spectrogram

contour_sp = zeros(length((parm.fp1.freqDelimitLo:parm.SampleFreq/parm.fftl:parm.fp1.freqDelimitHi)),parm.NumTimeBins);
for n = 1:length(start)
    [cm] = GPL_full('cm',n,GPL_struct);
    contour_sp((parm.FreqBinLo-parm.fp1.ContextFreqBinLo+1):(parm.fp1.ContextFreqBinHi-(parm.fp1.ContextFreqBinHi-parm.FreqBinHi)-parm.fp1.ContextFreqBinLo+1),start(n):finish(n)) = cm;
    GPL_struct(n).julian_start_time = window_jst+datenum([0 0 0 0 0 start(n)*parm.fftOverlap/parm.SampleFreq]);
    GPL_struct(n).julian_end_time = window_jst+datenum([0 0 0 0 0 finish(n)*parm.fftOverlap/parm.SampleFreq]);
end

% Cm: All Contour window
fetched_call.gpl_contour_window = contour_sp;


% Cm_max: Single strongest contour window (if requested)
if parm.fp1.cm_maxsubplot == 1
    cm_max_win = zeros(length((parm.fp1.freqDelimitLo:parm.SampleFreq/parm.fftl:parm.fp1.freqDelimitHi)),parm.NumTimeBins);
    for n = 1:length(start)
        [cm_max] = GPL_full('cm_max',n,GPL_struct);
        cm_max_win((parm.FreqBinLo-parm.fp1.ContextFreqBinLo+1):(parm.fp1.ContextFreqBinHi-(parm.fp1.ContextFreqBinHi-parm.FreqBinHi)-parm.fp1.ContextFreqBinLo+1),start(n):finish(n)) = cm_max;
    end
    fetched_call.cm_max_window = cm_max_win;
end

% Julian time params
if cross_flag == 0
    fetched_call.window_start_time = datestr(window_jst);
else
    fetched_call.window_start_time = datestr(window_jst+datenum([0 0 0 0 0 parm.nrec/parm.SampleFreq]));
end


%% Plot

% SI unit conversion
time_range = (1:length(baseline0))*parm.fftOverlap/parm.SampleFreq; % time in seconds
freq_range = parm.freq_range; % freq range
freq_range = freq_range(parm.fp1.ContextFreqBinLo:parm.fp1.ContextFreqBinHi);

% convert GPL/manual start times to seconds
start_sec = start*parm.fftOverlap/parm.SampleFreq;
finish_sec = finish*parm.fftOverlap/parm.SampleFreq;

%% FP1
if process_num == 1 % GPL-Manual Matching
    mstart_sec = seconds(datetime(mst_j,'convertFrom','datenum') - datetime(window_jst,'ConvertFrom','datenum'));
    mend_sec = seconds(datetime(met_j,'convertFrom','datenum') - datetime(window_jst,'ConvertFrom','datenum'));

    % Correct box locations for second half window overlap cases 
    if cross_flag > 0 
        mstart_sec = mstart_sec - cross_flag*parm.window_length;
        mend_sec = mend_sec - cross_flag*parm.window_length;
    end

    fig = resetFigure();
    t_s1 = sprintf('Manual Detection %i',mst_num); sgtitle(t_s1)

    % Full screen figure
    if parm.fp1.fullscreen == 1
        set(fig,'units','normalized','outerposition',[0 0 1 1])
    end
    
    % Plot the unchanged spectrogram OR Triton spectrogram
    subplot(3,1,1); hold on
    if parm.fp1.triton_sp == 1
        plot_tritonSp
    else
        image(time_range, freq_range, fetched_call.spectrogram)
        if cross_flag == 0
            [rectHandle] = det_boxes(start_sec, finish_sec, mstart_sec, mend_sec, parm, fetched_call, 0, 0, fig);
        else
            [rectHandle] = det_boxes(start_sec, finish_sec, mstart_sec, mend_sec, parm, fetched_call, cross_flag, 0, fig);
        end
        title('Spectrogram [GPL]') 
        xlabel('Time [s]'); ylabel('Frequency [Hz]')
        set(gca,'YDir','normal')
        text('Units', 'normalized', 'Position', [3/parm.window_length, -0.18, 0], ... % Window start time text
            'String', sprintf(datestr(window_jst+datenum([0 0 0 0 0 cross_flag*parm.window_length]))), 'HorizontalAlignment', 'right');
        axis([time_range(1) time_range(end) parm.fp1.freqDelimitLo parm.fp1.freqDelimitHi])
     end
    
    hold off

    % Plot energy sum, detection thresholds, and detection regions OR
    % cm_max contours
    subplot(3,1,2); hold on

    if parm.fp1.cm_maxsubplot == 1 % Cm_max Contours
        imagesc(time_range,freq_range,fetched_call.cm_max_window)
        if cross_flag == 0
            [rectHandle] = det_boxes(start_sec, finish_sec, mstart_sec, mend_sec, parm, fetched_call, 0, 0, fig);
        else
            [rectHandle] = det_boxes(start_sec, finish_sec, mstart_sec, mend_sec, parm, fetched_call, cross_flag, 0, fig);
        end
        xlabel('Time [s]'); ylabel('Frequency [Hz]'); 
        t_s = sprintf('Cm max: Single Strongest GPL Contour\n'); title(t_s)
        plot3([0 parm.window_length],[parm.freq_lo parm.freq_lo],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'--w','LineWidth',1)
        plot3([0 parm.window_length],[parm.freq_hi parm.freq_hi],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'--w','LineWidth',1)   
        axis([0 time_range(end) parm.fp1.freqDelimitLo-1 parm.fp1.freqDelimitHi+1])
        % Plot numbers for GPL detections
        for k = 1:length(start_sec)
            text('Units', 'normalized', 'Position', [((start_sec(k)+finish_sec(k))/2)/parm.window_length, 1.06, 0], ... 
            'String', num2str(k), 'HorizontalAlignment', 'center');
        end
        text('Units', 'normalized', 'Position', [0, 1.06, 0], ... % "Detection #" Text
            'String', 'Detection #:', 'HorizontalAlignment', 'right');
        % axis tight; 
        box on; grid off
    else  % GPL Energy Summation
        plot(time_range, log10(abs(baseline0)),'k'); % Energy Sum as black line
        for j = 1:length(start) % Plot detection regions as red line
            plot(time_range(start(j):finish(j)),log10(abs(baseline0(start(j):finish(j)))),'r');
        end
        noise_ceil_line = log10(parm.noise_ceiling*noise_floor); % Find noise_floor
        det_thresh_line = log10(parm.thresh*noise_floor); % find the noise ceiling
        plot(time_range,noise_ceil_line+0*time_range,'k--'); % plot noise floor as black dashed
        plot(time_range,det_thresh_line+0*time_range,'g--'); % plot noise ceiling as green dashed
        xlabel('Time [s]'); ylabel('Energy'); 
        t_s = sprintf('GPL Energy Summation (%i-%iHz)',parm.freq_lo,parm.freq_hi); title(t_s)
        if parm.fp1.legend == 1
            n_dum = plot(nan,nan,'-k'); g_dum = plot(nan,nan,'-r'); nc_dum = plot(nan,nan,'--k'); t_dum = plot(nan,nan,'--g');
            legend([n_dum g_dum nc_dum t_dum],'Energy Sum','GPL Detections','Noise Ceiling','Threshold','Location','northeastoutside','Orientation','vertical')
        end
        axis tight; box on; grid off
    end
    hold off

    % Plot GPL slate
    haxes = subplot(3,1,3); hold on
    imagesc(time_range, freq_range, fetched_call.gpl_contour_window);
    
    if cross_flag == 0
        [rectHandle] = det_boxes(start_sec, finish_sec, mstart_sec, mend_sec, parm, fetched_call, 0, 0, fig);
    else
        [rectHandle] = det_boxes(start_sec, finish_sec, mstart_sec, mend_sec, parm, fetched_call, cross_flag, 0, fig);
    end
    xlabel('Time [s]'); ylabel('Frequency [Hz]'); 
    t_s = sprintf('Cm: Full GPL Contours\n'); title(t_s)
    % Plot numbers for GPL detections
    for k = 1:length(start_sec)
        text('Units', 'normalized', 'Position', [((start_sec(k)+finish_sec(k))/2)/parm.window_length, 1.06, 0], ... 
        'String', num2str(k), 'HorizontalAlignment', 'center');
    end
    text('Units', 'normalized', 'Position', [0, 1.06, 0], ... % "Detection #" Text
        'String', 'Detection #:', 'HorizontalAlignment', 'right');
    text('Units', 'normalized', 'Position', [3/parm.window_length, -0.22, 0], ... % Window start time text
        'String', sprintf(datestr(window_jst+datenum([0 0 0 0 0 cross_flag*parm.window_length]))), 'HorizontalAlignment', 'right');
    axis(haxes,[0 time_range(end) parm.fp1.freqDelimitLo-1 parm.fp1.freqDelimitHi+1])
    % Adhoc indicator 
    if allowAdhocFlag == 0
         t_handle = text('Units', 'normalized', 'Position', [0.8, -0.22, 0], ... % Adhoc status
        'String', 'Adhoc Status: WAIT', 'HorizontalAlignment', 'right','BackgroundColor',[240 86 86]./255);
    else
        text('Units', 'normalized', 'Position', [0.8, -0.22, 0], ... % Adhoc status
        'String', 'Adhoc Status: Proceed', 'HorizontalAlignment', 'right','BackgroundColor',[86 240 109]./255);
    end
    % Frequency delimits
    plot3([0 parm.window_length],[parm.freq_lo parm.freq_lo],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'--w','LineWidth',1)
    plot3([0 parm.window_length],[parm.freq_hi parm.freq_hi],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'--w','LineWidth',1)   
    
    % Allocation
    fetched_call.gpl_match = [];
    fetched_call.note = cell(1,1);
    note_flag = 0;
    state(1:5) = {[]};
    endflag = [];
    userErrorFlag = [];

    % GPL Match Input
    specialModeIndicator('Mode: Quick Entry (Enter 0 for special)',1)
    specialModeIndicator('Window Flag: None',2)
    [call_id,~] = inputCallID('Enter match detection # or 0 for special case: ',GPL_struct,1);

    % Process input
    if call_id > 0 % Single detection paired 

        fetched_call.gpl_match = GPL_struct(call_id); % Write in call
        fetched_call.note = 'Single Match';
        start = GPL_struct(call_id).start_time/parm.SampleFreq;
        finish = GPL_struct(call_id).end_time/parm.SampleFreq;
        confirmBox(start, finish, parm, fetched_call, cross_flag, 1, fig) % Yellow Confirmation Box

    elseif call_id == 0 % Special Cases
        
        specialModeIndicator('Special case pairing: Idle',1)
        fprintf('\nSpecial Cases opened, enter cases, "c" to continue.\n')
        endloop = false;
        
        while ~endloop 
            
            special_ID = input('Enter special case: ','s');
            switch special_ID

                case 'r' % Regular Entry 
                    specialModeIndicator('Special case pairing: Regular entry',1)
                    [call_id,breakoutFlag] = inputCallID('Regular Entry, enter detection #: ',GPL_struct,0);
                    if breakoutFlag == 1
                        breakoutClear(1,state)
                        continue
                    end
                    [fetched_call] = removeCall(fetched_call);
                    removeBox(fig)
                    fetched_call.gpl_match = GPL_struct(call_id); % Write in call
                    fetched_call.note = {'Single Match'};
                    start = GPL_struct(call_id).start_time/parm.SampleFreq;
                    finish = GPL_struct(call_id).end_time/parm.SampleFreq;
                    confirmBox(start, finish, parm, fetched_call, cross_flag, 1, fig) % Yellow Confirmation Box
                    note_flag = 1;
                    state{1} = 'r';
                    
                case 'm' % Multiple GPL detections span single call
                    specialModeIndicator('Special Case Mode: Multiple GPL Calls',1)
                    [fetched_call] = removeCall(fetched_call);
                    [fetched_call,breakoutFlag] = multipleGPLMatch(GPL_struct,fetched_call); % Write in calls
                    if breakoutFlag == 1
                        breakoutClear(1,state)
                        continue
                    end
                    removeBox(fig) 
                    start = fetched_call.gpl_match(1).start_time / parm.SampleFreq;
                    finish = fetched_call.gpl_match(end).end_time / parm.SampleFreq;
                    confirmBox(start, finish, parm, fetched_call, cross_flag, 1, fig) % Yellow Confirmation Box
                    note_flag = 1;
                    state{1} = 'm';
                    
                case 'e' % GPL Miss: No detection is present (including no empty/poor contours)
                    specialModeIndicator('Special Case Mode: GPL Miss',1)
                    fprintf('GPL Miss: No detection is present.\n')
                    removeBox(fig) 
                    missIndicator(fig,mstart_sec,mend_sec,fetched_call,parm,1)
                    [fetched_call] = removeCall(fetched_call);
                    fetched_call.gpl_match = [];
                    fetched_call.note = {'GPL Miss'};
                    state{1} = 'e';
                    
                case 'a' % Adhoc: Allows for multiple
                    specialModeIndicator('Special Case Mode: Adhoc',3)
                    if allowAdhocFlag == 0
                        adhocWarning(fig,t_handle)
                        breakoutClear(1,state)
                        continue
                    end
                    [call_id,breakoutFlag] = inputCallID('Adhoc Detection, Enter detection number: ',GPL_struct,0);
                    if breakoutFlag == 1
                        breakoutClear(3,state)
                        continue
                    end
                    adhoc = GPL_struct(call_id);
                    adhoc.note = {'On Effort'};
                    [rectHandle] = det_boxes([], [], start_sec(call_id), finish_sec(call_id), parm, fetched_call, 0, -1, fig);
                    adhoc_detection = [adhoc_detection; adhoc];
                    note_flag = 2;
                    pause(1)
                    [strings] = IDstate(state); % Return text back to pairing pick
                    specialModeIndicator(strings,1)
                    
                case 's' % Switch: Pair manual call to a different GPL call (bc it is better) 
                    specialModeIndicator('Special Case Mode: Detection Swap',1)
                    [call_id,breakoutFlag] = inputCallID('Swap Manual detection with another GPL, Enter Detection #: ',GPL_struct,0);
                    if breakoutFlag == 1
                        breakoutClear(1,state)
                        continue
                    end
                    removeBox(fig,0) 
                    if strcmp(fetched_call.note,'False Entry')
                        falseEntry = 1;
                    else
                        falseEntry = 0;
                    end
                    [fetched_call] = removeCall(fetched_call);
                    fetched_call.gpl_match = GPL_struct(call_id); % Write in call
                    if falseEntry == 1
                        fetched_call.note = {'False Entry Swapped'};
                        missIndicator(fig,mstart_sec,mend_sec,fetched_call,parm,2)
                    else
                        fetched_call.note = {'Detection Swap'};
                        missIndicator(fig,mstart_sec,mend_sec,fetched_call,parm,3)
                    end
                    start = GPL_struct(call_id).start_time/parm.SampleFreq;
                    finish = GPL_struct(call_id).end_time/parm.SampleFreq;
                    confirmBox(start, finish, parm, fetched_call, cross_flag, 1, fig) % Yellow Confirmation Box
                    note_flag = 1;
                    state{1} = 's';
                    
                case 'f' % Omit manual detection due to false log entry 
                    specialModeIndicator('Special Case Mode: False Log Entry',1)
                    fprintf('Detection marked as false log entry.\n')
                    removeBox(fig)
                    missIndicator(fig,mstart_sec,mend_sec,fetched_call,parm,2)
                    [fetched_call] = removeCall(fetched_call);
                    fetched_call.gpl_match = [];
                    fetched_call.note = {'False Entry'};
                    state{1} = 'f';
                    
                case 'rp1' % Reprocess: GPL detector missed part of the call (after call is selected)
                    if note_flag == 0 % Check a call has been paired
                        dispErrorFlag(8)
                        continue
                    end
                    if note_flag == 1
                        fetched_call.note{end+1} = {'Reprocess: Incomplete GPL'};
                        fprintf('Pairing flagged for reprocess: GPL missed part of the call\n')
                        updatePairBox(fig,1)
                    elseif note_flag == 2
                        adhoc_detection(end).note{end+1} = {'Reprocess: Incomplete GPL'};
                        fprintf('Adhoc flagged for reprocess: GPL missed part of the call\n')
                        removeBox(fig,2,rectHandle)
                        [rectHandle] = det_boxes([], [], adhoc_detection(end).start_time/parm.SampleFreq, adhoc_detection(end).end_time/parm.SampleFreq, parm, fetched_call, 0, -3,fig);
                    end
                    if note_flag == 1
                        state{2} = 'rp1';
                    end
                    [strings] = IDstate(state);
                    specialModeIndicator(strings,1)
                    
                case 'rp2' % Reprocess: Weak/Empty Contour (after call is selected)
                    if ~isfield(fetched_call,'gpl_match') % check a call has been paired
                        dispErrorFlag(8)
                        continue
                    end
                    if note_flag == 1
                        fetched_call.note{end+1} = {'Reprocess: Weak/Empty Contour'};
                        fprintf('Pairing flagged for reprocess: Weak/Empty Contour\n')
                        updatePairBox(fig,1)
                    elseif note_flag == 2
                        adhoc_detection(end).note{end+1} = {'Reprocess: Weak/Empty Contour'};
                        fprintf('Adhoc flagged for reprocess: Weak/Empty Contour\n')
                        removeBox(fig,2,rectHandle)
                        [rectHandle] = det_boxes([], [], adhoc_detection(end).start_time/parm.SampleFreq, adhoc_detection(end).end_time/parm.SampleFreq, parm, fetched_call, 0, -3,fig);
                    end
                    state{2} = 'rp2';
                    [strings] = IDstate(state);
                    specialModeIndicator(strings,1)
                    
                case 'fn' % Flag window: GPL missed a different, good call
                    specialModeIndicator('Window Flag: False Negative',2,state)
                    fprintf('Flagging window for a seperate call that GPL missed\n')
                    window.spectrogram = fetched_call.spectrogram;
                    window.gpl_contour_window = fetched_call.gpl_contour_window;
                    window.note = {'Seperate call missed by GPL'};
                    window.win_start_time = fetched_call.window_start_time;
                    window_flag = [window_flag; window];
                    state{3}{end+1} = 'fn';
                    
                case 'new' % Flag a new call and add it to adhoc
                    specialModeIndicator('Window Flag: New Call Type',2)
                    [call_id,breakoutFlag] = inputCallID('Adhoc: New call type, Enter Detection #: ',GPL_struct,0);
                    if breakoutFlag == 1
                         breakoutClear(3,state)
                        continue
                    end
                    adhoc = GPL_struct(call_id);
                    name = input('Give the call a name, or enter "c" to continue: ','s');
                    switch name % Name the new call type
                        case 'c'
                            adhoc.note = {'Off Effort'};
                        otherwise
                            adhoc.note = {name};
                    end
                    [rectHandle] = det_boxes([], [], start_sec(call_id), finish_sec(call_id), parm, fetched_call, 0, -1, fig);
                    adhoc_detection = [adhoc_detection; adhoc];
                    
                case 'c' % Terminate special case and continue to next window
                    fprintf('Continuing to next window.\n')
                    if isempty(state{1}) % No pairing done yet
                        dispErrorFlag(21)
                        continue
                    end
                    endloop = true;
                    
                case 'x' % Flag a user error that was made
                    specialModeIndicator('Window Flag: User Error',2,state)
                    cont = false;
                    while ~cont
                        when = input('Enter "c" to mark current window for error, or "a" to mark previous window: ','s');
                        switch when
                            case 'c'
                                userErrorFlag = when;
                                cont = true;
                            case 'a'
                                userErrorFlag = when;
                                cont = true;
                            otherwise
                                dispErrorFlag(23)
                        end
                    end
                    state{3}{end+1} = 'x';
                    
                case 'k' % Acknowledge TP but skip window
                    specialModeIndicator('Special Case Mode: TP acknowledged, continuing',1)
                    removeBox(fig,1)
                    [rectHandle] = det_boxes([], [], mstart_sec, mend_sec, parm, fetched_call, 0, -2,fig);
                    fprintf('TP acknowledged, skipping window.\n')
                    [fetched_call] = removeCall(fetched_call);
                    fetched_call.note = 'TP Acknowledged';
                    endloop = true;
                    state{1} = 'k';
                    pause(1)

                case 'und' % Undeterminable due to window crossing
                    specialModeIndicator('Special Case Mode: Undetermined',1)
                    fprintf('Call marked undetermined due to window crossing.\n')
                    [fetched_call] = removeCall(fetched_call);
                    removeBox(fig)
                    updatePairBox(fig,2)
                    fetched_call.gpl_match = [];
                    fetched_call.note = 'Undetermined: Window crossing';
                    state{1} = 'und';

                case 'end' % End session
                    go = false;
                    while ~go
                        fcheck = input('Saving final pick: Are you sure you want to end? Enter "y"/"n": ','s');
                        switch fcheck
                            case 'y' % Yes, end
                                go = true;
                                term = true;
                            case 'n' % Dont end
                                go = true;
                                term = false;
                            otherwise
                                dispErrorFlag(22) % bad input
                        end
                    end
                    if term == true % End
                        if isempty(state{1}) % No pairing done yet
                            dispErrorFlag(21)
                            continue
                        end
                        endflag  = 1;
                        [fetched_call] = clearSavedWindows(fetched_call,parm);
                        return  
                    elseif term == false % Dont end
                        continue
                    end

                case 'set' % Settings
                    if parm.fp1.triton_sp == 1
                        co = false;
                        while ~co
                            ctrl = input('Enter "1" to edit the spectrogram, "2" to edit the GPL contour window: ','s');
                            switch ctrl 
                                case '1'
                                    control_tritonSp
                                    co = true;
                                case '2'
                                    [fetched_call] = contourColorAdjust(fetched_call,parm,fig,time_range,freq_range);
                                    co = true;                                 
                                otherwise
                                    dispErrorFlag(29)
                            end
                        end
                        
                    else
                        dispErrorFlag(28)
                    end
           
                otherwise
                    dispErrorFlag(9)
                    
            end % Switch: Special cases

        end % While: Special case loop

    end % If: Regular/Special case 

    hold off % Subplot 3: Contour Window
    clf(fig) % Reset figure

    % Clear out saved windows in output if necessary
    [fetched_call] = clearSavedWindows(fetched_call,parm);
    
end % Process 1


%% FP2
if process_num == 2 % Single Window Evaluation
    
    figure(8); hold on
    
    % Plot the unchanged spectrogram 
    p1 = subplot(3,1,1); hold on
    imagesc(time_range, freq_range, fetched_call.spectrogram)
    t_s1 = sprintf('Window: %s',(string(datestr(window_jst)))); suptitle(t_s1); title('Spectrogram') 
    xlabel('Time [s]'); ylabel('Frequency [Hz]')
    set(gca,'YDir','normal')
    text(-5,parm.freq_lo-6,sprintf(datestr(window_jst)))
    axis([time_range(1) time_range(end) parm.freq_lo parm.freq_hi])
    [rectHandle] = det_boxes(start_sec, finish_sec, nan, nan, parm, fetched_call, 0, 0, fig);
    hold off

    % Plot energy sum, detection thresholds, and detection regions
    p2 = subplot(3,1,2); hold on
    plot(time_range, log10(abs(baseline0)),'k'); % Energy Sum as black line
    for j = 1:length(start) % Plot detection regions as red line
        plot(time_range(start(j):finish(j)),log10(abs(baseline0(start(j):finish(j)))),'r');
    end
    noise_ceil_line = log10(parm.noise_ceiling*noise_floor); % Find noise_floor
    det_thresh_line = log10(parm.thresh*noise_floor); % find the noise ceiling
    plot(time_range,noise_ceil_line+0*time_range,'k--'); % plot noise floor as black dashed
    plot(time_range,det_thresh_line+0*time_range,'g--'); % plot noise ceiling as green dashed
    xlabel('Time [s]'); ylabel('Energy'); title('GPL Energy Summation')
    n_dum = plot(nan,nan,'-k'); g_dum = plot(nan,nan,'-r'); nc_dum = plot(nan,nan,'--k'); t_dum = plot(nan,nan,'--g');
    legend([n_dum g_dum nc_dum t_dum],'Energy Sum','GPL Detections','Noise Ceiling','Threshold','Location','northeastoutside')
    axis tight; box on; grid off
    hold off

    % Plot GPL slate
    p3 = subplot(3,1,3); hold on
    imagesc(time_range, freq_range, fetched_call.gpl_contour_window);
    xlabel('Time [s]'); ylabel('Frequency [Hz]'); 
    t_s = sprintf('GPL Contours\n'); title(t_s)
    % Plot numbers for GPL detections
    for k = 1:length(start_sec)
        text(start_sec(k), parm.freq_hi+5, num2str(k))
    end
    text(-4,parm.freq_hi+5,'Detection #:')
    text(-5,parm.freq_lo-6,sprintf(datestr(window_jst)))
    [rectHandle] = det_boxes(start_sec, finish_sec, nan, nan, parm, fetched_call, 0, 0, fig);
    axis tight
    hold off
    
    hold off
    
end % Process 2


%% FP3
if process_num == 3 % GPLReview
    
    fetched_call = rmfield(fetched_call,'sample_start');
    fetched_call = rmfield(fetched_call,'sample_end');

    figure(8); hold on
    
    % Plot the unchanged spectrogram 
    p1 = subplot(3,1,1); hold on
    imagesc(time_range, freq_range, fetched_call.whitened_sp)
    t_s1 = sprintf('Window: %s',(string(datestr(window_jst)))); suptitle(t_s1); title('Spectrogram') 
    xlabel('Time [s]'); ylabel('Frequency [Hz]')
    set(gca,'YDir','normal')
    text(-5,parm.freq_lo-6,sprintf(datestr(window_jst)))
    axis([time_range(1) time_range(end) parm.freq_lo parm.freq_hi])
    [rectHandle] = det_boxes(start_sec, finish_sec, nan, nan, parm, fetched_call, 0, 0, fig);
    hold off

    % Plot energy sum, detection thresholds, and detection regions
    p2 = subplot(3,1,2); hold on
    plot(time_range, log10(abs(baseline0)),'k'); % Energy Sum as black line
    for j = 1:length(start) % Plot detection regions as red line
        plot(time_range(start(j):finish(j)),log10(abs(baseline0(start(j):finish(j)))),'r');
    end
    noise_ceil_line = log10(parm.noise_ceiling*noise_floor); % Find noise_floor
    det_thresh_line = log10(parm.thresh*noise_floor); % find the noise ceiling
    plot(time_range,noise_ceil_line+0*time_range,'k--'); % plot noise floor as black dashed
    plot(time_range,det_thresh_line+0*time_range,'g--'); % plot noise ceiling as green dashed
    xlabel('Time [s]'); ylabel('Energy'); title('GPL Energy Summation')
    n_dum = plot(nan,nan,'-k'); g_dum = plot(nan,nan,'-r'); nc_dum = plot(nan,nan,'--k'); t_dum = plot(nan,nan,'--g');
    legend([n_dum g_dum nc_dum t_dum],'Energy Sum','GPL Detections','Noise Ceiling','Threshold','Location','northeastoutside')
    axis tight; box on; grid off
    hold off

    % Plot GPL slate
    p3 = subplot(3,1,3); hold on
    imagesc(time_range, freq_range, fetched_call.gpl_contour_window);
    xlabel('Time [s]'); ylabel('Frequency [Hz]'); 
    t_s = sprintf('GPL Contours\n'); title(t_s)
    % Plot numbers for GPL detections
    for k = 1:length(start_sec)
        text(start_sec(k), parm.freq_hi+5, num2str(k))
    end
    text(-4,parm.freq_hi+5,'Detection #:')
    text(-5,parm.freq_lo-6,sprintf(datestr(window_jst)))
    [rectHandle] = det_boxes(start_sec, finish_sec, nan, nan, parm, fetched_call, 0, 0, fig);
    axis tight
    hold off
    
    hold off
    
    % Call Selection input
    ind = 1; uselect = [];
    while true
        
        fprintf('\nEvaluating window: %s\n',string(datestr(window_jst)))
        uselect(ind) = input('Enter: ');
        if uselect(ind) == 0
            break
        end
        ind = ind + 1;
        
    end
    % Load into output
    if sum(uselect > 0) > 0
        fetched_call.calls = GPL_struct(uselect > 0);
    else
        fetched_call.calls = [];
    end
    if sum(uselect < 0) > 0
        fetched_call.FN_flag = 1;
    else
        fetched_call.FN_flag = 0;
    end
    
end % Fetch Process 3
