%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPL_fetch is a process to extract GPL detections from a set of manual
% detections. The process will fetch the correct window that a manual
% detection falls within, apply the GPL algorithm to that window and allow
% the user to make a decision on which GPL detection number corresponds to
% their manual pick. There are several special options to handle different
% cases that may arise, see the readme for more information. 
%
% GPL_parameter_input_v3.m should be opened and configured to reflect the
% parameters that would be used to detect the signal of interest that will
% be fetched.
%
% A manual detection log must be loaded for this process, in the form of
% .xlsx file that is created by the Logger Remora on Triton. 
%
%% 
% Written by Ian Cosgrove, using GPL algorithm functions developed by Tyler
% Helble, and advised by Joshua Jones, Last updated 07/19/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clc


%% User setup

% Pair detlog detections
fetch_p1 = 1; % Manual Detection Loop (Off:0,On:1)
fetch_p1_detlog = 'NUNAT_SB_03_TritonLog_LowFreq_KH_30-14 Hz down.xlsx'; % Detlog .xlsx file

% Select Transfer function to load. Comment out if not used. 
tf_file = '976_210816_B_HARP.tf';

% If you have a fetch parameter file ready to go enter it here. If not
% leave empty ('')
parm_file = '';




% Single Window evaluation
fetch_p2 = 0; % Single Window Extraction (Off:0,On:1)
fetch_p2_datetime = '08-11-2021 19:58:02.0'; % Form: 'MM-DD-YYYY HH:MM:SS.m OR datenum

% GPL Review v2 (recommended to use original GPLReview)
fetch_p3 = 0; % GPL Detection Review (Off:0,On:1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load parameter file and xwavs - Tyler Helble

% Load xwav files 
fprintf('\nSelect folder containing .x.wav that pertain to manual detections.\n');
    [filename, pathname] = uigetfile('*.x.wav'); % Select File
    if filename == 0
        dispErrorFlag(17)
        return
    end
    cwd = pwd; 
    cd(pathname) % Set current directory to path containing xwav file
    file_dir = pwd; 
    addpath(pwd); 
    files = dir('*.x.wav'); % Variable files contains xwav data
    cd(cwd); % Set current directory back to current working directory
% End Helble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off')

%% Extract xwav headers - Ian Cosgrove

% xwav name handle
field_name = cell(length(files),1);

fprintf('Loading selected xwav files...\n')

% Loop through each xwav in 'files' and extract xwav headers and audio info
for n = 1:length(files) 
    xwav_name = files(n).name; 
    PARAMS = getxwavheaders(file_dir,xwav_name); % Retrieve xwav data
    field_name{n} = sprintf('xwav%i',n);
    xwav_struct.(field_name{n}).julian_start_time = PARAMS.ltsahd.dnumStart +...
        datenum([2000,0,0,0,0,0]); % xwav start time in julian time
%     xwav_struct.(field_name{n}).year = PARAMS.ltsahd.year; % date string for start times
%     xwav_struct.(field_name{n}).month = PARAMS.ltsahd.month; % " "
%     xwav_struct.(field_name{n}).day = PARAMS.ltsahd.day; % " "
%     xwav_struct.(field_name{n}).hour = PARAMS.ltsahd.hour; % " "
%     xwav_struct.(field_name{n}).minute = PARAMS.ltsahd.minute; % " "
%     xwav_struct.(field_name{n}).second = PARAMS.ltsahd.secs; % " "
%     xwav_struct.(field_name{n}).sample_freq = PARAMS.ltsahd.sample_rate; % xwav sampel frequency
%     xwav_struct.(field_name{n}).file_name = PARAMS.ltsahd.fname; % xwav origin filename
    audio_data = audioinfo(files(n).name);
    xwav_struct.(field_name{n}).TotalSamples = audio_data.TotalSamples; % Samples per file
    xwav_struct.(field_name{n}).nwav = PARAMS.ltsa.nrftot; % decimated xwavs per file
    xwav_start(n) = xwav_struct.(field_name{n}).julian_start_time(1);
end

fprintf('Files loaded.\n')

% Load manual parameter file or run parameter setup
if ~isempty(parm_file)
    load(parm_file)
else
    parm = [];
    [parm] = GPL_fetch_parameter_input(PARAMS,xwav_struct,parm);
end

% Load TF file if not done so already
if ~isfield(parm.fp1.tf,'file')
    parm.fp1.tf.file = tf_file;
end

% Adjust the start time to after deck test
xwav_start(1) = xwav_struct.xwav1.julian_start_time(parm.pre_offset);

% Check for timestamp issue
if ~isempty(find(diff(diff(xwav_struct.xwav1.julian_start_time)) > 0.1,1))
    fprintf('\nWarning: Possible timestamp inconsistency in xwav headers (deck test or other jump in time).\n\n')
end

% Configure time bins
if parm.fetch.time_delay_switch == 0
    parm.NumTimeBins = parm.NumTimeBins - parm.fftl/parm.fftOverlap + 1;
end



%% Fetch Process 1: Extract All windows from an exisitng Detlog

if fetch_p1 == 1
    
    try % Load detlog with path error catch
        manual_log = readtable(fetch_p1_detlog);
    catch ME
        if (strcmp(ME.identifier,'MATLAB:readtable:OpenFailed')) || (strcmp(ME.identifier,'MATLAB:textio:textio:FileNotFound'))
            disp('Error: Cannot locate Detlog file, ensure "detlogs" folder is added to path and try again.');
            return
        else
            rethrow(ME);
        end
    end

    % Request to keep working on workbook or start over
    proceed = false; go = false;
    while ~proceed
        fselect = input('Enter "n" for a new GPL_Fetch effort or "c" to continue an existing file: ', 's');
        if strcmp(fselect,'c') 
            proceed = true;
            fprintf('Continuing existing file, please select .mat file.\n')
            [file,path] = uigetfile('*.mat');
            if isequal(file,0) % Canceled selection
                dispErrorFlag(17)
                return
            else
                workingFile = load(fullfile(path,file));
                latest_pick = length(workingFile.hyd_fetch.calls);

                % Adjust latest pick back to correct start point
                lpick = workingFile.hyd_fetch.calls(end).manual_start_time;
                for k = 1:length(manual_log.StartTime)
                    if strcmp(datestr(manual_log.StartTime(k),'yyyy-mm-dd HH:MM:SS.FFF'),lpick)
                        latest_pick = k;
                    end
                end

                total_calls = length(manual_log.InputFile);
                while ~go
                    resetstring = sprintf('Enter "c" to continue on the SAME output file (Pick %i/%i) or "new" to start on a NEW output file: ',latest_pick+1,total_calls);
                    resetcheck = input(resetstring,'s');
                    if strcmp(resetcheck,'c') % continue from last pick
                        go = true;
                        hyd_fetch = workingFile.hyd_fetch;
                    elseif strcmp(resetcheck,'new') % reset
                        proceed = true;
                        file = input('Enter new output filename (exclude ".mat"): ','s');
                        file = strcat(file,'.mat');
                        hyd_fetch(1).calls = [];
                        hyd_fetch(1).adhoc_detections = [];
                        hyd_fetch(1).window_flag = [];
                        hyd_fetch(1).parm = parm;
                        hyd_fetch(1).manual_log = [];
                        go = true;
                    else
                        dispErrorFlag(19)
                    end
                end
            end
        elseif strcmp(fselect,'n') % new effort
            proceed = true;
            fprintf('Beginning a new file.\n')
            file = input('Enter output filename (exclude ".mat"): ','s');
            file = strcat(file,'.mat');
            hyd_fetch(1).calls = [];
            hyd_fetch(1).adhoc_detections = [];
            hyd_fetch(1).window_flag = [];
            hyd_fetch(1).parm = parm;
            hyd_fetch(1).manual_log = [];
            latest_pick = 0;
        else
            dispErrorFlag(14)
        end % If
    end % While

    % Allocation
    calls_append = [];
    fetched_callredo = [];
    adhoc_append = [];
    window_flag_append = [];
    endflag = [];
    dontEndFlag = 0;
    note = [];

    % Sorting detections
    mst_julian = sort(datenum(manual_log.StartTime)); % Start times 
    met_julian = sort(datenum(manual_log.EndTime)); % End times

    % Make sure the log isnt already done
    if (latest_pick + 1) > length(mst_julian)
        dispErrorFlag(20)
        return
    else % Text indicating decision
        if strcmp(fselect,'c') 
            if strcmp(resetcheck,'c')
                fprintf('Resuming on pick %i.\n',latest_pick+1)
            elseif strcmp(resetcheck,'new')
                fprintf('Starting a new file on pick %i.\n',latest_pick+1)
            end
        elseif strcmp(fselect,'n') 
            fprintf('Starting new pairing file.\n')
        end
    end
    
    % Loop over manual start times and process relevant window(s)
    for n = (latest_pick+1):length(mst_julian)

        % Check for warning of incorrect xwavs selected
        [dec] = xwavSelectWarning(xwav_start,mst_julian(n));
        if strcmp(dec,'b') % Break out if requested
            return
        end

        % Locate sample range that contains manual detection
        [SampleRangeStart, SampleRangeEnd, error_flag, file_string, sample_range, detStart, detEnd, ...
            detStart_rel, detEnd_rel, file_ident] = fetchSort(mst_julian, met_julian, n, xwav_start, parm, xwav_struct);

        % Identify if adhoc needs to wait for other calls in same window
        [allowAdhocFlag] = allowAdhoc(mst_julian,met_julian,sample_range,SampleRangeStart,SampleRangeEnd,xwav_start,parm,n,file_ident);
        
        % Create sample vector for xwav timestamps
        sample_range = 1:parm.nrec:parm.nrec*length(xwav_struct.(field_name{file_ident}).julian_start_time);

        if error_flag == 4 % Nearly impossible case: detection crosses xwav files
            dispErrorFlag(error_flag)
            continue
        end
        
        % Load audio data for window(s) containing detection.
        if SampleRangeStart == SampleRangeEnd

            % Account for beginning/end file cases
            [start_fix, end_fix] = sampleRangeFix(sample_range, SampleRangeStart, parm, xwav_struct, file_string);
            
            if end_fix > xwav_struct.(field_name{file_ident}).TotalSamples
                disp('Error: Over-indexed, continuing')
                continue
            end

            sub_data = double(audioread(files(file_ident).name,...
                [start_fix,end_fix],'native')); % Account for overlap for full window size
            
            parm = timeDelayConfig(parm,SampleRangeStart);
            
            fprintf('\nManual Detection: %s\n',datestr(mst_julian(n)))
            fprintf('Processing Window: %s\n',datestr(xwav_struct.(field_name{file_ident}).julian_start_time(SampleRangeStart)))

            [GPL_struct,fetched_call,error_flag,adhoc_detection,window_flag,endflag,userErrorFlag,parm] = GPL_v3_fetch(sub_data,detStart,detEnd,parm,xwav_struct.(file_string).julian_start_time(SampleRangeStart),mst_julian(n),met_julian(n),n,0,1,manual_log,allowAdhocFlag);
            
            % Check compatibility, make them redo if doesn't fit.
            if ~isempty(calls_append)
                areCompat = structcompatable(fetched_call,calls_append);
                if areCompat ~= 1
                     dispErrorFlag(18)
                    GPL_struct = []; fetched_call = [];
                    [GPL_struct,fetched_call,error_flag,adhoc_detection,window_flag,endflag,userErrorFlag,parm] = GPL_v3_fetch(sub_data,detStart,detEnd,parm,xwav_struct.(file_string).julian_start_time(SampleRangeStart),mst_julian(n),met_julian(n),n,0,1,manual_log,allowAdhocFlag);
                end
            end
            
            % Try again for error cases
            while error_flag ~= 0
                dispErrorFlag(error_flag)
                if error_flag == 3 % Redo previous call
                    fprintf('\nRe-doing the previous call\n')
                    calls_append(n-1) = []; % Remove previous call
                    [SampleRangeStart, SampleRangeEnd, ~, file_string, sample_range, detStart, detEnd, ...
                        detStart_rel, detEnd_rel, file_ident] = fetchSort(mst_julian, met_julian, n-1, xwav_start, parm, xwav_struct);
                    [start_fix, end_fix] = sampleRangeFix(sample_range, SampleRangeStart, parm, xwav_struct, file_string);
                    sub_data = audioread(files(file_ident).name,[start_fix,end_fix]);
                    [GPL_struct,fetched_call,error_flag,adhoc_detection,window_flag,endflag,userErrorFlag,parm] = GPL_v3_fetch(sub_data,detStart,detEnd,parm,xwav_struct.(file_string).julian_start_time(SampleRangeStart),mst_julian(n-1),met_julian(n-1),n-1,0,1,manual_log,allowAdhocFlag);
                else % Other error redo cases
                    GPL_struct = []; fetched_call = [];
                    [GPL_struct,fetched_call,error_flag,adhoc_detection,window_flag,endflag,userErrorFlag,parm] = GPL_v3_fetch(sub_data,detStart,detEnd,parm,xwav_struct.(file_string).julian_start_time(SampleRangeStart),mst_julian(n),met_julian(n),n,0,1,manual_log,allowAdhocFlag);
                end
            end
            
            % Load into primary output struct
            if isempty(fetched_callredo) % No re-do
                calls_append = [calls_append fetched_call];
                [adhoc_append] = appendAdhocDetection(adhoc_append,adhoc_detection);
                window_flag_append = [window_flag_append window_flag];
            else % Include re-do
                calls_append = [calls_append fetched_callredo fetched_call];
                fetched_callredo = [];
            end
            

            
%             % Re-do matching if user needs to
%             if redo_flag == 1
%                 GPL_struct = []; fetched_call = [];
%                 [GPL_struct,fetched_call,redo_flag] = GPL_v3_fetch(sub_data,detStart,detEnd,parm,xwav_struct.(file_string).julian_start_time(SampleRangeStart),mst_julian(n),met_julian(n),n,0,1,manual_log);
%             end
            

        else % Detection crosses between windows

            clear sub_data
            calls_TEMP = []; window_TEMP = [];
            for k = 1:(SampleRangeEnd - SampleRangeStart + 1)
                fprintf('\nProcessing: %s\n',datestr(xwav_struct.(field_name{file_ident}).julian_start_time(SampleRangeStart+k-1)))
                parm = timeDelayConfig(parm,SampleRangeStart);
                if parm.fetch.time_delay_switch == 1
                    if parm.fp1.triton_sp == 1 % Triton sp, using future/past samples (likely unused)
                        sub_data = double(audioread(files(file_ident).name,[sample_range(SampleRangeStart)-(parm.time_delay*parm.SampleFreq)+(k-1)*parm.nrec, ...
                            sample_range(SampleRangeStart)-(parm.time_delay*parm.SampleFreq)+k*parm.nrec],'native')); 
                    else % Old sp with future/past samples
                        sub_data = double(audioread(files(file_ident).name,[sample_range(SampleRangeStart)-(parm.time_delay*parm.SampleFreq)+(k-1)*parm.nrec, ...
                            sample_range(SampleRangeStart)-(parm.time_delay*parm.SampleFreq)+parm.fftOverlapPercent*0.01*parm.fftl+k*parm.nrec],'native')); 
                    end
                else
                    if parm.fp1.triton_sp == 1 % Triton sp with no future/past samples (most likely used) 
                        sub_data = double(audioread(files(file_ident).name,[sample_range(SampleRangeStart)+(k-1)*parm.nrec, ...
                            sample_range(SampleRangeStart)+k*parm.nrec-1],'native'));
                    else % Old sp with no future/past samples
                        sub_data = double(audioread(files(file_ident).name,[sample_range(SampleRangeStart)-(parm.SampleFreq)+(k-1)*parm.nrec, ...
                            sample_range(SampleRangeStart)-(parm.SampleFreq)+parm.fftOverlapPercent*0.01*parm.fftl+k*parm.nrec],'native')); 
                    end
                end
                [GPL_struct,fetched_call,error_flag,adhoc_detection,window_flag,endflag,userErrorFlag,parm] = GPL_v3_fetch(sub_data,detStart,detEnd,parm,xwav_struct.(file_string).julian_start_time(SampleRangeStart),mst_julian(n),met_julian(n),n,k-1,1,manual_log,allowAdhocFlag);
                calls_TEMP = [calls_TEMP fetched_call];
                [adhoc_append] = appendAdhocDetection(adhoc_append,adhoc_detection); % Adhoc added as usual
                window_TEMP = [window_TEMP window_flag];
                if k < (SampleRangeEnd - SampleRangeStart + 1)
                    disp('Please wait until window overlap case is handled before ending.')
                end
            end
            [call_ready, window_ready] = prepareWindowOutput(calls_TEMP,window_TEMP);
            calls_append = [calls_append call_ready];
            window_flag_append = [window_flag_append window_ready];
        end

        clear note
        % Load user error flags into either current or previous call. 
        if ~isempty(userErrorFlag)
            switch userErrorFlag
                case 'c'
                    calls_append.note{end+1} = 'User Error';
                    disp('Current window marked for error.')
                case 'a'
                    if ~isempty(hyd_fetch.calls)
                        note = hyd_fetch(1).calls(end).note;
                        note = {note 'User Error'};
                        hyd_fetch(1).calls(end).note = note;
                        disp('Previous window marked for error.')
                    end
            end
        end

        % Extract received levels for each pairing and adhoc
        flag.post = 0; flag.pre = 0;
        if (start_fix - parm.nrec) > 1 % Add padding except first/last window
            if (end_fix + parm.nrec) < xwav_struct.(file_string).TotalSamples
                rl_start = start_fix - parm.nrec; 
                rl_end = end_fix + parm.nrec;
            else
                flag.post = 1;
            end
        else
            flag.pre = 1;
        end
        rl_data = double(audioread(files(file_ident).name,...
                [rl_start,rl_end],'native'));
        [calls_append, adhoc_append] = estFetchRL(calls_append,adhoc_append,parm,rl_data,flag);

        % Save each output to the file as they are made
        hyd_fetch(1).calls = [hyd_fetch(1).calls calls_append];
        hyd_fetch(1).adhoc_detections = [hyd_fetch(1).adhoc_detections; adhoc_append];
        hyd_fetch(1).window_flag = [hyd_fetch(1).window_flag; window_flag_append];
        if isempty(hyd_fetch.parm)
            hyd_fetch(1).parm = parm;
        end
        if isempty(hyd_fetch.manual_log)
            hyd_fetch(1).manual_log = manual_log;
        end
        save(file,'hyd_fetch')

        if endflag == 1
            fprintf('Ending session at pick %i.\n',n)
            return
        end

        % Reset since previous is already saved
        calls_append = [];
        adhoc_append = [];
        window_flag_append = [];

        % Output verification
        if parm.fp1.outputVerify == 1
            verifyGPLFetchOutput(hyd_fetch,n)
        end
 
    end % Loop over all MST
    
    fprintf('\nMatching complete.\n')

end % Fetch Process 1




















%% Fetch Process 2: Single Window Evaluation

if fetch_p2 == 1

    % Locate sample range of requested window
    if ischar(fetch_p2_datetime) || isstring(fetch_p2_datetime)
        time_req = datenum(fetch_p2_datetime);
    else % already datenum
        time_req = fetch_p2_datetime;
    end

    % Sort request into an xwav
    Sort = sort([time_req xwav_start(1,:)]); % Sort
    FileLoc = find(Sort == time_req,1,'last')-1; % Locate index

    % Create sample range possibilities for the file
    file_string = sprintf('xwav%i',FileLoc);
    sample_range = 1:parm.nrec:xwav_struct.(file_string).TotalSamples;
    if sample_range(end) < xwav_struct.(file_string).TotalSamples
        sample_range(end+1) = xwav_struct.(file_string).TotalSamples; % Assuming total samples is int. multiple of nrec
    end
    
    % Locate sample range of request
    req_start = round(seconds(datetime(time_req,'ConvertFrom','datenum')-datetime(xwav_start(FileLoc),'ConvertFrom','datenum'))*parm.SampleFreq); 
    if req_start == 0 % Start of xwav
        req_start = 1;
    end
    if parm.fetch.time_delay_switch == 1
        parm = timeDelayConfig(parm,req_start);
    end
    SampleRangeindex = find(sort([sample_range req_start]) == req_start,1,'last' );
    if SampleRangeindex == 1 % Start of xwav
        start = 1;
        finish = parm.nrec+1e-2*parm.fftOverlapPercent*parm.fftl-1;
    elseif SampleRangeindex == length(sample_range)-1 % End of xwav
        start = sample_range(SampleRangeindex);
        finish = xwav_struct.(file_string).TotalSamples;
    else % Everything else
        start = sample_range(SampleRangeindex)-(parm.time_delay_switch*parm.SampleFreq);
        finish = start+parm.nrec+parm.fftOverlapPercent*1e-2*parm.fftl;
    end
    
    sub_data = audioread(files(FileLoc).name,[start finish]);
    [GPL_struct,fetched_call] = GPL_v3_fetch(sub_data,nan,nan,parm,xwav_struct.(file_string).julian_start_time(SampleRangeindex),nan,nan,nan,nan,2,nan);
    
end



%% Fetch Process 3: GPL Detection Evaluation (Incomplete)

if fetch_p3 == 1
    
    % Load detections
    fprintf('\nSelect folder containing detection .mat files to process\n');
    [filename, pathname] = uigetfile('*.mat'); % Select File
    cwd = pwd; 
    cd(pathname) % Set current directory to path containing xwav file
    file_dir = pwd; 
    addpath(pwd); 
    detfiles = dir('*.mat'); % Variable files contains xwav data
    cd(cwd); % Set current directory back to current working directory
    
    tp_calls = [];
    fn_flag = 0;
    
    % Sort .mat files by number
    order = 1;
    k = 0;
    if length(detfiles) > 1
        for n = 1:length(detfiles)
            fname = detfiles(n).name;
            while ~isnan(str2double(string(fname(end-4-k:end-4))))
                order(n) = str2double(string(fname(end-4-k:end-4)));
                k = k + 1;
            end
            k = 0;
        end
    end
    
    % Loop through detections files
    for a1 = 1:length(detfiles)
        
        % Sort file number correctly
        corr_fname = detfiles(order == a1);
        fname = sprintf('%s',corr_fname.name);
        fprintf('Loading File %i ...\n',a1); 
        load(fname)
        fprintf('File successfully loaded\n');
        
        start_time = [hyd.detection.calls.start_time];
        end_time = [hyd.detection.calls.end_time];
        
        % Isolate calls by xwav file
        for a2 = 1:length(files)
            
            current_calls = find(strcmp({hyd.detection.calls(:).fname}, files(a2).name));
            
            % Create sample range possibilities for the file
            file_string = sprintf('xwav%i',a2);
            sample_range = 1:parm.nrec:xwav_struct.(file_string).TotalSamples;
            if sample_range(end) < xwav_struct.(file_string).TotalSamples
                sample_range(end+1) = xwav_struct.(file_string).TotalSamples; % Assuming total samples is int. multiple of nrec
            end
            
            fprintf('\nOverview for evaluation:\nEnter Detection # for TP GPL detections (e.g. 3)\nEnter "0" for no TP/to break out of a window\nEnter "-1" (before entering 0) to flag a GPL missed detection\n\n')
            
            % Isolate & Process calls by sample range window
            if ~isempty(current_calls)
                
                for a3 = 1:length(sample_range)-1
                    
                    win_start = find(start_time >= sample_range(a3) & start_time < sample_range(a3 + 1));
                    win_end = find(end_time >= sample_range(a3) & end_time < sample_range(a3 + 1));

                    if isequal(win_start, win_end) % No overlap
                        if ~isempty(win_start) || ~isempty(win_end)
                        
                            if a3 == 1 % Start of xwav
                                start = 1;
                                finish = parm.nrec+1e-2*parm.fftOverlapPercent*parm.fftl;
                            elseif a3 == length(a3)-1 % End of xwav
                                start = sample_range(a3);
                                finish = xwav_struct.(file_string).TotalSamples;
                            else % Everything else
                                start = sample_range(a3)-(parm.time_delay*parm.SampleFreq);
                                finish = start+parm.nrec+parm.fftOverlapPercent*1e-2*parm.fftl;
                            end
                            
                            % Isolate calls for the given window and create
                            % plot
                            current_window = hyd.detection.calls(win_start);
                            sub_data = audioread(files(a2).name,[start finish]);
                            [~,fetched_window] = GPL_v3_fetch(sub_data,nan,nan,parm,xwav_struct.(file_string).julian_start_time(a3),nan,nan,nan,0,3,nan);
                            for a4 = 1:length(fetched_window.calls)
                                fetched_window.calls(a4).xwav_fname = files(a2).name;
                                fetched_window.calls(a4).det_fname = fname;
                            end              
                        end
                        
                    else % Window Overlap
                        
                        
                        
                    end
                        
                    tp_calls = [tp_calls fetched_window.calls]; 
                    if fetched_window.FN_flag == 1
                        fn_flag = fn_flag + 1;
                        fn_report(fn_flag).fn_window = fetched_call.window_start_time;
                    end

                end 
                    
                
            end % If: Processing calls for current window
        end % For: Loop over xwav files
        
            
    end % For: Loop over all det files
        
    hyd_tp(1).calls = tp_calls;
    hyd_tp(1).parm = hyd.detection.parm;
    
end % If: Run Fetch Process 3
