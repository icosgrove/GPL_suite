% Add RL measurements to Fetch outputs, for versions where the output
% doesn't contain them.
% Written: Ian 09/2024
tic
% Load Fetch output file
fname = 'bmfm_fetch_file#1.mat';
main_file = load(fname);
tf_file = '976_210816_B_HARP.tf'; % Required

% Choose xwav files
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

% Pre-processing on xwavs
field_name = cell(length(files),1);
fprintf('Loading selected xwav files...\n')
% Loop through each xwav in 'files' and extract xwav headers and audio info
for n = 1:length(files) 
    xwav_name = files(n).name; 
    PARAMS = getxwavheaders(file_dir,xwav_name); % Retrieve xwav data
    field_name{n} = sprintf('xwav%i',n);
    xwav_struct.(field_name{n}).julian_start_time = PARAMS.ltsahd.dnumStart +...
        datenum([2000,0,0,0,0,0]); % xwav start time in julian time
    audio_data = audioinfo(files(n).name);
    xwav_struct.(field_name{n}).TotalSamples = audio_data.TotalSamples; % Samples per file
    xwav_struct.(field_name{n}).nwav = PARAMS.ltsa.nrftot; % decimated xwavs per file
    xwav_start(n) = xwav_struct.(field_name{n}).julian_start_time(1);
end
fprintf('Files loaded.\n')

parm = main_file.hyd_fetch.parm;
xwav_start(1) = xwav_struct.xwav1.julian_start_time(parm.pre_offset);
parm.fp1.tf.file = tf_file;

% Primary loop over all pairings 
disp('Processing pairings...')
for n = 1:size(main_file.hyd_fetch.calls,2)

    if ~isempty(main_file.hyd_fetch.calls(n).gpl_match)
        for k = 1:size(main_file.hyd_fetch.calls(n).gpl_match,1)

            % Start/end julian time
            st_ju = main_file.hyd_fetch.calls(n).gpl_match(k).julian_start_time;
            en_ju = main_file.hyd_fetch.calls(n).gpl_match(k).julian_end_time;

            % Find xwav file 
            startLoc = find(sort([xwav_start st_ju])==st_ju);
            endLoc = find(sort([xwav_start en_ju])==en_ju);
            if startLoc ~= endLoc
                disp('Call crosses xwav, skipping') % Should be impossible
                fprintf('%i\n',n)
                continue
            else
                file_ident = startLoc-1;
                sample_range = 1:parm.nrec:parm.nrec*length(xwav_struct.(field_name{file_ident}).julian_start_time);
            end
            file_string = sprintf('xwav%i',file_ident);

            % Find window containing call
            detStart = round(seconds(datetime(st_ju,'ConvertFrom','datenum')-datetime(xwav_start(file_ident),'ConvertFrom','datenum'))*parm.SampleFreq); 
            detEnd = round(seconds(datetime(en_ju,'ConvertFrom','datenum')-datetime(xwav_start(file_ident),'ConvertFrom','datenum'))*parm.SampleFreq); 
            if file_ident == 1 % deck test offset file 1 only
                start_fix = detStart + (parm.pre_offset-1)*parm.nrec;
                end_fix = detEnd + (parm.pre_offset-1)*parm.nrec;
            else
                start_fix = detStart;
                end_fix = detEnd;
            end

            % Sort real start/end samples into sample range
            start_fix = min(find(sort([sample_range start_fix])==start_fix,1))-1;
            end_fix = min(find(sort([sample_range end_fix])==end_fix,1))-1;
            if start_fix ~= end_fix
                disp('Call crossed windows, skipping') % Should be impossible at this stage
                fprintf('%i\n',n)
                continue
            else
                start_fix = sample_range(start_fix);
                end_fix = sample_range(end_fix+1)-1;
            end

            % Do RL measurements 
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

            % Pre-processing on data chunk
            stime = main_file.hyd_fetch.calls(n).gpl_match(k).start_time;
            endtime = main_file.hyd_fetch.calls(n).gpl_match(k).end_time;
            if (flag.pre == 0) && (flag.post == 0)
                stime = stime + parm.nrec; 
                endtime = endtime + parm.nrec;
            end
            if parm.fp1.rl_time_fix == 1 % Apply time delay adjustment
                stime = stime + parm.fp1.tdelay;
                endtime = endtime + parm.fp1.tdelay;
            end
            % Apply large window that centers the call
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
            if ~isempty(main_file.hyd_fetch.calls(n).gpl_match(k).cm.values) % cm: all contours
                [p2p_RL,rms_RL,sel,peak_freq,RL_note] = computeFetchRL(main_file.hyd_fetch.calls(n).gpl_match(k).cm,parm,data,dur);
                main_file.hyd_fetch.calls(n).gpl_match(k).cm.peak_freq = peak_freq;
                main_file.hyd_fetch.calls(n).gpl_match(k).cm.p2p_RL = p2p_RL;
                main_file.hyd_fetch.calls(n).gpl_match(k).cm.rms_RL = rms_RL;
                main_file.hyd_fetch.calls(n).gpl_match(k).cm.SEL = sel;
                main_file.hyd_fetch.calls(n).gpl_match(k).cm.RL_note = RL_note;
            end
            if ~isempty(main_file.hyd_fetch.calls(n).gpl_match(k).cm_max.values) % cm_max: single strongest contour
                [p2p_RL,rms_RL,sel,peak_freq,RL_note] = computeFetchRL(main_file.hyd_fetch.calls(n).gpl_match(k).cm_max,parm,data,dur);
                main_file.hyd_fetch.calls(n).gpl_match(k).cm_max.peak_freq = peak_freq;
                main_file.hyd_fetch.calls(n).gpl_match(k).cm_max.p2p_RL = p2p_RL;
                main_file.hyd_fetch.calls(n).gpl_match(k).cm_max.rms_RL = rms_RL;
                main_file.hyd_fetch.calls(n).gpl_match(k).cm_max.SEL = sel;
                main_file.hyd_fetch.calls(n).gpl_match(k).cm_max.RL_note = RL_note;
            end
            if ~isempty(main_file.hyd_fetch.calls(n).gpl_match(k).cm_max2.values) % cm_max2: 2nd strongest contour
                [p2p_RL,rms_RL,sel,peak_freq,RL_note] = computeFetchRL(main_file.hyd_fetch.calls(n).gpl_match(k).cm_max2,parm,data,dur);
                main_file.hyd_fetch.calls(n).gpl_match(k).cm_max2.peak_freq = peak_freq;
                main_file.hyd_fetch.calls(n).gpl_match(k).cm_max2.p2p_RL = p2p_RL;
                main_file.hyd_fetch.calls(n).gpl_match(k).cm_max2.rms_RL = rms_RL;
                main_file.hyd_fetch.calls(n).gpl_match(k).cm_max2.SEL = sel;
                main_file.hyd_fetch.calls(n).gpl_match(k).cm_max2.RL_note = RL_note;
            end

        end
    end

end % Main calls

% Do adhocs
disp('Processing adhocs...')
for n = 1:length(main_file.hyd_fetch.adhoc_detections)

    % Start/end julian time
    st_ju = main_file.hyd_fetch.adhoc_detections(n).julian_start_time;
    en_ju = main_file.hyd_fetch.adhoc_detections(n).julian_end_time;

    % Find xwav file 
    startLoc = find(sort([xwav_start st_ju])==st_ju);
    endLoc = find(sort([xwav_start en_ju])==en_ju);
    if startLoc ~= endLoc
        disp('Call crosses xwav, skipping') % Should be impossible
        fprintf('%i\n',n)
        continue
    else
        file_ident = startLoc-1;
        sample_range = 1:parm.nrec:parm.nrec*length(xwav_struct.(field_name{file_ident}).julian_start_time);
    end
    file_string = sprintf('xwav%i',file_ident);

    % Find window containing call
    detStart = round(seconds(datetime(st_ju,'ConvertFrom','datenum')-datetime(xwav_start(file_ident),'ConvertFrom','datenum'))*parm.SampleFreq); 
    detEnd = round(seconds(datetime(en_ju,'ConvertFrom','datenum')-datetime(xwav_start(file_ident),'ConvertFrom','datenum'))*parm.SampleFreq); 
    if file_ident == 1 % deck test offset file 1 only
        start_fix = detStart + (parm.pre_offset-1)*parm.nrec;
        end_fix = detEnd + (parm.pre_offset-1)*parm.nrec;
    else
        start_fix = detStart;
        end_fix = detEnd;
    end

    % Sort real start/end samples into sample range
    start_fix = min(find(sort([sample_range start_fix])==start_fix,1))-1;
    end_fix = min(find(sort([sample_range end_fix])==end_fix,1))-1;
    if start_fix ~= end_fix
        disp('Call crossed windows, skipping') % Should be impossible at this stage
        fprintf('%i\n',n)
        continue
    else
        start_fix = sample_range(start_fix);
        end_fix = sample_range(end_fix+1)-1;
    end

    % Do RL measurements 
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

    % Pre-processing on data chunk
    stime = main_file.hyd_fetch.adhoc_detections(n).start_time;
    endtime = main_file.hyd_fetch.adhoc_detections(n).end_time;
    if (flag.pre == 0) && (flag.post == 0)
        stime = stime + parm.nrec; 
        endtime = endtime + parm.nrec;
    end
    if parm.fp1.rl_time_fix == 1 % Apply time delay adjustment
        stime = stime + parm.fp1.tdelay;
        endtime = endtime + parm.fp1.tdelay;
    end
    % Apply large window that centers the call
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
    if ~isempty(main_file.hyd_fetch.adhoc_detections(n).cm.values) % cm: all contours
        [p2p_RL,rms_RL,sel,peak_freq,RL_note] = computeFetchRL(main_file.hyd_fetch.adhoc_detections(n).cm,parm,data,dur);
        main_file.hyd_fetch.adhoc_detections(n).cm.peak_freq = peak_freq;
        main_file.hyd_fetch.adhoc_detections(n).cm.p2p_RL = p2p_RL;
        main_file.hyd_fetch.adhoc_detections(n).cm.rms_RL = rms_RL;
        main_file.hyd_fetch.adhoc_detections(n).cm.SEL = sel;
        main_file.hyd_fetch.adhoc_detections(n).cm.RL_note = RL_note;
    end
    if ~isempty(main_file.hyd_fetch.adhoc_detections(n).cm_max.values) % cm_max: single strongest contour
        [p2p_RL,rms_RL,sel,peak_freq,RL_note] = computeFetchRL(main_file.hyd_fetch.adhoc_detections(n).cm_max,parm,data,dur);
        main_file.hyd_fetch.adhoc_detections(n).cm_max.peak_freq = peak_freq;
        main_file.hyd_fetch.adhoc_detections(n).cm_max.p2p_RL = p2p_RL;
        main_file.hyd_fetch.adhoc_detections(n).cm_max.rms_RL = rms_RL;
        main_file.hyd_fetch.adhoc_detections(n).cm_max.SEL = sel;
        main_file.hyd_fetch.adhoc_detections(n).cm_max.RL_note = RL_note;
    end
    if ~isempty(main_file.hyd_fetch.adhoc_detections(n).cm_max2.values) % cm_max2: 2nd strongest contour
        [p2p_RL,rms_RL,sel,peak_freq,RL_note] = computeFetchRL(main_file.hyd_fetch.adhoc_detections(n).cm_max2,parm,data,dur);
        main_file.hyd_fetch.adhoc_detections(n).cm_max2.peak_freq = peak_freq;
        main_file.hyd_fetch.adhoc_detections(n).cm_max2.p2p_RL = p2p_RL;
        main_file.hyd_fetch.adhoc_detections(n).cm_max2.rms_RL = rms_RL;
        main_file.hyd_fetch.adhoc_detections(n).cm_max2.SEL = sel;
        main_file.hyd_fetch.adhoc_detections(n).cm_max2.RL_note = RL_note;
    end

end

% Save back to file
hyd_fetch = main_file.hyd_fetch;
save(fname,'hyd_fetch')
toc