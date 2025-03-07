%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process_HARP_v3 is the main script that will run the Generalized Power
% Law (GPL) Detector Algorithm on a set of acoustic data. The algorithm
% will parse through provided xwavs and perform statistical cmputations to
% extract cetacean vocalizations, or other non-biological signals of
% interest, from the data.

% This version was written by Ian Cosgrove, advised by Joshua Jones, starting
% 01/12/2024. The original versions (GPL_v2 (2020), GPL 2014) were created
% by Tyler Helble. See his publication: "A generalized power-law detection 
% algorithm for humpback whale vocalizations" (2012)
% (https://www.cetus.ucsd.edu/publications.html) for full documentation on
% the theory behind the GPL detector. 

% This script will perform the functions of loading the xwavs that are to
% be processed by the detector. The xwav headers will be extracted and the
% total amount of time will be separated into samller windows, these 
% windows will be processed one at a time in the script GPL_v3. After
% GPL_v3 and it's subfunctions have completed the detection process on a
% window, the output data will return to this script as GPL_struct, and it
% will be saved to the main output file. Once all windows have been
% processed, the final detection output file will be saved from here. 

% Documentation written by Ian Cosgrove
% Contact: icosgrove@ucsd.edu
% Last Updated: 03/13/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

total_elapsed_time = tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load parameter file and xwavs - Tyler Helble

% Load xwav files 
fprintf('\nSelect folder containing .x.wav to process\n');
    [filename, pathname] = uigetfile('*.x.wav'); % Select File
    cwd = pwd; 
    cd(pathname) % Set current directory to path containing xwav file
    file_dir = pwd; 
    addpath(pwd); 
    files = dir('*.x.wav'); % Variable files contains xwav data
    cd(cwd); % Set current directory back to current working directory


% Allocate matrix for detections
calls = [];


% End Helble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract xwav headers - Ian Cosgrove

% Pre-allocate
field_name = cell(length(files),1);

% Loop through each xwav in 'files' and extract xwav headers and audio info
for n = 1:length(files) 
    xwav_name = files(n).name; 
    PARAMS = getxwavheaders(file_dir,xwav_name); % Retrieve xwav data
    field_name{n} = sprintf('xwav%i',n);
    xwav_struct.(field_name{n}).julian_start_time = PARAMS.ltsahd.dnumStart; % xwav start time in julian time
    xwav_struct.(field_name{n}).year = PARAMS.ltsahd.year; % date string for start times
    xwav_struct.(field_name{n}).month = PARAMS.ltsahd.month; % " "
    xwav_struct.(field_name{n}).day = PARAMS.ltsahd.day; % " "
    xwav_struct.(field_name{n}).hour = PARAMS.ltsahd.hour; % " "
    xwav_struct.(field_name{n}).minute = PARAMS.ltsahd.minute; % " "
    xwav_struct.(field_name{n}).second = PARAMS.ltsahd.secs; % " "
    xwav_struct.(field_name{n}).sample_freq = PARAMS.ltsahd.sample_rate; % xwav sampel frequency
    xwav_struct.(field_name{n}).file_name = PARAMS.ltsahd.fname; % xwav origin filename
    audio_data = audioinfo(files(n).name);
    xwav_struct.(field_name{n}).TotalSamples = audio_data.TotalSamples; % Samples per file
    xwav_struct.(field_name{n}).nwav = PARAMS.ltsa.nrftot; % decimated xwavs per file
end


% Run Parameter setup
clear parm
[parm] = GPL_parameter_input_v3(PARAMS,xwav_struct);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear n

% Counter to separate output files once they reach a certain size
% Note that when naming the file, the number of the file must be the final 
% character before '.mat'. Example: det_sb_01_file1.mat, det_sb_01_file2.mat, etc...
if parm.output_size_separation == 1
    file_sep_counter = 0;
end

% Counter to track number of detections per window for plotting later.
if parm.det_frequency_plot == 1
    call_freq_counter = 0;
end

% Counter for tracking runtime through each window
if parm.call_time_tracker == 1
    time_counter = 0;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Primary Loop - GPL Algorithm
% Ian Cosgrove

if parm.pad == 0 % Run consecutive, even windows without padding.

    % Loop through all xwav files in selected folder.
    for n = 1:length(files)

        % Calculate number of sub-loops through each file. (Taken either by
        % how xwavs are broken up or calculation from param. input)
        if xwav_struct.(field_name{n}).nwav == 1
            num_wav_loops = xwav_struct.(field_name{n}).TotalSamples / parm.nrec;
            single_wav_flag = 1;
        else
            num_wav_loops = xwav_struct.(field_name{n}).nwav;
            single_wav_flag = 0;
        end

        temp_name = files(n).name;
        PARAMS = getxwavheaders(file_dir,temp_name); 
        
        first_call_flag = 1;

        % Find deck test jumps
        if parm.deck_test_adjustment == 1
            j_flag = findTimeStepJump(PARAMS,parm,num_wav_loops,parm.nrec);
            jst = xwav_struct.(field_name{n}).julian_start_time(1);
        end
        
        % Loop through all disk writes within current xwav
        for q = 1:num_wav_loops

            % Display current window being processed
            if single_wav_flag == 0
                fprintf('\n\nProcessing: %s\n',datestr(xwav_struct.(field_name{n}).julian_start_time(q) + ...
                    datenum([2000,0,0,0,0,0])));
            else
                fprintf('\n\nProcessing: %s\n',datestr(xwav_struct.(field_name{n}).julian_start_time + ...
                    datenum([2000,0,0,0,0,(q-1)*(parm.nrec/parm.SampleFreq)])));
            end
            
            
            % Determine sample range of current window
            if parm.future_sample_ext == 1
                SampleRangeStart = (parm.nrec*(q-1)) + 1;
                if single_wav_flag == 0
                    
                    if q ~= xwav_struct.(field_name{n}).nwav
                        SampleRangeEnd = parm.nrec*(q) + 1 + parm.fftOverlapPercent*1e-2*parm.fftl;
                        timeBinReset_flag = 0;
                    else
                        SampleRangeEnd = parm.nrec*q; % Ensure requested samples ~> Total samples
                        parm.NumTimeBins = (parm.nrec-parm.fftl)/parm.fftOverlap +1;
                        timeBinReset_flag = 1;
                    end
                    
                else % Single .wav case
                    
                    if q ~= num_wav_loops
                        SampleRangeEnd = parm.nrec*(q) + 1;
                    else
                        SampleRangeEnd = parm.nrec*q; % Ensure requested samples ~> Total samples
                    end
                    
                end

            else % Tyler's sample extraction

                sample_padding = parm.pad*parm.SampleFreq; 
                num_new_samples = parm.nrec - 2*sample_padding;
                offset = (q-1)*num_new_samples + sample_padding;
                SampleRangeStart = offset + 1 - sample_padding;
                SampleRangeEnd = offset + 1 + sample_padding + num_new_samples;

            end
            if parm.deck_test_adjustment == 1
                if SampleRangeStart >= j_flag
                    jst = xwav_struct.(field_name{n}).julian_start_time(q);
                end
            else
                jst = xwav_struct.(field_name{n}).julian_start_time(q);
            end
                
            % Extract audio data for the current window
            sub_data = audioread(files(n).name,[SampleRangeStart,SampleRangeEnd]);
            %%%%%%%%%%%%%%%%%%% audioread autonormalize fix: 'native' %%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Error check for corrupted HARP data - Tyler Helble
            error_check = sort(abs(sub_data)); 
            error_check = mode(error_check);
            if error_check == 0
                fileID = fopen('errors.txt','at'); % Create text file
                fprintf(fileID, files(q).name,'\n'); % File Location
                fprintf(fileID, datestr(xwav_struct.(field_name{n}).julian_start_time),'\n');
                fclose(fileID); % Close text file
                continue
            end

            % End Helble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Begin timer for individual window processing  
            if parm.call_time_tracker == 1
                v3_runtime = tic;
            end
            

            % Process current window through GPL algoirthm: identify calls,
            % extract call contours, take contour measurements, and perform
            % preliminary filtering.
            [GPL_struct] = GPL_v3(sub_data,parm);
            
            % Record the window number of first detection
            if first_call_flag == 1
                if ~isempty(GPL_struct)
                    if single_wav_flag == 0
                        hyd(1).detection.start.window = xwav_struct.(field_name{n}).julian_start_time(q) + ...
                            datenum([2000,0,0,0,0,0]);
                        first_call_flag = 0;
                    else
                        hyd(1).detection.start.window = xwav_struct.(field_name{n}).julian_start_time + ...
                            datenum([2000,0,0,0,0,0]);
                        first_call_flag = 0;
                    end
                end
            end
            
            % End timer for individual window and load it to output
            if parm.call_time_tracker == 1
                time_counter = time_counter + 1;
                call_time_array(time_counter) = toc(v3_runtime);
            end
            
            % Loop through detections just found, update detection times
            % relative to start of the xwav.

            for k = 1:length(GPL_struct)
                
                % Update sample time to be correct count from start of xwav
                GPL_struct(k).start_time = GPL_struct(k).start_time + ((q - 1)*parm.nrec);
                GPL_struct(k).end_time = GPL_struct(k).end_time + ((q - 1)*parm.nrec);
                
                % Update Julian time relative to start of xwav
                GPL_struct(k).julian_start_time = jst + ...
                    datenum(2000,0,0,0,0,GPL_struct(k).start_time/parm.SampleFreq);
                GPL_struct(k).julian_end_time = jst +...
                    datenum(2000,0,0,0,0,GPL_struct(k).end_time/parm.SampleFreq);
                
                % Save xwav name paired to each detection
                GPL_struct(k).fname = files(n).name;
                
            end % For: Detection timestamp update
            
            
            % Count number of calls in previous window
            if parm.det_frequency_plot == 1
                call_freq_counter = call_freq_counter + 1;
                call_count(call_freq_counter) = length(GPL_struct);
            end
            
            
            
            % Append detections and move on to the next window
            calls = [calls, GPL_struct];
            
            % Load potential final window for output, only saved if file
            % size is reached. 
            if ~isempty(GPL_struct)
                if single_wav_flag == 0
                    end_win = xwav_struct.(field_name{n}).julian_start_time(q) + ...
                            datenum([2000,0,0,0,0,0]);
                else
                   end_win = xwav_struct.(field_name{n}).julian_start_time + ...
                            datenum([2000,0,0,0,0,0]); 
                end
            end
            
            
            % Check size of the calls and separate once greater than size
            % set by parameters
            if parm.output_size_separation == 1
                
                calls_class = whos('calls');
                
                if calls_class.bytes > parm.max_file_size*(2^30)
                    
                    file_sep_counter = file_sep_counter + 1;
                    curr_filename = sprintf('detections_fetch_NUNAT_sb_01_disk01_file%i.mat',file_sep_counter);
                      
                    % Load outputs into hyd
                    hyd(1).detection.start.time = datestr(calls(1).julian_start_time); % File start time
                    hyd(1).detection.end.window = end_win;
                    hyd(1).detection.end.time = datestr(calls(end).julian_end_time); % File end time
                    hyd(1).detection.calls = calls; % Detections
                    hyd(1).detection.parm = parm; % Parameter file
                    if parm.call_time_tracker == 1
                        hyd(1).detection.runtime = call_time_array; % Runtime per window processed
                        call_time_array = [];
                        time_counter = 0;
                    end
                    if parm.det_frequency_plot == 1
                        hyd(1).detection.call_frequency = call_count; % Number of calls per window processed
                        call_count = [];
                        call_freq_counter = 0;
                    end
                    fprintf('\nSaving Output File %i ...\n',file_sep_counter);
                    save(curr_filename,'hyd','-v7.3')
                    fprintf('\nFile %i saved successfully.\n',file_sep_counter);
                    calls = [];
                    
                end % If: File size check
                
            end % If: Output size split is on
            
            if timeBinReset_flag == 1
                parm.NumTimeBins = floor(parm.nrec/parm.fftOverlap);
            end
            
        end % Individual xwav loop

    end % Files loop

end % No padding condition






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process audio data with extra padding - Tyler Helble
% This will use the same GPL algorithm steps as set in GPL_v3, but it will
% repeat samples from the previous window in the subsequent window in
% order to account for calls that cross slate boundaries. Excluding the
% first set of samples processed, every subsequent set will repeat the last
% 2*FFTL samples from the previous set. 

% For example: For an FFTL = 2000, and 100000 samples processed at a time,
% the first loop will read samples 1-100001, the second loop will read
% samples 96001-196001, the third loop 192001-292001, etc.


if parm.pad == 1
    
    % Calculate padding values based off of samplng frequency 
    sample_padding = parm.pad*parm.SampleFreq; % # Samples to be re-processed
    num_new_samples = parm.nrec - 2*sample_padding; % # of samples that are new to each window, not including repeats from padding

    % Loop through all .x.wav files in selected folder
    for q = 1:length(files) 
  
        if length(files(q).name) > 4 
            if files(q).name(end-5:end) == '.x.wav' % Verify audio is of form .x.wav
                file_size = audioinfo(files(q).name);
                file_size = file_size.TotalSamples; % Extract # of samples in .x.wav
            end

            scale_factor = 1;
            
            % Extract the Julian start date of the current xwav from the header
            temp_name = files(q).name;
            PARAMS_pad = getxwavheaders(file_dir,temp_name); 
            julian_start_date = PARAMS_pad.ltsahd.dnumStart(1);

            % Locate deck test jumps if they occur.
            if parm.deck_test_adjustment == 1
                [j_flag,new_startTime,ideal_idx] = findTimeStepJump(PARAMS_pad,parm,file_size/num_new_samples/scale_factor,num_new_samples);
                if isempty(j_flag) % No deck test, have the user turn off param
                    return
                end
            end

            % Loop through the .x.wav file the number of times necessary
            % for all samples to be processed including the overlap. (Notes
            % this will cause the last 2*FFTL samples of the file to be
            % ignored).
            for j = 1:file_size/num_new_samples/scale_factor
    
                % Calculate the sample offset due to padding
                offset = (j-1)*num_new_samples + sample_padding;

                % Check for jump in time (deck test) and update julian
                % start date
                if parm.deck_test_adjustment == 1
                    if j == new_startTime
                        julian_start_date = PARAMS_pad.ltsahd.dnumStart(ideal_idx - 1);
                    end
                end
    
                % Print the DateTime currently being processed
                datestr(julian_start_date + datenum(2000,0,0,0,0,(offset+1-sample_padding)/parm.SampleFreq))
 
                % Verify that the sample range about to be processed does
                % not exceed file size
                if offset + 1 + sample_padding + num_new_samples <= file_size 
   
                    % Read in current range of data as determined by the
                    % offset and number of samples processed per window
                    % (parm.nrec).
                    sub_data = audioread(files(q).name,[offset + 1 - sample_padding,...
                        offset + 1 + sample_padding + num_new_samples]);

  
                    % Check for corrupted data (majority zeros)
                    error_check = sort(abs(sub_data));
                    error_check = mode(error_check);
        
                    if error_check > 0
                    
                        % Run the GPL Algorithm on the current sample set
                        [GPL_struct] = GPL_v3(sub_data,parm);
             
                        % Loop through detections found in the current
                        % window, add atsrt/end times
                        for k = 1:length(GPL_struct)

                            % Start/End times in samples (counted from
                            % start of xwav)
                            GPL_struct(k).start_time = GPL_struct(k).start_time+((j-1)*num_new_samples);
                            GPL_struct(k).end_time=GPL_struct(k).end_time+((j-1)*num_new_samples);
                            
                            % Julian start/end time for each detection. (In
                            % real time, converted from sample time)
                            GPL_struct(k).julian_start_time = julian_start_date+datenum(2000,0,0,0,0,GPL_struct(k).start_time/parm.SampleFreq);
                            GPL_struct(k).julian_end_time = julian_start_date+datenum(2000,0,0,0,0,GPL_struct(k).end_time/parm.SampleFreq);
                            
                            % Pair xwav filename to detection
                            GPL_struct(k).fname = temp_name;
                        end

                    % Append detections from current window onto previous
                    % detections
                    calls = [calls,GPL_struct];
         
                    if parm.output_size_separation == 1
                
                        calls_class = whos('calls');
                
                        if calls_class.bytes > parm.max_file_size*(2^30)
                            
                            %%%%%%%%%%% USER INPUT
                            curr_filename = sprintf('detections_v3_NUNAT_sb_01_disk01_file%i.mat',file_sep_counter);
                            %%%%%%%%%%%
                            hyd(1).detection.start = datestr(calls(1).julian_start_time);
                            hyd(1).detection.end = datestr(calls(end).julian_end_time);
                            hyd(1).detection.calls = calls;
                            hyd(1).detection.parm = parm;
                            fprintf('\nSaving Output File %i ...\n',file_sep_counter);
                            save(curr_filename,'hyd','-v7.3')
                            fprintf('\nFile %i saved successfully.\n',file_sep_counter);
                            calls = [];
                    
                        end % If: File size check
                
                    end % If: Output size split is on
                    
                    
                    else % If error check fails, write report

                    fileID = fopen('errors.txt','at');
                    fprintf(fileID, files(q).name,'\n');
                    fprintf(fileID, datestr(julian_start_date+datenum(2000,0,0,0,0,(offset+1-sample_padding)/parm.sample_freq)),'/n');
                    fclose(fileID);
            
                    end % End processing of the current window
                    
                end % End wav length check
            end % End loop through current xwav

        end % End filename size check
    end % End loop through all xwav files in selected folder

end % End padding version check




%% Save detections into output structure 

% Split outputs into smaller files - Ian Cosgrove
if parm.output_size_separation == 1 
    
    if ~isempty(calls) % Save detections leftover from file split-up 
        
        %%%%%%%%%%%%% USER INPUT
        last_filename = sprintf('detections_fetch_NUNAT_sb_01_disk01_file%i.mat',file_sep_counter+1);
        %%%%%%%%%%%%%
        hyd(1).detection.start.time = datestr(calls(1).julian_start_time);
        hyd(1).detection.end.window = q;
        hyd(1).detection.end.time = datestr(calls(end).julian_end_time);
        hyd(1).detection.calls = calls;
        hyd(1).detection.parm = parm;
        fprintf('\nSaving Output File %i ...\n',file_sep_counter + 1);
        save(last_filename,'hyd','-v7.3')
        fprintf('\nFile %i saved successfully.\n',file_sep_counter + 1);
        
    end
    
else % Normal saving with no split up - Tyler Helble
    
    % Load detection and parameter data into 'hyd'
    hyd(1).detection.calls = calls;
    hyd(1).detection.parm = parm;

    % Set detection name
    detection_name = (strcat('v2_D_parm_test'));

    % Save final output
    save(strcat(detection_name,files(1).name(1:end-6),'.mat'),'hyd')
    
end


toc(total_elapsed_time) % Print runtime
