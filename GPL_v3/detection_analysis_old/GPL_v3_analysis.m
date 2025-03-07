%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script will perform detection performance analysis on GPL_v3
% detections.
%
% Functions:
%  - Determine if GPL_v3 successfully recovered a manual detection.
%    Criteria for a manual detection detection to be successfully detected is
%    that it must have a GPL detection start time within a margin of time
%    around the real start time OR there is a GPL detection start time that
%    is bounded by a manual start/end time.
% 
%  - Count how many Boundary cross cases appear
%
%  - Plot Detection Frequency and associated runtime
%
% Written by Ian Cosgrove 04/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load in file directory containing .mat detection files
%%% Tyler Helble's code

fprintf('\nSelect folder containing detection .mat files to process\n');
[filename, pathname] = uigetfile('*.mat'); % Select File
cwd = pwd; 
cd(pathname) % Set current directory to path containing xwav file
file_dir = pwd; 
addpath(pwd); 
files = dir('*.mat'); % Variable files contains xwav data
cd(cwd); % Set current directory back to current working directory

%%% End Tyler's code

%% Sort .mat files by file number 

order = 1;
k = 0;
if length(files) > 1
    
    % Determine the order of the files by name. Note that the number of the
    % file must be the final character before '.mat'. Example:
    % det_sb_01_file1.mat, det_sb_01_file2.mat, etc...
    for n = 1:length(files)
        fname = files(n).name;
        while ~isnan(str2double(string(fname(end-4-k:end-4))))
            order(n) = str2double(string(fname(end-4-k:end-4)));
            k = k + 1;
        end
        k = 0;
    end
end

%% Load in manual detlog 

manual_data = readtable('CORRECTED_NUNAT_SB_01_disk01_90_30_Hz_down.xls');
manual_start_times = manual_data.StartTime; % Julian start time
julian_mst = datenum(manual_start_times);

%% Select what analysis you want performed

boundary_cross_count = 1; % Count how many calls crossed boundary windows
precision_recall = 1; % Find the precision and recall based on provided manual detections
call_frequency = 1;
call_runtime = 1;


%% Analysis

% Allocation
bc_count = 0; % Boundary counter
TP = 0; % True positive counter
FN = 0; % False negative counter
FP = 0; % False positive counter
a1 = 1; % Counter
a2 = 1; % Counter
gpl_tp = []; % TP array
gpl_start_times_counter = 0; % GPL counter
num_addon = 0;

% Configuration
amount_of_time = 1; % Margin of time +/- a manual detection that a GPL detection must fall within
time_units = 'S'; % Units for that time
warning('off')

% Loop over all .mat files and perform requested analysis
for n = 11:length(files)
    
    corr_fname = files(order == n);
    fname = sprintf('%s',corr_fname.name);
    fprintf('\nLoading File %i ...\n',n); tic
    load(fname)
    fprintf('\nFile successfully loaded\n'); toc
    fprintf('\nAnalyzing File %i\n',n); tic
    
    % Count number of calls that cross window boundaries
    if boundary_cross_count == 1 
        
        crossInd = find([hyd.detection.calls.boundary_cross_flag] == 1); % ID crossing cases
        bc_count = bc_count + length(crossInd); % Add to running total
        
    end
    
    % Extract manual start times within the region of current file
    mst_a = find(julian_mst > datenum(hyd.detection.start.time));
    mst_b = find(julian_mst < datenum(hyd.detection.end.time));
    mst_ind = intersect(mst_a,mst_b);
    current_mst = manual_start_times(mst_ind);
    
    
    % Compute Precision/Recall based of manual detection set
    if precision_recall == 1
        
        gpl_start_times = [hyd.detection.calls.julian_start_time]';
        gpl_end_times = [hyd.detection.calls.julian_end_time]';
        gpl_start_times_counter = gpl_start_times_counter + length(gpl_start_times);
        
        for k = 1:length(current_mst)
            
            if k > 1
                num_addon = length(current_mst);
            end
        
            % Find the margin of julian time around current manual det
            [~, datep1, datem1] = julian_time_conversions(current_mst(k), amount_of_time, time_units);
            datep1 = datenum(datep1);
            datem1 = datenum(datem1);

            % Identify GPL detections that are a potential match
            x1 = find(gpl_start_times > datem1); % Locate dets above lower bound
            x2 = find(gpl_start_times < datep1); % Locate dets underneath upper bound
            potential_hit = intersect(x1,x2); % Find overlap between calls above/below (potential match)

            % If no hit within margin, do a second check to see if manual
            % start time is bounded by a GPL start/end time.
            if isempty(potential_hit) && ~isempty(x2)
                closest_det_below = max(x2); 
                cd_s = gpl_start_times(closest_det_below);
                cd_e = gpl_end_times(closest_det_below);
                det = datenum(current_mst(k));
                if (det > cd_s) && (det < cd_e) % Detection is bounded on both sides by the start/end of a GPL detection
                    potential_hit = max(x2);
                end
            end

            if ~isempty(potential_hit) % TP

                TP = TP + 1;
                
                tp_array{a1,1} = k + num_addon; % write the manual det # 
                tp_array{a1,2} = datestr(current_mst(k)); % Write manual start time in datetime
                
                if length(potential_hit) > 1 % Multiple GPL Dets within margin

                    q = length(potential_hit);

                    for s = 1:q % Loop over number of multiples

                        % Find absolute difference in time 
                        difference(s) = abs(datenum(current_mst(k)) - gpl_start_times(potential_hit(s))); % Find each time difference

                    end

                    % Consider the detection with the smallest time
                    % difference to be the TP
                    closest_det = find(difference == min(difference)); % Find closest GPL detection and save that
                    tp_array{a1,3} = datestr(gpl_start_times(potential_hit(closest_det)));
                    tp_array{a1,4} = potential_hit(closest_det); % GPL det # for hit

                else % Only one GPL detection inside margin

                    tp_array{a1,3} = datestr(gpl_start_times(potential_hit)); % Only 1 GPL hit
                    tp_array{a1,4} = potential_hit;

                end

                a1 = a1 + 1;

            else % FN because no detection in region

                FN = FN + 1;
                fn_array{a2,1} = k + num_addon; % manual det #
                fn_array{a2,2} = datestr(current_mst(k)); % manual det time
                closest_det_below = gpl_start_times(max(x2)); % Locate nearest GPL detection before manual
                fn_array{a2,3} = datestr(closest_det_below);
                fn_array{a2,4} = max(x2);
                closest_det_above = gpl_start_times(min(x1)); % Locate nearest GPL detections after manual
                fn_array{a2,5} = datestr(closest_det_above);
                fn_array{a2,6} = min(x1);
                fn_array{a2,7} = current_mst(k);

                a2 = a2 + 1;

            end

        end % End manual start time comparison

        % Calculate # of FP
        TP_out(n) = TP;
        FP_out(n) = length(gpl_start_times) - TP_out(n);
        FN_out(n) = FN;
        TP = 0; FP = 0; FN = 0;
        
    end % Precision/Recall is on

    fprintf('\nFile %i analysis complete.\n',n); toc
    
    % Create Call Frequency and Runtime plots
    if call_frequency == 1
        
        % Create vector of timestamps for x-axis
        ts = datenum(hyd.detection.start.time);
        tf = datenum(hyd.detection.end.time);
        space = datenum([2000 0 0 0 0 hyd.detection.parm.nrec/hyd.detection.parm.SampleFreq]) - ...
            datenum([2000 0 0 0 0 0]);
        time_vec = ts:space:tf;
        
        figure(10); hold on
        
        % Plot Manual detection locations
        max_det_count = max(hyd.detection.call_frequency)+5;
        det_rep = repmat(datenum(current_mst),1,2)';
        count_rep = repmat([0 max_det_count],length(current_mst),1)';
        plot(det_rep,count_rep,'-r','LineWidth',1)
        legend('Detection Count','Manual Detections')
        if length(time_vec) ~= length(hyd.detection.call_frequency)
            if length(hyd.detection.call_frequency) == length(time_vec)+1
                time_vec = [time_vec time_vec(end)+datenum([0 0 0 0 0 hyd.detection.parm.nrec/hyd.detection.parm.SampleFreq])];
            else
                disp('Size Error between number of windows expected and received in freq. counter')
                return
            end
        end
        plot(time_vec,hyd.detection.call_frequency,'-b','LineWidth',1)
        xlabel('Date'); ylabel('Number of Detections')
        title('Detection Count per Window')
        box on; grid on
        axis([min(time_vec) max(time_vec) 0 max(hyd.detection.call_frequency)+5])
        
        % Create xticks dynamically with figure display size
        fig1 = gcf;
        fig1_size = fig1.Position(3:4);
        label_count = 8; % Number of x-axis tick markers
        label_calc = min(label_count, round(fig1_size(1)/100));
        xticks(linspace(min(time_vec), max(time_vec), label_calc))
        date_vec = datestr(linspace(min(time_vec), max(time_vec), label_calc), 0);
        xticklabels(date_vec)
        xtickangle(20)
        
        pause
        clear cfg
        
    end % If: Call Frequency Plot
    
end % For: loop over all .mat files



%% Results

if precision_recall == 1
    
    % Verification array formatting
    addon = {1 1 1 1 };
    tp_array = [addon; tp_array];
    tp_array{1,1} = 'Manual Det #';
    tp_array{1,2} = 'Manual Start Time';
    tp_array{1,3} = 'GPL Start Time';
    tp_array{1,4} = 'GPL Det #';
    addon = {1 1 1 1 1 1 1};
    fn_array = [addon; fn_array];
    fn_array{1,1} = 'Manual Det #';
    fn_array{1,2} = 'Manual Start Time';
    fn_array{1,3} = 'Previous GPL Det';
    fn_array{1,4} = 'Previous GPL #';
    fn_array{1,5} = 'Next GPL Det';
    fn_array{1,6} = 'Next GPL Det #';
    fn_array{1,7} = 'Julian Time';

    % Missed detection times output
    fn_times = fn_array(2:end,:);
    fn_times(:,4) = []; % Clear Previous GPL det#
    for m = 2:size(fn_array,1)
        fn_j_times{m-1,1} = fn_array{m,end};
        fn_j_times{m-1,2} = fn_array{m,3};
        fn_j_times{m-1,3} = fn_array{m,5};
    end
    %save('gpl_fn.mat','fn_j_times');

    % Calculate Precision and Recall:
    TP = sum(TP_out); FP = sum(FP_out); FN = sum(FN_out);
    P = TP / (TP + FP);
    R = TP / (TP + FN);

    % Print Results:
    fprintf('\nMargin: %i%s\nTotal GPL Detections: %i\nTotal Manual Detections: %i\nTrue Positives: %i\nFalse Negatives: %i\nFalse Positives: %i\nPrecision: %f\nRecall: %f\n',...
        amount_of_time,time_units,gpl_start_times_counter, length(manual_start_times), TP, FN, FP, P, R);
    
end

if boundary_cross_count == 1
    fprintf('\nNumber of calls crossing window boundaries: %i\n',bc_count);
end



%% Plot the detection frequency and associated runtime vs. deployment time















