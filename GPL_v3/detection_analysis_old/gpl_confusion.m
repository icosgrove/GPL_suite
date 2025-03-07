% Compare GPL detections to manual timeseries 
% Written by Ian Cosgrove Nov. 2023

clear fn_array

%% Load Mannual/GPL Data

if ~exist('hyd','var')
    load('detections_NUNAT_SB_01_disk01_r10.mat'); % GPL data
end

gpl_start_times = [hyd.detection.calls.julian_start_time]';
gpl_end_times = [hyd.detection.calls.julian_end_time]';

manual_data = readtable('CORRECTED_NUNAT_SB_01_disk01_90_30_Hz_down.xls');

manual_start_times = manual_data.StartTime; % Julian start time



%% Comparison Loop

% Compare results within 2 seconds of each other
amount_of_time = 1;
time_units = 'S';

% Pre-allocate
TP = 0;
FN = 0;
m = 1;
gpl_comp = zeros(length(gpl_start_times),2);
missed_det = zeros(length(manual_start_times),2);

a1 = 1;
a2 = 1;
gpl_tp = [];

for h = 1:length(manual_start_times) % Run by processing each manual start time and checking for GPL dets
    
    % Find the margin of julian time around current manual det
    [~, datep1, datem1] = julian_time_conversions(manual_start_times(h), amount_of_time, time_units);

    datep1 = datenum(datep1);
    datem1 = datenum(datem1);
    
    x1 = find(gpl_start_times > datem1); % Locate dets above lower bound
    x2 = find(gpl_start_times < datep1); % Locate dets underneath upper bound
    x3 = intersect(x1,x2); % Find overlap between calls above/below (potential match)

    if isempty(x3) % Check if manual lies within gpl, not necessarily within start time region
        closest_det_below = max(x2);
        cd_s = gpl_start_times(closest_det_below);
        cd_e = gpl_end_times(closest_det_below);
        det = datenum(manual_start_times(h));
        if (det > cd_s) && (det < cd_e) % Detection is bounded on both sides by the start/end of a GPL detection
            x3 = max(x2);
        end
    end
    
    if ~isempty(x3) % TP
        
        TP = TP + 1;
        tp_array{a1,1} = h; % write the manual det # 
        tp_array{a1,2} = datestr(manual_start_times(h)); % Write manual start time in datetime
        
        if length(x3) > 1 % Multiple GPL Dets within margin
            
            q = length(x3);
            
            for s = 1:q % Loop over number of multiples
                
                difference(s) = abs(datenum(manual_start_times(h)) - gpl_start_times(x3(s))); % Find each time difference
                
            end
            
            closest_det = find(difference == min(difference)); % Find closest GPL detection and save that
            tp_array{a1,3} = datestr(gpl_start_times(x3(closest_det)));
            tp_array{a1,4} = x3(closest_det); % GPL det # for hit
            gpl_tp = [gpl_tp gpl_start_times(x3)']; % record gpl hits for FP array
            
        else
            
            tp_array{a1,3} = datestr(gpl_start_times(x3)); % Only 1 GPL hit
            tp_array{a1,4} = x3;
            gpl_tp = [gpl_tp gpl_start_times(x3)]; % Record gpl hits for FP array
            
        end
        
        a1 = a1 + 1;
        
    else % x3 is empty, so FN
        
        FN = FN + 1;
        fn_array{a2,1} = h; % manual det #
        fn_array{a2,2} = datestr(manual_start_times(h)); % manula det time
        
        closest_det_below = gpl_start_times(max(x2)); % Locate nearest GPL detection before manual
        fn_array{a2,3} = datestr(closest_det_below);
        fn_array{a2,4} = max(x2);
        closest_det_above = gpl_start_times(min(x1)); % Locate nearest GPL detections after manual
        fn_array{a2,5} = datestr(closest_det_above);
        fn_array{a2,6} = min(x1);
        fn_array{a2,7} = manual_start_times(h);

        a2 = a2 + 1;
            
    end
    
end % End manual start time comparison

% Calculate # of FP
FP = length(gpl_start_times) - TP;


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

clear a b cm
%% Uncomment to see the sp of the TP 
% for r = 2:length(tp_array)
%     
%         call = cell2mat(tp_array(r,4));
%         cm_vals = single(hyd.detection.calls(call).cm.values);
%         cm_index = hyd.detection.calls(call).cm.index;
%         peak_energy1 = hyd.detection.calls(call).cm.scale;
%         contour_size = hyd.detection.calls(call).cm.size;
% 
%         spec = zeros(contour_size(1), contour_size(2)); % allocate matrix the size of th e
%         scaled_vals1 = cm_vals*peak_energy1/(2^16-1); % scale values from 0-2^16-1 back to 0-peak_energy
% 
%         cm = spec;
%         cm(cm_index) = scaled_vals1;
%     
%         cm_slope = hyd.detection.calls(call).cm.slope;
%         cm_slope_y_int = hyd.detection.calls(call).cm.slope_y_int;
%         cm_slope_data = (1:contour_size(2))*cm_slope + cm_slope_y_int + hyd.detection.parm.bin_lo;
% 
%     
%         y_start = hyd.detection.parm.bin_lo; y_end = hyd.detection.parm.bin_hi; % frequency bin range 
%         x_start = 1; x_end = contour_size(2); % time bin range 
%         
%         x = x_start:x_end;
%         y = y_start:y_end;
%         
%         freq=(0:hyd.detection.parm.fftl/2) /hyd.detection.parm.fftl *hyd.detection.parm.sample_freq; % 0-1000Hz, 1025 indices
%         freq=freq(hyd.detection.parm.bin_lo:hyd.detection.parm.bin_hi); % restrict to processed frequency bins 
%         
%         time = (hyd.detection.parm.skip/hyd.detection.parm.sample_freq) * x;
%         
%         
%         
%         f = figure(9); hold on
%         box on; grid on
%         for n = 1:contour_size(2)
%             for m = 1:length(y)
%                 a(n) = n*(hyd.detection.parm.skip/hyd.detection.parm.sample_freq); % convert time bins to seconds
%                 b(m) = freq(m); % convert frequency bins to Hz
%             end
%         end
%         surf(a,b,cm); shading faceted
%         c = (hyd.detection.parm.skip/hyd.detection.parm.sample_freq)*(1:contour_size(2)); 
%         d = hyd.detection.parm.freq_lo*ones(length(c));
%         e = hyd.detection.parm.freq_hi*ones(length(c));
%         plot(c,d,'r--','LineWidth',2); % plot dashed line at minimum frequency
%         plot(c,e,'r--','LineWidth',2); % plot dashed line at maximum frequency
%         axis([a(1) a(end) b(1) b(end)]);
%         xlabel('Time (s)'); ylabel('Frequency (Hz)'); 
%         title_string = sprintf('Manual Det#%i GPL Det#%i ', tp_array{r,1}, tp_array{r,4});
%         title(title_string);
%         hold off
%         
% 
%         pause;
%         clf
%         clear a b cm
% end
%%







fn_times = fn_array(2:end,:);
fn_times(:,4) = []; % Clear Previous GPL det#
for m = 2:length(fn_array)
    fn_j_times{m-1,1} = fn_array{m,end};
    fn_j_times{m-1,2} = fn_array{m,3};
    fn_j_times{m-1,3} = fn_array{m,5};
end
save('gpl_fn.mat','fn_j_times');






  

% Calculate Precision and Recall:
P = TP / (TP + FP);
R = TP / (TP + FN);


% Print Results:
fprintf('\nMargin: %i%s\nTotal GPL Detections: %i\nTotal Manual Detections: %i\nTrue Positives: %i\nFalse Negatives: %i\nFalse Positives: %i\nPrecision: %f\nRecall: %f\n',amount_of_time,time_units,length(gpl_start_times), length(manual_start_times), TP, FN, FP, P, R);

    
    
    
    
    
    















