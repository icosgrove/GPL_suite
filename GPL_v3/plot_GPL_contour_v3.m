% Plot call contours from GPL detector output

% This script plots 'cm', 'cm_max', and 'cm_max2' from the GPL detector
% output .mat file. 'Cm' is all the islands of the call (what the detector
% thinks is actually the part of the spectrogram that is the call). Cm_max
% is the single strongest contour of the call (island of the call with the most
% energy). Cm_max2 is the single second strongest island of the call.

% In the section of this script named 'Load', add your detection file and
% which detection # to plot. 

% If you plot with seconds/Hz, black dashed lines will be plotted showing
% your frequency range that was processed from parameters. If you plot
% using bins, the frequency bin range will be plotted as a red dashed line.

% Setting waveform to 0/1 will not plot/plot the waveform of cm_max, the
% strongest contour.

% Written by Ian Cosgrove Sept. 2023, updated to v3 Apr. 2024


%% Close previous figure
close(figure(8))

%% Load Detection file

load v2_D_parm_test90_30_primary_NUNAT_SB01_20210816_060302.wav.mat

call = 1; % detection number to be plotted

% Switches:  1: on, 0: off
waveform = 0;   % plot cm_max waveform
display_measurements = 1;   % show contour measurements






%% Recreate Cm, Cm_max, Cm_max2
% Load cm: full contour characteristics
cm_vals = single(hyd.detection.calls(call).cm.values);
cm_index = hyd.detection.calls(call).cm.index;
peak_energy1 = hyd.detection.calls(call).cm.scale;
contour_size = hyd.detection.calls(call).cm.size;

spec = zeros(contour_size(1), contour_size(2)); % allocate matrix the size of th e
scaled_vals1 = cm_vals*peak_energy1/(2^16-1); % scale values from 0-2^16-1 back to 0-peak_energy

cm = spec;
cm(cm_index) = scaled_vals1; % write contour values to there respective indices in the spectrogram

% Load cm_max: strongest contour characteristics
max_vals = single(hyd.detection.calls(call).cm_max.values);
max_index = hyd.detection.calls(call).cm_max.index;
peak_energy2 = hyd.detection.calls(call).cm_max.scale;

scaled_vals2 = max_vals*peak_energy2/(2^16-1); % scale values from 0-2^16-1 back to 0-peak_energy

cm_max = spec;
cm_max(max_index) = scaled_vals2;


% Load cm_max2: second strongest contour characteristics
max2_vals = single(hyd.detection.calls(call).cm_max2.values);
max2_index = hyd.detection.calls(call).cm_max2.index;
peak_energy3 = hyd.detection.calls(call).cm_max2.scale;

scaled_vals3 = max2_vals*peak_energy3/(2^16-1); % scale values from 0-2^16-1 back to 0-peak_energy

cm_max2 = spec;
cm_max2(max2_index) = scaled_vals3;



call_title = sprintf('Detection #%i', call);
s_title = {call_title,' '};

% cm_max Slope 
if isfield(hyd.detection.calls,'cm_max')
    cm_max_slope = hyd.detection.calls(call).cm_max.slope;
    cm_max_slope_y_int = hyd.detection.calls(call).cm_max.slope_y_int;
    cm_max_slope_data = (1:contour_size(2))*cm_max_slope + cm_max_slope_y_int + hyd.detection.parm.FreqBinLo;
end

% cm_max2 slope
if isfield(hyd.detection.calls,'cm_max2')
    cm_max2_slope = hyd.detection.calls(call).cm_max2.slope;
    cm_max2_slope_y_int = hyd.detection.calls(call).cm_max2.slope_y_int;
    cm_max2_slope_data = (1:contour_size(2))*cm_max2_slope + cm_max2_slope_y_int + hyd.detection.parm.FreqBinLo;
end

% cm slope
if isfield(hyd.detection.calls,'cm')
    cm_slope = hyd.detection.calls(call).cm.slope;
    cm_slope_y_int = hyd.detection.calls(call).cm.slope_y_int;
    cm_slope_data = (1:contour_size(2))*cm_slope + cm_slope_y_int + hyd.detection.parm.FreqBinLo;
end
    
%% Plot





if isfield(hyd.detection.calls,'cm')

    freq=(0:hyd.detection.parm.fftl/2) /hyd.detection.parm.fftl *hyd.detection.parm.SampleFreq; % 0-1000Hz, 1025 indices
    freq=freq(hyd.detection.parm.FreqBinLo:hyd.detection.parm.FreqBinHi); % restrict to processed frequency bins 
    [x,y] = size(cm); % contour size

    figure(8); hold on
    suptitle(s_title)
    subplot(1,3,1); hold on % cm
    box on; grid on
    time = (1:y)*(hyd.detection.parm.fftOverlap/hyd.detection.parm.SampleFreq);
    surf(time,freq,cm); view(2)
    
    % Plot Slope
    slope_y = time*cm_slope + cm_slope_y_int;
    slope_z = max(max(cm))*ones(1,length(slope_y));
    plot3(time,slope_y,slope_z,'--r','LineWidth',2);
    
%     freq_x = (hyd.detection.parm.fftOverlap/hyd.detection.parm.SampleFreq)*(1:contour_size(2)); 
%     min_freq = hyd.detection.parm.freq_lo*ones(length(freq_x));
%     max_freq = hyd.detection.parm.freq_hi*ones(length(min_freq));
%     plot(freq_x,min_freq,'k--','LineWidth',2); % plot dashed line at minimum frequency
%     plot(freq_x,max_freq,'k--','LineWidth',2); % plot dashed line at maximum frequency
    axis([time(1) time(end) hyd.detection.parm.freq_lo ...
        hyd.detection.parm.freq_hi]) % set axis to time range and the frequency range +/- 5 Hz
    xlabel('Time [s]'); ylabel('Frequency [Hz]'); title('Cm: All Contours');
    hold off
end

if isfield(hyd.detection.calls,'cm_max')
    subplot(1,3,2); hold on % cm_max
    box on; grid on
    surf(time,freq,cm_max); view(2)
    
    slope_y = time*cm_max_slope + cm_max_slope_y_int;
    slope_z = max(max(cm_max))*ones(1,length(slope_y));
    plot3(time,slope_y,slope_z,'--r','LineWidth',2);
    
%     plot(freq_x,min_freq,'k--','LineWidth',2); plot(min_freq,max_freq,'k--','LineWidth',2);
     axis([time(1) time(end) hyd.detection.parm.freq_lo ...
        hyd.detection.parm.freq_hi])
    xlabel('Time [s]'); ylabel('Frequency [Hz]'); title('Cm Max: Strongest Contour');
    hold off
end

if isfield(hyd.detection.calls,'cm_max2')
    subplot(1,3,3); hold on % cm_max2
    box on; grid on
    surf(time,freq,cm_max2); view(2)
    
    slope_y = time*cm_max2_slope + cm_max2_slope_y_int;
    slope_z = max(max(cm_max2))*ones(1,length(slope_y));
    plot3(time,slope_y,slope_z,'--r','LineWidth',2);
    
%     plot(freq_x,min_freq,'k--','LineWidth',2); plot(min_freq,max_freq,'k--','LineWidth',2);
     axis([time(1) time(end) hyd.detection.parm.freq_lo ...
        hyd.detection.parm.freq_hi])
    xlabel('Time [s]'); ylabel('Frequency [Hz]');
    title('Cm Max2: Second Strongest Contour');
    hold off
    hold off

end
    

            
if waveform == 1
    
    figure(11)
    plot(hyd.detection.calls(call).cm_max_waveform);
    t_string = sprintf('Detection #%i Waveform of Cm max', call);
    xlabel('Samples'); ylabel('Amplitude'); title(t_string);
    axis([0 length(hyd.detection.calls(call).cm_max_waveform)...
        0.1*min(hyd.detection.calls(call).cm_max_waveform)+...
        min(hyd.detection.calls(call).cm_max_waveform) ...
        0.1*max(hyd.detection.calls(call).cm_max_waveform)+...
        max(hyd.detection.calls(call).cm_max_waveform)])
    box on; grid on
            
end
            
            

if display_measurements == 1
    
      fprintf('\nDetection #%i\n', call);
      fprintf('Start Time: %s\n',datestr(hyd.detection.calls(call).julian_start_time));
      fprintf('End Time: %s\n\n',datestr(hyd.detection.calls(call).julian_end_time));
      
      if isfield(hyd.detection.calls,'cm')
        fprintf('Cm:\n')
        fprintf('\tDuration [s]: %f\n\tSlope [Hz/s]: %f\n\tFrequency Bandwidth [Hz]: %f\n\tStart Freq. [Hz]: %f\n\tEnd freq. [Hz]: %f\n\tMin. Freq. [Hz]: %f\n\tMax. Freq. [Hz]: %f\n',...
            hyd.detection.calls(call).cm.duration_sec, hyd.detection.calls(call).cm.slope, ...
            hyd.detection.calls(call).cm.freq_bandwidth_hz, hyd.detection.calls(call).cm.start_freq_hz, ...
            hyd.detection.calls(call).cm.end_freq_hz, hyd.detection.calls(call).cm.min_freq_hz, ...
            hyd.detection.calls(call).cm.max_freq_hz);
      end
      
      if isfield(hyd.detection.calls,'cm_max')
          fprintf('\nCm_max:\n')
          fprintf('\tDuration [s]: %f\n\tSlope [Hz/s]: %f\n\tFrequency Bandwidth [Hz]: %f\n\tStart Freq. [Hz]: %f\n\tEnd freq. [Hz]: %f\n\tMin. Freq. [Hz]: %f\n\tMax. Freq. [Hz]: %f\n',...
              hyd.detection.calls(call).cm_max.duration_sec, hyd.detection.calls(call).cm_max.slope, ...
              hyd.detection.calls(call).cm_max.freq_bandwidth_hz, hyd.detection.calls(call).cm_max.start_freq_hz, ...
              hyd.detection.calls(call).cm_max.end_freq_hz, hyd.detection.calls(call).cm_max.min_freq_hz, ...
              hyd.detection.calls(call).cm_max.max_freq_hz);
      end
      
      if isfield(hyd.detection.calls,'cm_max2')
          fprintf('\nCm_max2:\n')
          fprintf('\tDuration [s]: %f\n\tSlope [Hz/s]: %f\n\tFrequency Bandwidth [Hz]: %f\n\tStart Freq. [Hz]: %f\n\tEnd freq. [Hz]: %f\n\tMin. Freq. [Hz]: %f\n\tMax. Freq. [Hz]: %f\n',...
              hyd.detection.calls(call).cm_max2.duration_sec, hyd.detection.calls(call).cm_max2.slope, ...
              hyd.detection.calls(call).cm_max2.freq_bandwidth_hz, hyd.detection.calls(call).cm_max2.start_freq_hz, ...
              hyd.detection.calls(call).cm_max2.end_freq_hz, hyd.detection.calls(call).cm_max2.min_freq_hz, ...
              hyd.detection.calls(call).cm_max2.max_freq_hz);
      end
      
end   
            
            