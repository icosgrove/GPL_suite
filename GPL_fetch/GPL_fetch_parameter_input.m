function [parm] = GPL_fetch_parameter_input(PARAMS,xwav_struct,parm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPL_parameter_input_v3 is a script to set fixed input parameters for the
% GPL detector algorithm. Some parameters are deterined by xwav header 
% information, others are manual input Some of these parameters will change 
% during the GPL process. 

% Written by Ian Cosgrove 01/26/2024
% Lsst Updated 07/23/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Multiplier for Tyler's spectrogram padding. Set to ~=0 to use
% Tyler's original spectrogram padding process. The value of the variable
% will determine how many seconds are padded around each spectrogram.
parm.pad = 0;

% Sample Frequency [samples/s]: Value straight from xwav header
parm.SampleFreq = PARAMS.ltsahd.sample_rate(1);

% Number of seconds to be processed at one time [seconds]: Ensure that your
% call type will not be longer in duration than the window size.
parm.window_length = 75;

% Number of seconds per disk write regardless of selected window size
parm.disk_write_length = round((xwav_struct.xwav1.julian_start_time(end) - xwav_struct.xwav1.julian_start_time(end-1))*24*60*60);
if parm.disk_write_length ~= 75
    fprintf('\nWarning: Disk write length not 75s, check for errors in data and adjust parm.disk_write_length\n')
    pause
end


% Number of samples per wav 
if xwav_struct.xwav1.nwav == 1 % xwav downloaded directly from Triton, only one wav
    parm.nrec = parm.window_length*parm.SampleFreq; 
else % xwav data in general, already split up every 75s by disk writes
    parm.nrec = (xwav_struct.xwav1.TotalSamples/xwav_struct.xwav1.nwav)*(parm.window_length/parm.disk_write_length);
end


% FIFO Removal switch: Remove n*xHz bins and replace them with average of
% the adacent bins. Done to remove 50Hz (or xHz) harmonics.
parm.fifo_removal = 0; % 1:on, 0:off
parm.fund_fifo_freq = 50; % Fundamental frequency of FIFO noise, rest are harmonics (integer multiples)
parm.fifo3bin = 1; % Average and replace the n*(49:51) frequency bins 


% Helble et al 2012 eq 6 (nu1 and nu2)
% GPL Algorithm Parameters
parm.xp1 = 2.5;     % Exponent for u: Normed over frequency
parm.xp2 = 0.5;     % Exponent for y: Normed over time
parm.whiten = 1;    % 1: Whiten in time only. 2: Whiten in time and frequency
parm.white_x = 1.5; % Enhanced mean noise multiplier
parm.NumLoops = 5;  % max number of loops each window is processed for calls

% Frequency range for GPL to process for calls [Hz]
parm.freq_lo = 15;
parm.freq_hi = 30;

% Frequency range to display for context in GPL Fetch [Hz]. (Set to be the
% same as parm.freq_lo/hi if you want no extra context). 
parm.fp1.freqDelimitLo = 10;
parm.fp1.freqDelimitHi = 35;

% Input warnings
if (parm.fp1.freqDelimitLo > parm.freq_lo) || (parm.fp1.freqDelimitHi < parm.freq_hi)
    disp('Warning: spectrogram context size misaligned with GPL window')
    pause
end
if parm.freq_lo > parm.freq_hi
    fprintf('\nWarning: Frequency limits possibly mismatched\n')
    pause
end

% Call duration bounds [s]
parm.MinCallDuration = 0.1;
parm.MaxCallDuration = 10;

% Spectral Parameters
parm.fftl = 2000; % FFT Length
parm.fftOverlapPercent = 90; % FFT overlap percentage [Format: xx% NOT 0.xx]
if parm.fftOverlapPercent == 0
    parm.fftOverlap = round(parm.fftl);
else
    parm.fftOverlap = round(parm.fftl - (parm.fftOverlapPercent*1e-2)*parm.fftl); % Calculate number of samples to step forward based on overlap %
end

% Time bin parameters:
parm.NumTimeBins = floor(parm.nrec/parm.fftOverlap);
    
% Frequency bin parameters:
parm.freq_range = 0:parm.SampleFreq/parm.fftl:(parm.fftl/2)/parm.fftl*parm.SampleFreq; % Define frequency range from FFTL

% GPL frequency bins to process
[~,parm.FreqBinLo] = min(abs(parm.freq_range - parm.freq_lo));
[~,parm.FreqBinHi] = min(abs(parm.freq_range - parm.freq_hi));
parm.NumFreqBins = parm.FreqBinHi - parm.FreqBinLo + 1;

% Context frequency bins
[~,parm.fp1.ContextFreqBinLo] = min(abs(parm.freq_range - parm.fp1.freqDelimitLo));
[~,parm.fp1.ContextFreqBinHi] = min(abs(parm.freq_range - parm.fp1.freqDelimitHi));
parm.ContextNumFreqBins = parm.fp1.ContextFreqBinHi - parm.fp1.ContextFreqBinLo + 1;

% Frequency bin range for summation
parm.SumFreqBinLo = parm.FreqBinLo;
parm.SumFreqBinHi = parm.FreqBinHi;

% 1-D
% Detection/Contour Threshold Parameters
parm.noise_ceiling = 20; % Upper limit for backgorund noise, multiple of noise floor
parm.thresh = 90; % Threshold multiplier of noise floor for an energy spike to be considered a detection

if parm.noise_ceiling > parm.thresh % Check threshold
    disp('Warning: Noise Ceiling is larger than threshold')
    pause
end

% 2-D
parm.ContourCutoff = 6; % Energy threshold for time,frequency coordinate to be included to the call contour
parm.MinIslandSize = 10; % Minimum number of pixels above cutoff for a detection to be considered a contour
parm.smaller_island_cutoff = 0.05; % All islands smaller than the largest island must be at least this percentage size of the largest island.


%% Measurement On/Off Switches
parm.waveform_on = 0; % Extract the original time signal waveform of the call
parm.template_on = 1; % Extract call contour and measurements (0: only start/end time is recorded)
parm.cm_on = 1; % Extract the call contour 
parm.cm_max_on = 1; % Extract the strongesr contour of the call
parm.cm_max2_on = 1; % Extract the second strongest contour of the call
parm.measurements_on = 1; % Obtain measurements of the contours
% parm.filter_on = 1; % Run primary filters 


    % Measurement on/off swtiches
    parm.measure.cm_freq_limits = 1; % Measure frequency start/end for Cm
    parm.measure.cm_max_freq_limits = 1; % " " for cm_max
    parm.measure.cm_max2_freq_limits = 1; % " " for cm_max2

    parm.measure.cm_freq_bandwidth = 1; % Measure freq bandwidth for cm
    parm.measure.cm_max_freq_bandwidth = 1; % " " for cm_max
    parm.measure.cm_max2_freq_bandwidth = 1; % " ' for cm_max2

    parm.measure.cm_duration = 1; % Measure call duration for cm
    parm.measure.cm_max_duration = 1; % " " for cm_max
    parm.measure.cm_max2_duration = 1; % " " for cm_max2

    parm.measure.cm_slope = 1; % Measure slope for cm
    parm.measure.cm_max_slope = 1; % " " for cm_max
    parm.measure.cm_max2_slope = 1; % " " for cm_max2

    parm.measure.spec_noise = 1; % Measure background noise of the call's xwav
    parm.measure.RL_min_size = 5; % Minimum contour size to allow RL to be found
    parm.measure.spec_rl = 1; % Measure the calls recieved level relative to background noise




%% Preliminary Filtration Switches and Values:

% On/off: Merge calls within time bin radius of another call (handle
% intra-call energy drops).
parm.filter_parm.switches.AdjacentCallMerger = 1;

% Set bin number for merging:
parm.filter_parm.values.AdjacentCallBinNum = 3;


    % % Contour measurement value filters:
    % % Contour Duration:
    % parm.filter_parm.values.cm_min_duration_s = 0.3; % Minimum time length (seconds) for cm
    % parm.filter_parm.values.cm_max_min_duration_s = 0.3; % " " for cm_max
    % parm.filter_parm.values.cm_max2_min_duration_s = 1; % " " for cm_max2
    % 
    % % Contour Frequency Bandwidth:
    % parm.filter_parm.values.cm_min_freq_bandwidth_hz = 3; % Minimum frequency bandwidth (Hz) for cm
    % parm.filter_parm.values.cm_max_min_freq_bandwidth_hz = 3; % " " for cm_max
    % parm.filter_parm.values.cm_max2_min_freq_bandwidth_hz = 3; % " " cm_max2
    % 
    % % Contour Frequncy start:
    % parm.filter_parm.values.cm_start_freq_hz = 55;
    % parm.filter_parm.values.cm_max_start_freq_hz = 55;
    % parm.filter_parm.values.cm_max2_start_freq_hz = 55;
    % 
    % % Contour Frequency End:
    % parm.filter_parm.values.cm_end_freq_hz = 40;
    % parm.filter_parm.values.cm_max_end_freq_hz = 40;
    % parm.filter_parm.values.cm_max2_end_freq_hz = 40;
    % 
    % % Contour Minimum Frequency:
    % parm.filter_parm.values.cm_min_freq_hz = 40;
    % parm.filter_parm.values.cm_max_min_freq_hz = 40;
    % parm.filter_parm.values.cm_max2_min_freq_hz = 40;
    % 
    % % Contour Maximum Frequency:
    % parm.filter_parm.values.cm_max_freq_hz = 55;
    % parm.filter_parm.values.cm_max_max_freq_hz = 55;
    % parm.filter_parm.values.cm_max2_max_freq_hz = 55;
    % 
    % % Contour slope range:
    % parm.filter_parm.values.cm_slope_lower = -20; % Minimum slope for cm 
    % parm.filter_parm.values.cm_max_slope_lower = -20; % " ' cm_max
    % parm.filter_parm.values.cm_max2_slope_lower = -20; % " " cm_max2
    % parm.filter_parm.values.cm_slope_upper = -5; % Maximum slope for cm
    % parm.filter_parm.values.cm_max_slope_upper = -5; % " " for cm_max
    % parm.filter_parm.values.cm_max2_slope_upper = -5; % " " cm_max2


% % Plot the detection window (will stop the detection process after each
% % window containing a call until user manually progresses.
% parm.plot = 0;


% Island radius: Number of bins to check around islands for others (connect
% the dots for contours)
parm.island_radius = 1; 

% % Max file size: When the output structure 'hyd' reaches a maximum size, it
% % will be saved as a .mat file and a new structure will be started. 
% parm.output_size_separation = 1;
% parm.max_file_size = 1; % Max file size in GB
% 
% % Plot detection frequency versus deployment date. 
% parm.det_frequency_plot = 1;
% 
% % Identify and plot the amount of time spent locating detections in each
% % window
% parm.call_time_tracker = 1;

% Deck test offset: The index of the first timestamp that the recording
% actually starts
parm.pre_offset = 6;




%% GPL Fetch Parameters

% Padding on GPL call extraction. Number of extra bins around the extracted
% call to include.
% parm.fetch.extract_gpl_padding = 5;

% Turn time delay on:1, off:0 (inconclusive reasoning)
parm.fetch.time_delay_switch = 1;

% Turn legends on:1, off:0 for plots
parm.fp1.legend = 0; 

% Save spectrogram/whitened spectrogram/GPL contour window in output, on:1, off:0
% Note that the spectrogram and contour window are plotted regardless of
% the values of these two switches, they just turn on/off saving to output
% file.
parm.fp1.spectrogram = 0;
parm.fp1.whitened_sp = 0;
parm.fp1.gpl_contour_window = 0;

% Wait time (seconds) of yellow confirm box
parm.fp1.box_wt = 0.5;

% Height the detection numbers are off of the contour subplot
% parm.fp1.numHe2ight = 3;

% Automatically full-screen the figure? (1:on,0:off)
parm.fp1.fullscreen = 0;

% Switch second subplot from default (GPL Energy Summation) to Cm-max
% Contours (Single Strongest island of each contour) (1:on,0:off)
parm.fp1.cm_maxsubplot = 1;
parm.fp1.savecm_max = 0; % Save cm_max matrix to output (1:on,0:off)

% Create main spectrogram using the methods from Triton: This allows the
% application of a transfer function and brightness/contrast control. Set
% brightness=0,contrast=100 for no change.
% (0:off,1:on)
parm.fp1.triton_sp = 1; % Triton Spectrogram
parm.fp1.tf.switch = 0; % Apply Transfer function
parm.fp1.brightness = 30; % Brightness (dB)
parm.fp1.contrast = 645; % Contrast (%dB) (In xx% form e.g. 125 for 125%)
parm.fp1.colorbar.switch = 1; % Turn colorbar on/off
parm.fp1.colormap = 'jet'; % Colormap Type
parm.fp1.sp_det_boxes = 1; % Turn on/off GPL detection boxes on spectrogram
parm.fp1.noise_floor = -40; % Attenuation level [dB] 


% Verify outputs one by one. This will plot the paired detections contour
% and provide measurements. (1:on,0:off)
parm.fp1.outputVerify = 0;

% Fetch P2P and rms RL measurements for paired calls and adhocs
% Recommended to leave all 5 RL measurement methods on for post-processing
% (1:on,0:off)
parm.fp1.p2p_rms_RL = 1;
parm.fp1.rl_time_fix = 1; % Apply time adjustment to timestamps analytically determined by %OL and FFTL
parm.fp1.rms_super_window_length = 3; % Size of the window of which the call is in the center. Value is a multiple of total duration of the call added to both ends. (e.g. '1.5' would add 1.5*T samples to the beginning and end, where T is the duration of the call)
parm.fp1.rms_time_padding = 0.5; % noise padding for rms window. Value is a multiplier of call length. 
parm.fp1.bandpass = 1; % Apply bandpass filter to RL measurements?
parm.fp1.bandpass_corners = 2; % Amount of Hz above & below the contours min/max frequency that the waveform will be bandlimited to (e.g. '2' for a 15-30Hz contour will make the passband 13-32Hz)
parm.fp1.single_tf_val = 1; % Use a single TF value at the peak frequency for conversion to uPa (No other method currently supported)
parm.fp1.rl.m3db = 1; % Use -3dB method to set RL start/end time (Madsen 2005)
parm.fp1.rl.m10db = 1; % ... -10dB method (Madsen 2005)
parm.fp1.rl.e90 = 1; % ... 90% cumulative energy method (Madsen 2005)
parm.fp1.rl.e97 = 1; % ... 97% cumulative energy method (Madsen 2005)
parm.fp1.rl.gpl_SE = 1; % Use GPL start/end times to find RL measurements

% Apply time delay sample calculation
if parm.fp1.rl_time_fix == 1
    parm.fp1.tdelay = parm.fftOverlapPercent*1e-2*parm.fftl/2;
end


%% Save Output

save('Bb90_30_parm.mat','parm');

end % Function
