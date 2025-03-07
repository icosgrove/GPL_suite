function [parm] = GPL_parameter_input_v3(PARAMS,xwav_struct)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPL_parameter_input_v3 is a script to set fixed input parameters for the
% GPL detector algorithm. Some parameters are deterined by xwav header 
% information, others are manual input Some of these parameters will change 
% during the GPL process. 

% Written by Ian Cosgrove 01/26/2024
% Lsst Updated 04/03/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Multiplier for Tyler's spectrogram padding. Set to ~=0 to use
% Tyler's original spectrogram padding process. The value of the variable
% will determine how many seconds are padded around each spectrogram.
% Tyler's sample range is recommended to use, since it won't skip any
% samples.
parm.pad = 1;

% Sample Frequency [samples/s]: Value straight from xwav
parm.SampleFreq = PARAMS.ltsahd.sample_rate(1);

% Number of seconds to be processed at one time [seconds]:
parm.window_length = 75;

% Number of samples per wav 
if xwav_struct.xwav1.nwav == 1 % xwav downloaded directly from Triton, only one wav
    parm.nrec = parm.window_length*parm.SampleFreq; 
else % xwav data in general, already split up every 75s by disk writes
    parm.nrec = xwav_struct.xwav1.TotalSamples/xwav_struct.xwav1.nwav;
end


% FIFO Removal switch: Remove n*50Hz bins and replace them with average of
% the adacent bins. Done to remove 50Hz harmonics.
parm.fifo_removal = 0;
parm.fund_fifo_freq = 50; % Fundamental frequency of FIFO noise, rest are harmonics (integer multiples)
parm.fifo3bin = 1; % Average and replace the n*(49:51) frequency bins 


% Helble et al 2012 eq 6 (nu1 and nu2)
% GPL Algorithm Parameters
parm.xp1 = 2.5;     % Exponent for u: Normed over frequency
parm.xp2 = 0.5;     % Exponent for y: Normed over time
parm.whiten = 1;    % 1: Whiten in time only. 2: Whiten in time and frequency
parm.white_x = 1.5; % Enhanced mean noise multiplier
parm.NumLoops = 5;  % max number of loops each window is processed for calls



% Frequency Range [Hz]
parm.freq_lo = 30;
parm.freq_hi = 90;

% Call duration bounds [s]
parm.MinCallDuration = 0.5;
parm.MaxCallDuration = 10;

% Spectral Parameters
parm.fftl = 2000; % FFT Length
parm.fftOverlapPercent = 90; % FFT overlap percentage (xx%, e.g. '50' for 50%)
if parm.fftOverlapPercent == 0
    parm.fftOverlap = parm.fftl;
else
    parm.fftOverlap = parm.fftl - (parm.fftOverlapPercent*1e-2)*parm.fftl; % Calculate number of samples to step forward based on overlap %
end

% Use future samples to create spectrogram? This will index forward to fill
% in the overlap. Otherwise use Tyler's sample extraction (also must have
% parm.pad=1)
parm.future_sample_ext = 1;
if parm.future_sample_ext == 1
    parm.pad = 0;
end


% Time bin parameters:
if parm.future_sample_ext == 1
    parm.NumTimeBins = floor(parm.nrec/parm.fftOverlap);
else
    parm.NumTimeBins = (parm.nrec-parm.fftl)/parm.fftOverlap + 1;
end
    

% Frequency bin parameters:
freq_range = 0:parm.SampleFreq/parm.fftl:(parm.fftl/2)/parm.fftl*parm.SampleFreq; % Define frequency range from FFTL

[~,parm.FreqBinLo] = min(abs(freq_range - parm.freq_lo));
[~,parm.FreqBinHi] = min(abs(freq_range - parm.freq_hi));
parm.NumFreqBins = parm.FreqBinHi - parm.FreqBinLo + 1;

% Frequency bin range for summation
parm.SumFreqBinLo = parm.FreqBinLo;
parm.SumFreqBinHi = parm.FreqBinHi;

% Detection/Contour Threshold Parameters
% 1-D
parm.noise_ceiling = 20; % Upper limit for backgorund noise, multiple of noise floor
parm.thresh = 90; % Call threhsold: minimum multiple of noise floor for a call to pass

% 2-D
parm.ContourCutoff = 6; % Minimum value for time,frequency coordinate to be included to the call contour
parm.MinIslandSize = 10; % Minimum number of pixels above cutoff for a detection to be considered a contour
parm.smaller_island_cutoff = 0.05; % All islands smaller than the largest island must be at least this percentage size of the largest island.


%% Measurement On/Off Switches
parm.waveform_on = 0; % Extract the original time signal waveform of the call
parm.template_on = 1; % Extract call contour and measurements (0: only start/end time is recorded)
parm.cm_on = 1; % Extarct the call contour 
parm.cm_max_on = 1; % Extract the strongesr contour of the call
parm.cm_max2_on = 1; % Extract the second strongest contour of the call
parm.measurements_on = 1; % Obtain measurements of the contours
parm.filter_on = 1; % Run primary filters 


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


    % Contour measurement value filters:
    % Contour Duration:
    parm.filter_parm.values.cm_min_duration_s = 0.3; % Minimum time length (seconds) for cm
    parm.filter_parm.values.cm_max_min_duration_s = 0.3; % " " for cm_max
    parm.filter_parm.values.cm_max2_min_duration_s = 1; % " " for cm_max2

    % Contour Frequency Bandwidth:
    parm.filter_parm.values.cm_min_freq_bandwidth_hz = 3; % Minimum frequency bandwidth (Hz) for cm
    parm.filter_parm.values.cm_max_min_freq_bandwidth_hz = 3; % " " for cm_max
    parm.filter_parm.values.cm_max2_min_freq_bandwidth_hz = 3; % " " cm_max2
    
    % Contour Frequncy start:
    parm.filter_parm.values.cm_start_freq_hz = 55;
    parm.filter_parm.values.cm_max_start_freq_hz = 55;
    parm.filter_parm.values.cm_max2_start_freq_hz = 55;
    
    % Contour Frequency End:
    parm.filter_parm.values.cm_end_freq_hz = 40;
    parm.filter_parm.values.cm_max_end_freq_hz = 40;
    parm.filter_parm.values.cm_max2_end_freq_hz = 40;
    
    % Contour Minimum Frequency:
    parm.filter_parm.values.cm_min_freq_hz = 40;
    parm.filter_parm.values.cm_max_min_freq_hz = 40;
    parm.filter_parm.values.cm_max2_min_freq_hz = 40;
    
    % Contour Maximum Frequency:
    parm.filter_parm.values.cm_max_freq_hz = 55;
    parm.filter_parm.values.cm_max_max_freq_hz = 55;
    parm.filter_parm.values.cm_max2_max_freq_hz = 55;

    % Contour slope range:
    parm.filter_parm.values.cm_slope_lower = -20; % Minimum slope for cm 
    parm.filter_parm.values.cm_max_slope_lower = -20; % " ' cm_max
    parm.filter_parm.values.cm_max2_slope_lower = -20; % " " cm_max2
    parm.filter_parm.values.cm_slope_upper = -5; % Maximum slope for cm
    parm.filter_parm.values.cm_max_slope_upper = -5; % " " for cm_max
    parm.filter_parm.values.cm_max2_slope_upper = -5; % " " cm_max2


% Plot the detection window (will stop the detection process after each
% window containing a call until user manually progresses.
parm.plot = 0;


% Island radius: Number of bins to check around islands for others (connect
% the dots for contours)
parm.island_radius = 1; 

% Max file size: When the output structure 'hyd' reaches a maximum size, it
% will be saved as a .mat file and a new structure will be started. 
parm.output_size_separation = 1;
parm.max_file_size = 0.5; % Max file size in GB

% Plot detection frequency versus deployment date. 
parm.det_frequency_plot = 1;

% Identify and plot the amount of time spent locating detections in each
% window
parm.call_time_tracker = 1;

% Apply deck test adjustment for timestamps
parm.deck_test_adjustment = 1;



%% Save Output

%save('sei.mat','parm');

end % Function

