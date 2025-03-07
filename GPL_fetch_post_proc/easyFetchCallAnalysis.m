% Statistical Analysis for GPL Fetch output calls.
% Written: Ian 09/2024

clear; clc; close

%% User Setup
% Load file with all calls
load('all_bb_90_fetch_SB_03.mat');
tf = load('976_210816_B_HARP.tf');

% Turn on/off histogram plots. Navigate to histogram section to specify
% parameters.
param.hist_on = 1;

% Turn on/off stats table. Navigate to stats table section to specify metrics
param.stat_table = 1;

% Turn on/off random example figure generation (requires stat table on)
param.fig_gen = 0;



%% Processing

% Separate into ready, extra, reprocess bins
[plot_calls, adhoc_off, reprocess_calls] = bin1Calls(hyd_fetch.calls,hyd_fetch.adhoc_detections);

% Extract bulk measurements
[STATS] = unloadStats(plot_calls);

%% Measurement Histograms 
% Figure # start
v=10; 

if param.hist_on == 1
    % Plot Duration/Slope/Frequency Bandwidth: cm
    histPARAMS.cm.dur.binWidth = 0.1;
    histPARAMS.cm.slope.binWidth = 1;
    histPARAMS.cm.nonzFB.binWidth = 1;
    histPARAMS.cm.fullFB.binWidth = 1;
    [v] = plotDSFB_cm(STATS,v,histPARAMS);
    
    
    % Plot Duration/Slope/Frequency Bandwidth: cm_max
    histPARAMS.cm_max.dur.binWidth = 0.1;
    histPARAMS.cm_max.slope.binWidth = 1;
    histPARAMS.cm_max.nonzFB.binWidth = 1;
    histPARAMS.cm_max.fullFB.binWidth = 1;
    [v] = plotDSFB_cm_max(STATS,v,histPARAMS);
    
    
    % Plot start/end/min/max frequency: cm
    histPARAMS.cm.minf.binWidth = 1;
    histPARAMS.cm.maxf.binWidth = 1;
    histPARAMS.cm.startf.binWidth = 1;
    histPARAMS.cm.endf.binWidth = 1;
    [v] = plotSEMM_cm(STATS,v,histPARAMS);
    
    
    % Plot start/end/min/max frequency: cm_max
    histPARAMS.cm_max.minf.binWidth = 1;
    histPARAMS.cm_max.maxf.binWidth = 1;
    histPARAMS.cm_max.startf.binWidth = 1;
    histPARAMS.cm_max.endf.binWidth = 1;
    [v] = plotSEMM_cm_max(STATS,v,histPARAMS);
    
    
    % Plot Peak Frequency & P2P RL (cm & cm_max)
    histPARAMS.cm.pf.binWidth = 1;
    histPARAMS.cm_max.pf.binWidth = 1;
    histPARAMS.cm.p2pRL.binWidth = 1;
    histPARAMS.cm_max.p2pRL.binWidth = 1;
    [v] = plotPFP2P(STATS,v,histPARAMS);
    
    
    % Plot RMS RL (cm)
    histPARAMS.cm.rms3.binWidth = 1;
    histPARAMS.cm.rms10.binWidth = 1;
    histPARAMS.cm.rms90.binWidth = 1;
    histPARAMS.cm.rms97.binWidth = 1;
    histPARAMS.cm.rmsGPL.binWidth = 1;
    [v] = plotrmsRL_cm(STATS,v,histPARAMS);
    
    
    % Plot RMS RL (cm_max)
    histPARAMS.cm_max.rms3.binWidth = 1;
    histPARAMS.cm_max.rms10.binWidth = 1;
    histPARAMS.cm_max.rms90.binWidth = 1;
    histPARAMS.cm_max.rms97.binWidth = 1;
    histPARAMS.cm_max.rmsGPL.binWidth = 1;
    [v] = plotrmsRL_cm_max(STATS,v,histPARAMS);
    
    
    % Plot SEL (cm)
    histPARAMS.cm.sel3.binWidth = 1;
    histPARAMS.cm.sel10.binWidth = 1;
    histPARAMS.cm.sel90.binWidth = 1;
    histPARAMS.cm.sel97.binWidth = 1;
    histPARAMS.cm.selGPL.binWidth = 1;
    [v] = plotSEL_cm(STATS,v,histPARAMS);
    
    
    % Plot SEL (cm_max)
    histPARAMS.cm_max.sel3.binWidth = 1;
    histPARAMS.cm_max.sel10.binWidth = 1;
    histPARAMS.cm_max.sel90.binWidth = 1;
    histPARAMS.cm_max.sel97.binWidth = 1;
    histPARAMS.cm_max.selGPL.binWidth = 1;
    [v] = plotSEL_cm_max(STATS,v,histPARAMS);
end

rawSTATS = STATS;


%% Stats Table

if param.stat_table == 1

    % Remove NaN's from STATS
    [STATS] = removeNaN(STATS);
    
    % Find mean, standard deviation, and range for each metric (cm & cm_max)
    % Switches for which measurements to include
    tablePARAMS.cm.dur_s = 1;
    tablePARAMS.cm.slope = 1;
    tablePARAMS.cm.band_hz = 1;
    tablePARAMS.cm.start_hz = 1;
    tablePARAMS.cm.end_hz = 1;
    tablePARAMS.cm.min_hz = 1;
    tablePARAMS.cm.max_hz = 1;
    tablePARAMS.cm.abs_band_hz = 1;
    tablePARAMS.cm.peak_freq = 1;
    tablePARAMS.cm.p2p_rl = 1;
    tablePARAMS.cm.rms_rl_m3 = 1;
    tablePARAMS.cm.rms_rl_m10 = 0;
    tablePARAMS.cm.rms_rl_e90 = 0;
    tablePARAMS.cm.rms_rl_e97 = 0;
    tablePARAMS.cm.rms_rl_gplSE = 0;
    tablePARAMS.cm.SEL_m3 = 1;
    tablePARAMS.cm.SEL_m10 = 0;
    tablePARAMS.cm.SEL_e90 = 0;
    tablePARAMS.cm.SEL_e97 = 0;
    tablePARAMS.cm.SEL_gplSE = 0;
    tablePARAMS.cm_max.dur_s = 1;
    tablePARAMS.cm_max.slope = 1;
    tablePARAMS.cm_max.band_hz = 1;
    tablePARAMS.cm_max.start_hz = 1;
    tablePARAMS.cm_max.end_hz = 1;
    tablePARAMS.cm_max.min_hz = 1;
    tablePARAMS.cm_max.max_hz = 1;
    tablePARAMS.cm_max.abs_band_hz = 1;
    tablePARAMS.cm_max.peak_freq = 1;
    tablePARAMS.cm_max.p2p_rl = 1;
    tablePARAMS.cm_max.rms_rl_m3 = 1;
    tablePARAMS.cm_max.rms_rl_m10 = 0;
    tablePARAMS.cm_max.rms_rl_e90 = 0;
    tablePARAMS.cm_max.rms_rl_e97 = 0;
    tablePARAMS.cm_max.rms_rl_gplSE = 0;
    tablePARAMS.cm_max.SEL_m3 = 1;
    tablePARAMS.cm_max.SEL_m10 = 0;
    tablePARAMS.cm_max.SEL_e90 = 0;
    tablePARAMS.cm_max.SEL_e97 = 0;
    tablePARAMS.cm_max.SEL_gplSE = 0;
    [TABLE,cm_table,cm_max_table] = makeTable(STATS,tablePARAMS);

end


%% Random Example Generation

if param.fig_gen == 1

    % Parameter Setup
    figPARAMS.numFigs = 2; % Number of random examples to generate
    figPARAMS.subplotOn = 1;
    figPARAMS.subplotSize = [2 1]; % Set to [1 1] for single figure plots
    figPARAMS.padding = 2; % [s] Space before and after the call

    figPARAMS.sampleFreq = 2000;
    figPARAMS.fftl = 2500;
    figPARAMS.overlapPerc = 95;
    figPARAMS.window = 'hanning';
    figPARAMS.freqLo = 10;
    figPARAMS.freqHi = 100;
    figPARAMS.brightness = 40;
    figPARAMS.contrast = 500;
    figPARAMS.colormap = 'jet';
    figPARAMS.tf_on = 0;
    figPARAMS.tf = tf;
    figPARAMS.colorbar.on = 1;

    figPARAMS.xlabel = 'Time [s]';
    figPARAMS.labelLoc = 'All'; % Options: 'All' , 'One' , 'Outside'
    figPARAMS.xFontSize = 10;
    figPARAMS.xlabelBold  = 1;
    figPARAMS.y1label = 'Frequency [Hz]';
    figPARAMS.y1FontSize = 10;
    figPARAMS.y1labelBold  = 1;
    figPARAMS.y2FontSize = 10;
    figPARAMS.y2labelBold  = 1;
    figPARAMS.title = 'Test';
    figPARAMS.supTitle = 'Super Test';
    figPARAMS.box = 1;
    figPARAMS.grid = 1;
    figPARAMS.tickDir = 'out'; % 'out' or 'in'

    % Filter params 
    figPARAMS.filter.on = 1;

    % Metrics & ranges. Copy and paste the 2 lines for metric & range 
    % and add on more metrics with respective ranges. Adjust the index of the
    % figPARAMS.filter.metric/stdDevRange as needed. If not they will
    % overwrite each other. Options for metrics in Readme. For < or > set
    % lower limit to 0 or upper limit to inf (regardless of sign on numeric
    % term)
    figPARAMS.filter.metric{1} = 'RMS RL [-3dB] [cm]';
    figPARAMS.filter.stdDevRange{1} = [1 1.5]; 
    figPARAMS.filter.metric{2} = 'Slope [cm_max]';
    figPARAMS.filter.stdDevRange{2} = [1 1.5];

    % Filter for calls to plot
    [figCalls] = filterCalls(figPARAMS,plot_calls,rawSTATS);

    % Pick n random from the filtered pile
    if length(figCalls) >= figPARAMS.numFigs
        figCalls = figCalls(randi([1,length(figCalls)],1,figPARAMS.numFigs));
    else
        figCalls = figCalls(randi([1,length(figCalls)],1,length(figCalls)));
        fprintf('\nNumber of calls plotted limited to %i (originally set to %i)\n',length(figCalls),figPARAMS.numFigs);
    end

    % Load xwavs
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
    parm = hyd_fetch.parm;

    % Configure plotting data
    if figPARAMS.subplotOn == 1
        figPARAMS.subsPerFig = (figPARAMS.subplotSize(1)*figPARAMS.subplotSize(2));
        figPARAMS.numTotalFigs = floor(figPARAMS.numFigs / figPARAMS.subsPerFig);
        if length(figCalls) < (figPARAMS.subsPerFig*figPARAMS.numTotalFigs)
            fprintf('\nRequesting %i examples, but there are only %i possible with the given filters\n',(figPARAMS.subsPerFig*figPARAMS.numTotalFigs),length(figCalls))
            return
        end
        for n = 1:figPARAMS.numTotalFigs
            fig = figure(n);
            [figPARAMS.row,figPARAMS.col] = getSPIndex(figPARAMS);
            if ~isempty(figPARAMS.supTitle)
                sgtitle(figPARAMS.supTitle)
            end
            if strcmp(figPARAMS.labelLoc,'One')
                t = tiledlayout(figPARAMS.subplotSize(1),figPARAMS.subplotSize(2));
            else
                t = [];
            end
            for n1 = 1:figPARAMS.subsPerFig
                
                if strcmp(figPARAMS.labelLoc,'One')
                    nexttile(t,n1)
                end

                % Find sample start/end in xwav
                % Start/end julian time
                st_ju = figCalls(n1).julian_start_time;
                en_ju = figCalls(n1).julian_end_time;
    
                % Find xwav file 
                startLoc = find(sort([xwav_start st_ju])==st_ju);
                endLoc = find(sort([xwav_start en_ju])==en_ju);
                if startLoc ~= endLoc
                    disp('Call crosses xwav, remove it and try again') % Should be impossible
                    fprintf('%i\n',n)
                    return
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
                    disp('Call crossed windows, remove it and try again') % Should be impossible at this stage
                    fprintf('%i\n',n)
                    return
                else
                    start_fix = sample_range(start_fix);
                    end_fix = sample_range(end_fix+1)-1;
                end

                % Pre-processing on data chunk
                stime = figCalls(n1).start_time + start_fix - (figPARAMS.padding*figPARAMS.sampleFreq);
                endtime = figCalls(n1).end_time + start_fix + (figPARAMS.padding*figPARAMS.sampleFreq);
                if parm.fp1.rl_time_fix == 1 % Apply time delay adjustment
                    stime = stime + parm.fp1.tdelay;
                    endtime = endtime + parm.fp1.tdelay;
                end
                % Check endpoints
                if stime < 1 % ensure no 0/negative samples
                    stime = 1;
                end
                if endtime > xwav_struct.(file_string).TotalSamples % ..and no over-indexing
                    endtime = xwav_struct.(file_string).TotalSamples;
                end
                call_data = double(audioread(files(file_ident).name,...
                    [stime,endtime],'native'));
                [freq,time,sp,figPARAMS] = makeSp(call_data,figPARAMS); % Create spectrogram
                freq = freq(figPARAMS.freq_loc_min:figPARAMS.freq_loc_max);
                
                % Add to plot
                if ~strcmp(figPARAMS.labelLoc,'One')
                    subplot(figPARAMS.subplotSize(1),figPARAMS.subplotSize(2),n1)
                end
                image(time,freq,sp)
                set(gca,'YDir','normal') % y-tick value
                setLabel(figPARAMS,'x',n1,fig,t) % x-axis label
                setLabel(figPARAMS,'y1',n1,fig,t) % y-axis label
                setLabel(figPARAMS,'y2',n1,fig,t) % second y-axis label
                ax = gca;
                ax.TickDir = figPARAMS.tickDir; % Tick direction
                colormap(figPARAMS.colormap) % colormap
                if ~isempty(figPARAMS.title) % ind. title
                    title(figPARAMS.title)
                end

            
            end

        end

    end





end



































