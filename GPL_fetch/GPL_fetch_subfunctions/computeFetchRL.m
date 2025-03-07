function [p2p_RL,rms_RL,sel,peak_freq,RL_note] = computeFetchRL(cm,parm,data,dur)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the received levels of specific contours. It is
% used for cm,cm_max,cm_max2 in the same way, just different sized
% contours.
% Written: Ian 08/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allocate
RL_note = [];
p2p_RL = [];
rms_RL.m3db= [];
rms_RL.m10db = [];
rms_RL.e90 = [];
rms_RL.e97 = [];
rms_RL.gpl_SE = [];
sel.m3db= [];
sel.m10db = [];
sel.e90 = [];
sel.e97 = [];
sel.gpl_SE = [];
peak_freq = [];

% Find corner frequencies for bandpass of the given contour
parm.bp_max_f = parm.freq_hi + parm.fp1.bandpass_corners;
parm.bp_min_f = parm.freq_lo - parm.fp1.bandpass_corners;

% Ensure the corners fall within valid frequency range
if parm.bp_max_f > parm.SampleFreq / 2 % Limit to Nyquist frequency
    parm.bp_max_f = parm.SampleFreq / 2;
end
if parm.bp_min_f < 1 % Don't do 0 or negative frequencies
    parm.bp_min_f = 1;
end

% Find the peak frequency of the waveform using methods from Triton
peak_data = data(parm.fp1.rms_super_window_length*dur+1 - parm.fp1.rms_time_padding*dur:parm.fp1.rms_super_window_length*dur+dur+parm.fp1.rms_time_padding*dur);
[peak_freq] = computePeakFreq(peak_data,parm);

% Find TF value at peak frequency
if ~isfield(parm.fp1.tf,'file')
    p2p_RL = nan;
    rms_RL = nan;
    sel = nan;
    disp('No TF file loaded, no RL values calculated.')
    return
else
    tf = load(parm.fp1.tf.file); % Load
    parm.fp1.tf.f = tf(:,1); % frequency
    parm.fp1.tf.vals = tf(:,2); % dB re 1 uPa/count
end

if sum(peak_freq == parm.fp1.tf.f) > 0
    tf_val = parm.fp1.tf.vals(find(peak_freq==parm.fp1.tf.f));
else
    loc = find(sort([peak_freq; parm.fp1.tf.f]) == peak_freq);
    tf_val = mean([parm.fp1.tf.vals(loc-1) parm.fp1.tf.vals(loc)]);
end

% ft = fft(data,2000);
% figure(50); %plot(10*log10(abs(ft)))
% semilogx((abs(ft)))

% Apply elliptical bandpass filter to the frequency range
if parm.fp1.bandpass == 1
    [b,a] = ellip(4,0.1,40,[parm.bp_min_f parm.bp_max_f]*2/parm.SampleFreq);
    data = filtfilt(b,a,data);
end

% Restrict window to call +/- small noise margins
% call = data(parm.fp1.rms_super_window_length*dur+1:parm.fp1.rms_super_window_length*dur+dur);
data = data(parm.fp1.rms_super_window_length*dur+1 - parm.fp1.rms_time_padding*dur:parm.fp1.rms_super_window_length*dur+dur+parm.fp1.rms_time_padding*dur);
st = parm.fp1.rms_time_padding*dur+1;
fi = length(data) - parm.fp1.rms_time_padding*dur;

% Find P2P RL
p2p_counts = abs(max(data)) + abs(min(data));
max_loc = find(max(data)==data);
min_loc = find(min(data)==data);
if ((max_loc < st) || (max_loc > fi)) || ((min_loc < st) || (min_loc > fi))
    p2p_RL = nan; % Outside of contour
    RL_note.p2p = 'Min or max waveform value outside of contour';
else
    p2p_RL = 10*log10(p2p_counts); % dB re 1 count
    p2p_RL = p2p_RL + tf_val; % dB re 1 uPa
end

% Compute rms/SEL RL's using selected method(s)
% Envelope cutoffs using -xdBs
if (parm.fp1.rl.m3db == 1) || (parm.fp1.rl.m10db == 1) 

    % Create amplitude envelope using Hilbert Transform
    analytic_signal = abs(hilbert(data));

    % Find max & location
    peak_amp = max(analytic_signal);
    peak_loc = find(analytic_signal == peak_amp);

    % If peak amplitude is outside of start/end, invalid measurement
    if (peak_loc < st) || (peak_loc > fi)
        if parm.fp1.rl.m3db == 1
            rms_RL.m3db = nan;
            sel.m3db = nan;
            RL_note.m3db = 'Peak amplitude outside of contour';
        end
        if parm.fp1.rl.m10db == 1
            rms_RL.m10db = nan;
            sel.m10db = nan;
            RL_note.m10db = 'Peak amplitude outside of contour';
        end

    else % Proceed
        
        % -3dB RMS & SEL RL
        if parm.fp1.rl.m3db == 1
            m3dbAmp = 10^(-3/20)*peak_amp;
            t1 = find(analytic_signal(1:peak_loc) <= m3dbAmp,1,'last');
            t2 = find(analytic_signal(peak_loc:end) <= m3dbAmp,1,'first') + length(analytic_signal(1:peak_loc));
            if (t1 < st) || (t2 > fi) % markers outside of contour
                rms_RL.m3db = nan;
                sel.m3db = nan;
                RL_note.m3db = '-3dB timestamps fell outside of contour';
            else % Proceed
                [rms_RL.m3db,sel.m3db] = computeRMSandSEL(data(t1:t2),parm);
                rms_RL.m3db = rms_RL.m3db + tf_val; % dB re 1 uPa
                sel.m3db = sel.m3db + tf_val; %
            end
        end

        % -10dB RMS & SEL RL
        if parm.fp1.rl.m10db == 1
            m10dbAmp = 10^(-10/20)*peak_amp;
            t3 = find(analytic_signal(1:peak_loc) <= m10dbAmp,1,'last');
            t4 = find(analytic_signal(peak_loc:end) <= m10dbAmp,1,'first') + length(analytic_signal(1:peak_loc));
            if (t1 < st) || (t2 > fi) % markers outside of contour
                rms_RL.m10db = nan;
                sel.m10db = nan;
                RL_note.m10db = '-10dB timestamps fell outside of contour';
            else % Proceed
                [rms_RL.m10db,sel.m10db] = computeRMSandSEL(data(t3:t4),parm);
                rms_RL.m10db = rms_RL.m10db + tf_val; % dB re 1 uPa
                sel.m10db = sel.m10db + tf_val; % 
            end
        end
    end
end

% Cumulative energy sum using x% of total energy
if (parm.fp1.rl.e90 == 1) || (parm.fp1.rl.e97)

    t = 1:length(data); % Time vector 
    e = data.^2; % Energy 
    csum = cumtrapz(t,e); % energy integral
    total_energy = trapz(t,e);
    csum = csum ./ total_energy; % Normalized for %

    % 90% Total energy
    if parm.fp1.rl.e90 == 1  
        e90_s = find(csum >= 0.05,1,'first');
        e90_f = find(csum <= 0.95,1,'last');
        if (e90_s < st) || (e90_f > fi)
            RL_note.e90 = 'Timestamp(s) are outside of call';
        end
        [rms_RL.e90,sel.e90] = computeRMSandSEL(data(e90_s:e90_f),parm);
        rms_RL.e90 = rms_RL.e90 + tf_val; % dB re 1 uPa
        sel.e90 = sel.e90 + tf_val; % 
    end

    % 97% Total energy
    if parm.fp1.rl.e97 == 1
        e97_s = find(csum >= 0.015,1,'first');
        e97_f = find(csum <= 0.985,1,'last');
        if (e97_s < st) || (e97_f > fi)
            RL_note.e97 = 'Timestamp(s) are outside of call';
        end
        [rms_RL.e97,sel.e97] = computeRMSandSEL(data(e97_s:e97_f),parm);
        rms_RL.e97 = rms_RL.e97 + tf_val; % dB re 1 uPa
        sel.e97 = sel.e97 + tf_val; % 
    end
end

% Find rms RL and SEL using GPL start/end markers
if parm.fp1.rl.gpl_SE == 1
    [rms_RL.gpl_SE,sel.gpl_SE] = computeRMSandSEL(data(st:fi),parm);
    rms_RL.gpl_SE = rms_RL.gpl_SE + tf_val; % dB re 1 uPa
    sel.gpl_SE = sel.gpl_SE + tf_val;
end


% Test Plotting

% figure(11); 
% subplot(2,1,1); hold on
% plot(data)
% xlabel('Time [seconds]'); ylabel('Counts');
% new_xticks = xticks/2000;
% xticklabels(arrayfun(@(x) num2str(x, '%.1f'), new_xticks, 'UniformOutput', false));
% box on; grid on
% plot([st st],[min(data) max(data)],'-r')
% plot([fi fi],[min(data) max(data)],'-r')
% title(sprintf('Waveform Bandpassed: [%i-%iHz]',parm.bp_min_f,parm.bp_max_f))
% axis tight
% p1 = plot(nan,nan,'-r');
% legend(p1,'Call Start/End')
% hold off
% 
% subplot(2,1,2); hold on
% plot(csum)
% plot([e90_s e90_s],[0 1],'-k')
% plot([e90_f e90_f],[0 1],'-k')
% plot([0 e90_s],[0.05 0.05],'--k')
% plot([0 e90_f],[0.95 0.95],'--k')
% plot([e97_s e97_s],[0 1],'-k')
% plot([e97_f e97_f],[0 1],'-k')
% plot([0 e97_s],[0.015 0.015],'--k')
% plot([0 e97_f],[0.985 0.985],'--k')
% plot([st st],[0 1],'-r')
% plot([fi fi],[0 1],'-r')
% xlabel('Time [seconds]'); ylabel('Energy %');
% new_xticks = xticks/2000;
% xticklabels(arrayfun(@(x) num2str(x, '%.1f'), new_xticks, 'UniformOutput', false));
% title('Relative Cumulative Energy')
% box on; grid on
% axis tight
% hold off
% 
% 
% figure(9); subplot(2,1,1); hold on
% plot(data)
% xlabel('Time [seconds]'); ylabel('Counts');
% new_xticks = xticks/2000;
% xticklabels(arrayfun(@(x) num2str(x, '%.1f'), new_xticks, 'UniformOutput', false));
% box on; grid on
% plot([st st],[min(data) max(data)],'-r')
% plot([fi fi],[min(data) max(data)],'-r')
% title(sprintf('Waveform Bandpassed: [%i-%iHz]',parm.bp_min_f,parm.bp_max_f))
% axis tight
% p1 = plot(nan,nan,'-r');
% legend(p1,'Call Start/End')
% hold off
% 
% subplot(2,1,2); hold on
% plot(analytic_signal); xlim([0 length(data)]); title(sprintf('Amplitude Envelope'))
% xlabel('Time [seconds]'); ylabel('Amplitude');
% new_xticks = xticks/2000;
% xticklabels(arrayfun(@(x) num2str(x, '%.1f'), new_xticks, 'UniformOutput', false));
% box on; grid on
% plot([st st],[0 max(analytic_signal)+1],'-r')
% plot([fi fi],[0 max(analytic_signal)+1],'-r')
% plot([1 length(analytic_signal)],[m3dbAmp m3dbAmp],'--k')
% plot([t1 t1],[0 max(analytic_signal)+1],'-k')
% plot([t2 t2],[0 max(analytic_signal)+1],'-k')
% plot([1 length(analytic_signal)],[m10dbAmp m10dbAmp],'--k')%,'Color',[22 115 47]./256)
% plot([t3 t3],[0 max(analytic_signal)+1],'--k')%,'Color',[240 168 24]./256)
% plot([t4 t4],[0 max(analytic_signal)+1],'--k')%,'Color',[240 168 24]./256)
% axis tight
% % p2 = plot(nan, nan, '-r', 'DisplayName', 'Call Start/End');
% % p3 = plot(nan, nan, '--b', 'DisplayName', '-3dB Amplitude');
% % p4 = plot(nan, nan, '--m', 'DisplayName', '-3dB Start/End');
% % legend([p2, p3, p4], 'Call Start/End', '-3dB Amplitude', '-3dB Start/End');
% hold off
% 
% figure(13)
% subplot(3,1,1); hold on
% plot(data)
% xlabel('Time [seconds]'); ylabel('Counts');
% new_xticks = xticks/2000;
% xticklabels(arrayfun(@(x) num2str(x, '%.1f'), new_xticks, 'UniformOutput', false));
% box on; grid on
% plot([st st],[min(data) max(data)],'-r')
% plot([fi fi],[min(data) max(data)],'-r')
% title(sprintf('Waveform Bandpassed: [%i-%iHz]',parm.bp_min_f,parm.bp_max_f))
% axis tight
% p1 = plot(nan,nan,'-r');
% legend(p1,'Call Start/End')
% hold off
% 
% subplot(3,1,2); hold on
% plot(analytic_signal); xlim([0 length(data)]); title(sprintf('Amplitude Envelope'))
% xlabel('Time [seconds]'); ylabel('Amplitude');
% new_xticks = xticks/2000;
% xticklabels(arrayfun(@(x) num2str(x, '%.1f'), new_xticks, 'UniformOutput', false));
% box on; grid on
% plot([st st],[0 max(analytic_signal)+1],'-r')
% plot([fi fi],[0 max(analytic_signal)+1],'-r')
% plot([1 length(analytic_signal)],[m3dbAmp m3dbAmp],'--k')
% plot([t1 t1],[0 max(analytic_signal)+1],'-k')
% plot([t2 t2],[0 max(analytic_signal)+1],'-k')
% plot([1 length(analytic_signal)],[m10dbAmp m10dbAmp],'--k')%,'Color',[22 115 47]./256)
% plot([t3 t3],[0 max(analytic_signal)+1],'--k')%,'Color',[240 168 24]./256)
% plot([t4 t4],[0 max(analytic_signal)+1],'--k')%,'Color',[240 168 24]./256)
% axis tight
% hold off
% 
% subplot(3,1,3); hold on
% plot(csum)
% plot([e90_s e90_s],[0 1],'-k')
% plot([e90_f e90_f],[0 1],'-k')
% plot([0 e90_s],[0.05 0.05],'--k')
% plot([0 e90_f],[0.95 0.95],'--k')
% plot([e97_s e97_s],[0 1],'-k')
% plot([e97_f e97_f],[0 1],'-k')
% plot([0 e97_s],[0.015 0.015],'--k')
% plot([0 e97_f],[0.985 0.985],'--k')
% plot([st st],[0 1],'-r')
% plot([fi fi],[0 1],'-r')
% xlabel('Time [seconds]'); ylabel('Energy %');
% new_xticks = xticks/2000;
% xticklabels(arrayfun(@(x) num2str(x, '%.1f'), new_xticks, 'UniformOutput', false));
% title('Relative Cumulative Energy')
% box on; grid on
% axis tight
% hold off