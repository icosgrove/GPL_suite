function [STATS] = unloadStats(plot_calls)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will extract measurements from GPL_fetch outputs into a new
% struct
% Written: Ian 09/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:length(plot_calls)
    STATS.cm.dur_s(k) = plot_calls(k).cm.duration_sec;
    STATS.cm.slope(k) = plot_calls(k).cm.slope;
    STATS.cm.band_hz(k) = plot_calls(k).cm.freq_bandwidth_hz;
    STATS.cm.start_hz(k) = plot_calls(k).cm.start_freq_hz;
    STATS.cm.end_hz(k) = plot_calls(k).cm.end_freq_hz;
    STATS.cm.min_hz(k) = plot_calls(k).cm.min_freq_hz;
    STATS.cm.max_hz(k) = plot_calls(k).cm.max_freq_hz - 2 ;
    STATS.cm.abs_band_hz(k) = plot_calls(k).cm.abs_bandwidth_hz;
    STATS.cm.peak_freq(k) = plot_calls(k).cm.peak_freq;
    STATS.cm.p2p_rl(k) = plot_calls(k).cm.p2p_RL;
    STATS.cm.rms_rl_m3(k) = plot_calls(k).cm.rms_RL.m3db;
    STATS.cm.rms_rl_m10(k) = plot_calls(k).cm.rms_RL.m10db;
    STATS.cm.rms_rl_e90(k) = plot_calls(k).cm.rms_RL.e90;
    STATS.cm.rms_rl_e97(k) = plot_calls(k).cm.rms_RL.e97;
    STATS.cm.rms_rl_gplSE(k) = plot_calls(k).cm.rms_RL.gpl_SE;
    STATS.cm.SEL_m3(k) = plot_calls(k).cm.SEL.m3db;
    STATS.cm.SEL_m10(k) = plot_calls(k).cm.SEL.m10db;
    STATS.cm.SEL_e90(k) = plot_calls(k).cm.SEL.e90;
    STATS.cm.SEL_e97(k) = plot_calls(k).cm.SEL.e97;
    STATS.cm.SEL_gplSE(k) = plot_calls(k).cm.SEL.gpl_SE;
end

for k = 1:length(plot_calls)
    STATS.cm_max.dur_s(k) = plot_calls(k).cm_max.duration_sec;
    STATS.cm_max.slope(k) = plot_calls(k).cm_max.slope;
    STATS.cm_max.band_hz(k) = plot_calls(k).cm_max.freq_bandwidth_hz;
    STATS.cm_max.start_hz(k) = plot_calls(k).cm_max.start_freq_hz;
    STATS.cm_max.end_hz(k) = plot_calls(k).cm_max.end_freq_hz;
    STATS.cm_max.min_hz(k) = plot_calls(k).cm_max.min_freq_hz;
    STATS.cm_max.max_hz(k) = plot_calls(k).cm_max.max_freq_hz;
    STATS.cm_max.abs_band_hz(k) = plot_calls(k).cm_max.abs_bandwidth_hz;
    STATS.cm_max.peak_freq(k) = plot_calls(k).cm_max.peak_freq;
    STATS.cm_max.p2p_rl(k) = plot_calls(k).cm_max.p2p_RL;
    STATS.cm_max.rms_rl_m3(k) = plot_calls(k).cm_max.rms_RL.m3db;
    STATS.cm_max.rms_rl_m10(k) = plot_calls(k).cm_max.rms_RL.m10db;
    STATS.cm_max.rms_rl_e90(k) = plot_calls(k).cm_max.rms_RL.e90;
    STATS.cm_max.rms_rl_e97(k) = plot_calls(k).cm_max.rms_RL.e97;
    STATS.cm_max.rms_rl_gplSE(k) = plot_calls(k).cm_max.rms_RL.gpl_SE;
    STATS.cm_max.SEL_m3(k) = plot_calls(k).cm_max.SEL.m3db;
    STATS.cm_max.SEL_m10(k) = plot_calls(k).cm_max.SEL.m10db;
    STATS.cm_max.SEL_e90(k) = plot_calls(k).cm_max.SEL.e90;
    STATS.cm_max.SEL_e97(k) = plot_calls(k).cm_max.SEL.e97;
    STATS.cm_max.SEL_gplSE(k) = plot_calls(k).cm_max.SEL.gpl_SE;
end