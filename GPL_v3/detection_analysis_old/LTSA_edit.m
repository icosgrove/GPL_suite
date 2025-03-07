%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script will add plotted elements to the Triton LTSA. It is assumed
% that the LTSA is plotted in Figure 1. Intended use is for analyzing
% detection frequency plots from the GPL detector. Access individual code
% blocks to change appearance parameters (External toggles not currently
% provided).
%
% Written by Ian Cosgrove 05/26/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in GPL output detection file
curr_folder = pwd;
[fname, pname] = uigetfile('*.mat', 'Select .mat file', curr_folder);

% Load detection or display selection cancellation error
if isequal(fname, 0)
    disp('File selection cancelled');
else
    fname_full = fullfile(pname, fname);
    load(fname_full); % Load
end


%% User Input

% LTSA Manual Inputs
x_start = 0; x_end = 179; % Hour number for LTSA start/end

% Select what you want added to LTSA
freqLimits = 1; % GPL detector freq. min and max
detFreq_subplot = 1; % Add detection frequency plot as a subplot below the LTSA

f1 = figure(1); hold on

% Plot Frequency Limits
if freqLimits == 1
    plot([x_start x_end], hyd.detection.parm.freq_lo*ones(1,2), '-r', 'LineWidth',2)
    plot([x_start x_end], hyd.detection.parm.freq_hi*ones(1,2), '-r', 'LineWidth',2)
end


hold off