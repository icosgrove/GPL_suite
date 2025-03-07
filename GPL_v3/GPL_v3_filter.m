% GPL_v3_filter
% This script will apply measurement fitlers to GPL_v3 detection output
% files.
% Ranges for filtering should be identified before running this filter, for
% example, a slope range of -10 to -5 Hz/s should be identified using GPL
% fetch measurement statistics or some other method, then manually input
% here for filtering. 
% Any combination of filters may be set.
% Written: Ian 10/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% USER SETUP

% Combine filtered detections into one file?
setup.combine = 0;

% Plot measurement statistics of raw detections before filtering?
setup.prePlot = 0;


%% Processing

% Load detection files:
disp('Select folder containing GPL_v3 detection files.')
[filename, pathname] = uigetfile('*.mat'); % Select File
cwd = pwd; 
cd(pathname) % Set current directory to path containing xwav file
file_dir = pwd; 
addpath(pwd); 
files = dir('*.mat'); % Variable files contains xwav data
cd(cwd); % Set current directory back to current working directory