% Code to separate calls by a species specification in the Comments or
% Parameter6 column of the detlog.
% Written by Ian 08/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User Input

% Output filename
outfnam = 'NUNAT_SB_03_TritonLog_KH_up.xlsx';

% Enter call type and filter column
opts.call_type = 'up'; % Expected entry in the column
opts.comments = 1; % 1:on, 0:off
opts.parameter6 = 1; % 1:on, 0:off


%% Analysis
% Pick files
[filename, pathname] = uigetfile('*.xlsx'); % Select File
cwd = pwd; 
cd(pathname) % Set current directory to path containing xwav file
file_dir = pwd; 
addpath(pwd); 
files = dir('*.xlsx'); % Variable files contains xwav data
cd(cwd); % Set current directory back to current working directory

% Allocate
cData = [];

% switch to cell if needed
if ischar(files)
    files = {files};
end

% Loop through each selected file
for n = 1:length(files)
    
    % read
    fname = files(n).name;
    opts = detectImportOptions(fname, 'PreserveVariableNames', true);
    opts = setvartype(opts, 'Parameter 6', 'char');
    data = readtable(fname,opts);
    data.Properties.VariableNames{'Parameter 6'} = 'Parameter6';
    data.Properties.VariableNames{'Start time'} = 'StartTime';
    data.Properties.VariableNames{'End time'} = 'EndTime';

    % find species type and append
    if opts.comments == 1
        filteredData = data(contains(data.Comments,opts.call_type),:);
    end

    % Parameter6 labels
    if opts.parameter6 == 1
        data_copy = data;
        data_copy = data_copy(contains(data.Parameter6,opts.call_type),:);
        filteredData1 = data(contains(data.Parameter6,opts.call_type),:);
        filteredData = [filteredData; filteredData1];
        filteredData = unique(filteredData);
    end

    cData = [cData; filteredData];

    clear filteredData
end

% switch to julian time & sort
cData.StartTime = datenum(datetime(cData.StartTime, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.FFF'));
cData.EndTime = datenum(datetime(cData.EndTime, 'InputFormat', 'yyyy-MM-dd HH:mm:ss.FFF'));
sortedData = sortrows(cData, 'StartTime');
sortedData(end-1,:) = [];

% save

writetable(sortedData, outfnam);