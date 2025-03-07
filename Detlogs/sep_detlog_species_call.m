%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script creates separate excel detlogs for each unique call type in 
% an exising detlog. The user may specify to separate by call type or by
% species code in user setup. 
% Multiple workbooks may be loaded as long as they are in the same folder. 
%
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User Input

% Output name for filtered detlog
outputfname_front = 'NUNAT_SB_03_TritonLog_LowFreq_KH';

% Separate by 'call' or 'species'?
sep_type = 'species'; % Enter 'call' or 'species'

warning('off')



%% Analysis
% Load folder containing .xlsx files
fprintf('\nSelect folder containing Excel .xlsx to process\n');
[filename, pathname] = uigetfile('*.xlsx'); % Select File
cwd = pwd; 
cd(pathname) % Set current directory to path containing xwav file
file_dir = pwd; 
addpath(pwd); 
files = dir('*.xlsx'); % Variable files contains xwav data
cd(cwd); % Set current directory back to current working directory

cdata = table();

% Append workbooks together
for n = 1:length(files)
    tdata = readtable(files(n).name); 
    
    % Convert to julian time to preserve decimals
    tdata.StartTime = datenum(tdata.StartTime);
    tdata.EndTime = datenum(tdata.EndTime);
    if ismember('Comments',tdata.Properties.VariableNames) % Handle comments not existing
        tdata.Comments = []; % They are fully removed from output 
    end
    cdata = [cdata; tdata];
end

% Create excel files for each call type
switch sep_type
    case 'call'
        callTypes = unique(cdata.Call);
        for k = 1:length(callTypes)
            sepData = cdata(strcmp(cdata.Call,callTypes{k}),:);
            outputName = strcat(outputfname_front,'_',callTypes{k},'.xlsx');
            writetable(sepData,outputName);  
        end
    case 'species'
        sepciesTypes = unique(cdata.SpeciesCode);
        for k = 1:length(sepciesTypes)
            sepData = cdata(strcmp(cdata.SpeciesCode,sepciesTypes{k}),:);
            outputName = strcat(outputfname_front,'_',sepciesTypes{k},'.xlsx');
            writetable(sepData,outputName);  
        end
end