function [TABLE,cm_table,cm_max_table] = makeTable(STATS,tablePARAMS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the mean, std. dev., and range of values for each requested
% measurement
% Written: Ian 09/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All metrics
all_metrics = {'Duration','Slope','Nonzero Frequency Bandwidth','Start Frequency','End Frequency',...
    'Minimum Frequency','Maximum Frequency','Absolute Frequency Bandwidth','Peak Frequency',...
    'Peak-to-Peak Received Level','RMS Received Level [-3dB]','RMS Received Level [-10dB]',...
    'RMS Received Level [90% Energy]','RMS Received Level [97% Energy]','RMS Received Level [GPL]',...
    'Sound Exposure Level [-3dB]','Sound Exposure Level [-10dB]','Sound Exposure Level [90% Energy]',...
    'Sound Exposure Level [97% Energy]','Sound Exposure Level [GPL]'};

% Find selected metrics
cmfields = fields(STATS.cm);
k1 = 1;
for n1 = 1:length(cmfields)
    if tablePARAMS.cm.(cmfields{n1}) == 1
        chosen_metrics_cm{k1,1} = all_metrics{n1};
        k1=k1+1;
    end
end
cm_maxfields = fields(STATS.cm_max);
k2 = 1;
for n2 = 1:length(cm_maxfields)
    if tablePARAMS.cm_max.(cm_maxfields{n2}) == 1
        chosen_metrics_cm_max{k2,1} = all_metrics{n2};
        k2=k2+1;
    end
end

% Create cm table
m = 1;
for n = 1:length(cmfields)
    if tablePARAMS.cm.(cmfields{n}) == 1
        stat = STATS.cm.(cmfields{n});
        TABLE.cm.metric{m} = chosen_metrics_cm{m};
        TABLE.cm.mean(m) = mean(stat);
        TABLE.cm.stdDev(m) = std(stat);
        TABLE.cm.min(m) = min(stat);
        TABLE.cm.max(m) = max(stat);
        % TABLE.cm.(cmfields{n}) = [mean(stat) std(stat) min(stat) max(stat)];
        m=m+1;
    end
end
cm_table = table(TABLE.cm.metric',TABLE.cm.mean',TABLE.cm.stdDev',TABLE.cm.min',TABLE.cm.max',...
    'VariableNames',{'Metric','Mean','StdDev','Min','Max'});

% cm_max table
m = 1;
for n = 1:length(cm_maxfields)
    if tablePARAMS.cm_max.(cm_maxfields{n}) == 1
        stat = STATS.cm_max.(cm_maxfields{n});
        TABLE.cm_max.metric{m} = chosen_metrics_cm_max{m};
        TABLE.cm_max.mean(m) = mean(stat);
        TABLE.cm_max.stdDev(m) = std(stat);
        TABLE.cm_max.min(m) = min(stat);
        TABLE.cm_max.max(m) = max(stat);
        % TABLE.cm_max.(cm_maxfields{n}) = [mean(stat) std(stat) min(stat) max(stat)];
        m=m+1;
    end
end
cm_max_table = table(TABLE.cm_max.metric',TABLE.cm_max.mean',TABLE.cm_max.stdDev',TABLE.cm_max.min',TABLE.cm_max.max',...
    'VariableNames',{'Metric','Mean','StdDev','Min','Max'});