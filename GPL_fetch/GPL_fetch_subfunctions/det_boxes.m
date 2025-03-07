function [rectHandle] = det_boxes(gpl_start, gpl_end, manual_start, manual_end, parm, fetched_call, cross_flag, confirmFlag,fig)
% This function will create boxes around GPL detections and a manual
% detection for a figure initialized in the primary script. Inputs:
% gpl_start: Start times (seconds) of GPL detections, gpl_end:
% Same for GPL end times, manual_start: Manual detection start time
% (seconds), manual_end: Manual detection end time (seconds), parm: GPL
% parameter structure, fetched_call: Manual detection structure 

flag = [];
rectHandle = [];

if confirmFlag == -1 % Add Adhoc and immediately end
    
    % Identify the number of subplots
    figure(fig);
    figgcf = gcf;
    axH = findall(figgcf,'Type','axes'); % Axis handles
    subplotH = axH(arrayfun(@(h) isempty(get(h,'Tag')),axH)); % Locate handles that may only be plots
    nSubplots = length(subplotH);
    
    % Add box
    if nSubplots == 3
        subplot(3,1,3); hold on
        h1 = plot3([manual_start manual_start],[parm.freq_lo-1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-m','LineWidth',1.5);
        h2 = plot3([manual_end manual_end],[parm.freq_lo-1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-m','LineWidth',1.5);
        h3 = plot3([manual_start manual_end],[parm.freq_lo-1 parm.freq_lo-1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-m','LineWidth',1.5);
        h4 = plot3([manual_start manual_end],[parm.freq_hi+1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-m','LineWidth',1.5);
        hold off
    elseif nSubplots == 2
        subplot(3,1,2); hold on
        h1 = plot3([manual_start manual_start],[parm.freq_lo-1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-m','LineWidth',1.5);
        h2 = plot3([manual_end manual_end],[parm.freq_lo-1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-m','LineWidth',1.5);
        h3 = plot3([manual_start manual_end],[parm.freq_lo-1 parm.freq_lo-1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-m','LineWidth',1.5);
        h4 = plot3([manual_start manual_end],[parm.freq_hi+1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-m','LineWidth',1.5);
        hold off
    end
    rectHandle = [h1 h2 h3 h4];
    return
end

if confirmFlag == -2 % TP acknowledge box
    
    % Identify the number of subplots
    figure(fig);
    figgcf = gcf;
    axH = findall(figgcf,'Type','axes'); % Axis handles
    subplotH = axH(arrayfun(@(h) isempty(get(h, 'Tag')), axH)); % Locate handles that may only be plots
    nSubplots = length(subplotH);
    
    % Add box
    if nSubplots == 3
        subplot(3,1,3); hold on
        plot3([manual_start manual_start],[parm.freq_lo-1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'--g','LineWidth',1.5);
        plot3([manual_end manual_end],[parm.freq_lo-1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'--g','LineWidth',1.5);
        plot3([manual_start manual_end],[parm.freq_lo-1 parm.freq_lo-1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'--g','LineWidth',1.5);
        plot3([manual_start manual_end],[parm.freq_hi+1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'--g','LineWidth',1.5);
        hold off
    elseif nSubplots == 2
        subplot(3,1,2); hold on
        plot3([manual_start manual_start],[parm.freq_lo-1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'--g','LineWidth',1.5);
        plot3([manual_end manual_end],[parm.freq_lo-1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'--g','LineWidth',1.5);
        plot3([manual_start manual_end],[parm.freq_lo-1 parm.freq_lo-1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'--g','LineWidth',1.5);
        plot3([manual_start manual_end],[parm.freq_hi+1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'--g','LineWidth',1.5);
        hold off
    end
    return
end

if confirmFlag == -3 % Adhoc reprocess dashed box
    
    % Identify the number of subplots
    figure(fig);
    figgcf = gcf;
    axH = findall(figgcf,'Type','axes'); % Axis handles
    subplotH = axH(arrayfun(@(h) isempty(get(h, 'Tag')), axH)); % Locate handles that may only be plots
    nSubplots = length(subplotH);
    
    % Add box
    if nSubplots == 3
        subplot(3,1,3); hold on
        plot3([manual_start manual_start],[parm.freq_lo-1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'--m','LineWidth',2);
        plot3([manual_end manual_end],[parm.freq_lo-1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'--m','LineWidth',2);
        plot3([manual_start manual_end],[parm.freq_lo-1 parm.freq_lo-1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'--m','LineWidth',2);
        plot3([manual_start manual_end],[parm.freq_hi+1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'--m','LineWidth',2);
        hold off
    elseif nSubplots == 2
        subplot(3,1,2); hold on
        plot3([manual_start manual_start],[parm.freq_lo-1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'--m','LineWidth',2);
        plot3([manual_end manual_end],[parm.freq_lo-1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'--m','LineWidth',2);
        plot3([manual_start manual_end],[parm.freq_lo-1 parm.freq_lo-1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'--m','LineWidth',2);
        plot3([manual_start manual_end],[parm.freq_hi+1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'--m','LineWidth',2);
        hold off
    end
    return
end

% GPL boxes
if ~isempty(gpl_start)
    for n = 1:length(gpl_start) % plot GPL detections
        plot3([gpl_start(n) gpl_start(n)],[parm.freq_lo-1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-r');
        plot3([gpl_end(n) gpl_end(n)],[parm.freq_lo-1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-r');
        plot3([gpl_start(n) gpl_end(n)],[parm.freq_lo-1 parm.freq_lo-1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-r');
        plot3([gpl_start(n) gpl_end(n)],[parm.freq_hi+1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-r');
        flag = 1;
    end
end

% Manual Box
if confirmFlag == 0
    % if cross_flag == 0
        plot3([manual_start manual_start],[parm.freq_lo-1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-g','LineWidth',1.5);
        plot3([manual_end manual_end],[parm.freq_lo-1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-g','LineWidth',1.5);
        plot3([manual_start manual_end],[parm.freq_lo-1 parm.freq_lo-1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-g','LineWidth',1.5);
        plot3([manual_start manual_end],[parm.freq_hi+1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-g','LineWidth',1.5);
        flag = [flag 2];
    % else
        % plot3([manual_start-parm.window_length manual_start-parm.window_length],[parm.freq_lo-1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-g','LineWidth',1.5);
        % plot3([manual_end-parm.window_length manual_end-parm.window_length],[parm.freq_lo-1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-g','LineWidth',1.5);
        % plot3([manual_start-parm.window_length manual_end-parm.window_length],[parm.freq_lo-1 parm.freq_lo-1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-g','LineWidth',1.5);
        % plot3([manual_start-parm.window_length manual_end-parm.window_length],[parm.freq_hi+1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-g','LineWidth',1.5);
        % flag = [flag 2];
    % end
else
    plot3([manual_start manual_start],[parm.freq_lo-1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-y','LineWidth',3);
    plot3([manual_end manual_end],[parm.freq_lo-1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-y','LineWidth',3);
    plot3([manual_start manual_end],[parm.freq_lo-1 parm.freq_lo-1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-y','LineWidth',3);
    plot3([manual_start manual_end],[parm.freq_hi+1 parm.freq_hi+1],[max(max(fetched_call.spectrogram)) max(max(fetched_call.spectrogram))],'-y','LineWidth',3);
end

gpl_dummy = plot(nan, nan, '-r');
manual_dummy = plot(nan, nan, '-g');

% Legend
if (parm.fp1.legend == 1) && (confirmFlag == 0)
    if length(flag) == 2 % GPL & Manual
        legend([gpl_dummy manual_dummy],'GPL Detections','Manual Detection','Location','southoutside','Orientation','horizontal')
    else
        if flag == 1 % GPL only
            legend(gpl_dummy,'GPL Detections','Location','southoutside','Orientation','horizontal')
        elseif flag == 2 % Manual only
            legend(manual_dummy,'Manual Detection','Location','southoutside','Orientation','horizontal')
        else % No plot
            return
        end
    end
end