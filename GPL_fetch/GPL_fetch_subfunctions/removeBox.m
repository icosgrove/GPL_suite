function removeBox(fig,missI,rectHandle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will remove the yellow detection box that delimits the
% current pairing in GPL fetch process 1. This is done so that the user may
% visually confirm that the previous pairing selection has been
% overwritten.
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(fig); subplot(3,1,3); hold on

% Remove grey undetermined box
gbox = findobj(gca,'Type','line','Color',[179 179 179]./255);
if ~isempty(gbox)
    delete(gbox)
end

% Remove green box for Tp ack case
if nargin == 2
    if missI == 1
        mline = findobj(gca,'Type','Line','Color','g');
        delete(mline)
        return
    end
end

% Remove solid adhoc box for replacing
if nargin == 3
    if missI == 2
        delete(rectHandle)
        % aline = findobj(gca,'Type','Line','Color','m');
        % delete(aline)
        return
    end
end

% Locate yellow lines and remove them
ylines = findobj(gca,'Type','line','Color','y');
if ~isempty(ylines)
    delete(ylines)
end

% Remove switch indicator
missbox = findobj(gca,'Type','line','Color',[245 138 56]./255);
if ~isempty(missbox)
    if nargin == 2 % Skip removal for switch or related cases
        return
    else
        delete(missbox)
    end
end

% Remove miss indicator
switchI = findobj(gca,'Type','line','Color',[4 204 197]./255);
if ~isempty(switchI)
    if nargin == 2
        return
    else
        delete(switchI)
    end
end

% Remove the false flag
falseRem = findobj(gca,'Type','line','LineWidth',4);
if ~isempty(falseRem)
    delete(falseRem)
end

hold off