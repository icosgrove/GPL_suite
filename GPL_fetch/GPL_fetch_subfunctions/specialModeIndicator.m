function specialModeIndicator(string1,type,state)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will add text below the third subplot indicating which
% special case mode is currently chosen. The previous indicator will be
% removed, if applicable.
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

position1 = [0.35, -0.22, 0]; % pairing
position2 = [0.35, -0.40, 0]; % window

% Fix empty spaces 
if nargin == 3
    if ~isempty(state{3})
        state{3} = state{3}(~cellfun('isempty',state{3}));
    end
end

% Delete previous text
boxes1 = findobj(gca,'Type','text','Position',position1);
if ~isempty(boxes1) % clear previous pairing indicator
    if type == 1 || type == 3
        delete(boxes1)
    end
end
boxes2 = findobj(gca,'Type','text','Position',position2);
if ~isempty(boxes2) % clear previous window indicator
    if type == 2
        delete(boxes2)
    end
end

% Add window notes to each other
if nargin == 3
    for n = 1:length(state{3})
        if isempty(state{3}(n))
            continue
        else
            switch string(state{3}(n))
                case 'fn'
                    string1 = strcat(string1,' + False Negative');
                case 'x'
                    string1 = strcat(string1, ' + User Error');
            end
        end
    end
end

% Add new text in plot
if type == 1 % Pairing note
    text('Units', 'normalized', 'Position', position1, ...
            'String', string1, 'HorizontalAlignment', 'right','BackgroundColor',[252 224 144]./255);
elseif type == 2 % Window Indicator
    text('Units', 'normalized', 'Position', position2, ...
            'String', string1, 'HorizontalAlignment', 'right','BackgroundColor',[252 224 144]./255);
elseif type == 3 % Adhoc
    text('Units', 'normalized', 'Position', position1, ...
            'String', string1, 'HorizontalAlignment', 'right','BackgroundColor',[241 150 242]./255);
end