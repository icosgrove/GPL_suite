function [plot_calls, adhoc_off, reprocess_calls] = bin1Calls(calls,adhoc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will separate GPL Fetch output calls into 3 bins: calls
% ready to got statistical analysis, off-effort adhocs to be handled
% differently, and calls requiring reprocessing.
% Written: Ian 09/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a pile of calls with 1 note and 1+ notes
single_note = []; reprocess_calls = [];
for n = 1:length(calls)
    if length(strcmp(calls(n).note,'')) == 1
        single_note = [single_note n];
    else
        reprocess_calls = [reprocess_calls n];
    end
end
single_note = calls(single_note);
reprocess_calls = calls(reprocess_calls);

% Pull out singles & swaps
easy_calls = []; 
for n = 1:length(single_note)
    if strcmp(single_note(n).note,'Single Match')
        easy_calls = [easy_calls n];
    elseif strcmp(single_note(n).note,'Detection Swap')
        easy_calls = [easy_calls n];
    else
        reprocess_calls = [reprocess_calls single_note(n)];
    end
end
easy_calls = single_note(easy_calls);

% Pull out on effort adhocs 
adhoc_ready = []; adhoc_off = [];
for n = 1:length(adhoc)
    if strcmp(adhoc(n).note,'On Effort')
        adhoc_ready = [adhoc_ready n];
    else
        adhoc_off = [adhoc_off n];
    end
end
adhoc_ready = adhoc(adhoc_ready);
adhoc_off = adhoc(adhoc_off);

% Combine
for n = 1:length([easy_calls.gpl_match])
    easy_calls(n).gpl_match.note = [];
end
easy = [easy_calls.gpl_match];
easy = easy(:);
plot_calls = [easy; adhoc_ready];