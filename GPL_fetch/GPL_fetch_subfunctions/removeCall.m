function [fetched_call] = removeCall(fetched_call)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will remove the current paired GPL detection, if there is
% one to remove. Intended for use when pairings are being overwritten.
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(fetched_call.gpl_match)
    fprintf('Previous Pairing and note is being overwritten.\n')
    fetched_call.gpl_match = [];
    fetched_call.note = '';
end