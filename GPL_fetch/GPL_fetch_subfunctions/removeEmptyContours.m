function notEmpty = removeEmptyContours(c1,c2,parm,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will remove empty contours so they are not combined witha
% full one
% Written: Ian 07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if ~isempty(c1.index) 
%     if ~isempty(c2.index)
%         notEmpty = combineContours(c1,c2,parm,n); % Both nonz
%     else
%         notEmpty = [c1 zeros(c2.size(1),c2.size(2))]; % nonz: c1 only
%     end
% else
%     if ~isempty(c2.index) % nonz: c2 only
%         notEmpty = c2;
%     else
%         notEmpty = []; % Both fail
%     end
% end
notEmpty = combineContours(c1,c2,parm,n);