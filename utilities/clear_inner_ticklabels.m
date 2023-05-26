function clear_inner_ticklabels(sh)
%
% clear_inner_ticklabels -- clears the inner ticklabels in a HxV-subplot
%
% usage   : clear_inner_ticklabels(sh);
% example : clear_inner_ticklabels(sh);
%
% INPUT ARGUMENTS
%
% 1. sh     vector or matrix of the subhandles to the subplots,
%           e.g., returned from gensubpl.m
%
% OUTPUT ARGUMENTS
%
% none
%
% AUTHOR: Helmut Riedel, April 2000

if ~nargin help clear_inner_ticklabels; return; end 

for hv = 1:length(sh(:))
  pos(hv,:) = get(subplot(sh(hv)),'position'); % position = [hmin vmin width height]
end %hv

hsh = sh(find(min(pos(:,2)) ~= pos(:,2)));     % subplots with nonminimal y-position
vsh = sh(find(min(pos(:,1)) ~= pos(:,1)));     % subplots with nonminimal x-position

for h = 1:length(hsh)
  set(hsh(h),'xticklabel','');
end
for v = 1:length(vsh)
  set(vsh(v),'yticklabel','');
end
return
