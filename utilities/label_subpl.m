function [hsh,vsh] = label_subpl(sh,xlabels,ylabels)
%
% label_subpl(sh,xlabels,ylabels)
%
% label the bottom subplots with the xlabels, 
% label the left   subplots with the ylabels,
% no labels to the inner subplots, if any, they are cleared
%
% usage   : [hsh,vsh] = label_subpl(sh,xlabels,ylabels);
% example : [hsh,vsh] = label_subpl(sh,{'time [ms]'},{'as' 'time' 'goes' 'by'});
%
% INPUT ARGUMENTS
%
% 1. sh       vector or matrix of the subhandles to the subplots,
%             e.g., returned from gensubpl.m
% 2. xlabels  cell array of xlabels, if xlabels contains only one label,
%             this label will be assigned to all bottom subplots
% 3. ylabels  cell array of ylabels, if ylabels contains only one label,
%             this label will be assigned to all left subplots
%
% OUTPUT ARGUMENTS
%
% 1. hsh      horizontal subhandles at the bottom (with xlabels)
% 2. vsh      vertical   subhandles at the bottom (with ylabels)
%
% AUTHOR: Helmut Riedel, April 2000

if ~nargin help label_subpl; return; end 
if ~iscell(xlabels) | ~iscell(ylabels)
  error('label_subpl: xlabels and ylabels have to be cell arrays');
end

for hv = 1:length(sh(:))
  pos(hv,:) = get(subplot(sh(hv)),'position');     % pos = [hmin vmin width height]
  xlabel(''); ylabel('');                          % clear all existing labels
end %hv

hsh = sh(find(min(pos(:,2)) == pos(:,2)));         % subplot handles with minimal x
vsh = sh(find(min(pos(:,1)) == pos(:,1)));         % subplot handles with minimal y

H = length(hsh); V = length(vsh);
if H*V ~= prod(size(sh))
  error('label_subpl: number of subplots does not match length(sh), exiting');
end

if length(xlabels) == 1                            % only one xlabel provided
  xlabels = repmat(xlabels,[H 1]);                 % multiplex it for H subplots
end                                            
if length(ylabels) == 1                            % only one ylabel provided
  ylabels = repmat(ylabels,[V 1]);                 % multiplex it for V subplots
end                                            
                                               
if length(xlabels) < H                             % less than H xlabels provided ?
  xlabels{H} = '';                                 % set remaining xlabels to empty
end                                            
if length(ylabels) < V                             % less than V ylabels provided ?
  ylabels{V} = '';                                 % set remaining ylabels to empty
end

for h=1:H
  subplot(hsh(h)); 
  xlabel(xlabels{h},'VerticalAlignment','top','FontSize',14);    % set the xlabels,
end                                                % let some distance to the xticklabels
for v=1:V                                      
  subplot(vsh(v));
  ylabel(ylabels{v},'VerticalAlignment','bottom','FontSize',14); % set the ylabels,
end                                                % let some distance to the yticklabels
return
