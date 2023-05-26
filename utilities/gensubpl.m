function [sh,fh] = gensubpl(mode,H,V,hf,vf,hin,vin,hout,vout,htl,vtl)
%
% gensubpl.m - generates a subplot with or without outer axis,
%              returns a vector or matrix of subhandles sh to the individual subplots
%
% usage   : [sh,fh] = gensubpl(mode,H,V,hf,vf,hin,vin,hout,vout,htl,vtl);
% examples: [sh,fh] = gensubpl('hblo',3,4,0.98,0.98,'t [ms]','U [uV]','ITD [ms]', ...
%                         'ILD [dB]',{'-0.3' '0.0' '0.3'},{'2' '4' '6' '8'});
%           [sh,fh] = gensubpl('vtlm',5,3,i,i,'time [s]','voltage [nV]');
%           [sh,fh] = gensubpl('vbln',4,5);
%
% INPUT ARGUMENTS (input arguments not given or set to i are replaced by the defaults marked by '*')
%
% 1. mode   'h' or 'v' : count quicker in horizontal    (*) or vertical direction  
%           't' or 'b' : first subplot is at the bottom (*) or top row
%           'l' or 'r' : first subplot is at the left   (*) or right column
%           The 8 possibilities for the example H=3 and V=2 are:
%        
%         * 'hbl': 4 5 6    'hbr': 6 5 4    'htl': 1 2 3    'htr': 3 2 1    
%                  1 2 3           3 2 1           4 5 6           6 5 4    
%                                                                           
%           'vbl': 2 4 6    'vbr': 6 4 2    'vtl': 1 3 5    'vtr': 5 3 1    
%                  1 3 5           5 3 1           2 4 6           6 4 2    
%
%           'c' clear the ticklabels of the inner plots, i.e., only subplots with axis at
%               the left or bottom side of the whole plot maintain their ticklabels.
%               This is useful for densely placed subplots.
%               The function called is clear_inner_ticklabels(sh), it will be usually
%               called AFTER plotting some data into the subplots since this
%               usually generates some xticklabels.
%         * 'f' open a new figure, without 'f' an existing figure will be used, if there is one
%         * 'g' grid on in every subplot, (hold and zoom is on automatically)
%           'm' return a subplot-handle-matrix rather than a vector.
%               CAUTION: in v-mode sh is a normal     matrix Nrows x Ncolumns
%                        in h-mode sh is a transposed matrix Ncolums x Nrows
%         * 'n' enumarate the subplots to see how they are arranged, a test
%           'o' outer axis, only in this case hin, vin, hout, vout, htl, vtl
%               are used. This mode (introducing an outer axis) somehow disables 
%               zooming in the subplots. I don't know how to avoid that.
% 2. H      number of horizontal subplots (= number of columns), default 3
% 3. V      number of vertical   subplots (= number of rows   ), default 3
% 4. hf     horizontal visible fraction, 0 <= hf <= 1, default 0.98 
% 5. vf     vertical   visible fraction, 0 <= vf <= 1, default 0.98
% 6. hin    string with xlabel of the inner horizontal axis, default 'time [ms]'
% 7. vin    string with ylabel of the inner vertical   axis, default 'voltage [\muV]'
% 8. hout   string with xlabel of the outer horizontal axis, default 'ITD [ms]'
% 9. vout   string with ylabel of the outer vertical   axis, default 'ILD [dB]'
%10. htl    cell array (row) with xticklabels on the horizontal axis, only for the lower row
%           in 'r'-mode htl will be flipped left-right, default {'-0.4' '0.0' '0.4'}
%11. vtl    cell array (row) with yticklabels on the vertical axis, only for the left column
%           in 't'-mode vtl will be flipped up-down, default {'-12' '0' '12'}
%
% OUTPUT ARGUMENTS
%
% 1. sh     the subplot handles needed to plot something into the subplots
%           sh can be returned as 
%           - vector (default) , then, e.g., write subplot(sh(n))  ; plot(x);
%           - matrix (mode 'm'), then, e.g., write subplot(sh(j,k)); plot(x);
%             in h-mode, the first index runs horizontal (along the columns)
%             in v-mode, the first index runs vertical   (along the rows)
% 2. fh     the figure handle of a new figure, e.g., needed by hgsave
%
% AUTHOR: Helmut Riedel, October 2000

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. default handling: set to default if not given or set to i %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 11    vtl  = {'-12' '0' '12'}    ; end
if nargin < 10    htl  = {'-0.4' '0.0' '0.4'}; end
if nargin <  9    vout = 'ILD [dB]'          ; end
if nargin <  8    hout = 'ITD [ms]'          ; end
if nargin <  7    vin  = 'voltage [\muV]'    ; end
if nargin <  6    hin  = 'time [ms]'         ; end
if nargin <  5    vf   = 0.98                ; end
if nargin <  4    hf   = 0.98                ; end
if nargin <  3    V    = 3                   ; end
if nargin <  2    H    = 3                   ; end
if nargin <  1    mode = 'hblfgn'            ; end

if isequal(vtl,i) vtl  = {'-12' '0' '12'}    ; end
if isequal(htl,i) htl  = {'-0.4' '0.0' '0.4'}; end
if vout   == i    vout = 'ILD [dB]'          ; end 
if hout   == i    hout = 'ITD [ms]'          ; end 
if vin    == i    vin  = 'voltage [\muV]'    ; end 
if hin    == i    hin  = 'time [ms]'         ; end 
if vf     == i    vf   = 0.98                ; end 
if hf     == i    hf   = 0.98                ; end 
if V      == i    V    = 3                   ; end 
if H      == i    H    = 3                   ; end 
if isempty(mode) | mode == i
  mode = 'x';                                         % avoid stupid warnings
end

%%
%% alternative handling of input arguments with varargin:
%% 'disadvantage': for default arguments [] has to be given instead of i
%%
%% function line:
% function [sh,fh] = gensubpl(varargin)
%argin = get_argin(varargin,'hblfgn',3,3,0.98,0.98,'time [ms]','voltage [\muV]', ...
%                        'ITD [ms]','ILD [dB]',{'-0.4' '0.0' '0.4'},{'-12' '0' '12'});
%
%mode = argin{1}; 
%H    = argin{2};  V    = argin{3};
%hf   = argin{4};  vf   = argin{5};
%hin  = argin{6};  vin  = argin{7};
%hout = argin{8};  vout = argin{9};
%htl  = argin{10}; vtl  = argin{11};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. check consistency of mode and set mode defaults, some input argument tests %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if find(mode=='h') & find(mode=='v')
  error('gensubpl: mode cannot contain h(orizontal) and v(ertical).');
end
if find(mode=='t') & find(mode=='b')
  error('gensubpl: mode cannot contain t(op) and b(ottom).');
end
if find(mode=='l') & find(mode=='r')
  error('gensubpl: mode cannot contain l(eft) and r(ight).');
end

if isempty(find(mode=='h')) & isempty(find(mode=='v'))  % no 'h' or 'v' ?
  mode = [mode 'h'];                                    % set default horizontal
end
if isempty(find(mode=='b')) & isempty(find(mode=='t'))  % no 'b' or 't' ?
  mode = [mode 'b'];                                    % set default bottom
end
if isempty(find(mode=='l')) & isempty(find(mode=='r'))  % no 'l' or 'r' ?
  mode = [mode 'l'];                                    % set default left
end

if size(htl,1) > 1                 % if htl is a colomn
  htl = htl';                      % convert to a row, because fliplr is used
end
if size(vtl,1) > 1                 % if vtl is a colomn
  vtl = vtl';                      % convert to a row, because fliplr is used
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. define positions using the default mode 'hbl': horizontal, bottom and left  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HV = H*V;                          % total number of subplots
sh = zeros(HV,1);                  % vector with the handles of the subplots

if find(mode=='o')                 % some space is needed for the outer axis
  h0 = 0.19;                         % global horizontal origin
  v0 = 0.19;                         % global vertical   origin
  hw = 0.74/H;  % 0.76/H;            % global horizontal width
  vw = 0.76/V;  % 0.76/V;            % global vertical   width
else                               % bigger subplots possible, no space for outer axis
  h0 = 0.12;                         % global horizontal origin
  v0 = 0.12;                         % global vertical   origin
  hw = 0.81/H;  % 0.83/H;            % global horizontal width
  vw = 0.83/V;  % 0.83/V;            % global vertical   width
end %'o'

if H==1
  h = h0 + hw*  rem((0:HV-1),H);                      % horizontal origins of subplots
else
  h = h0 + hw*  rem((0:HV-1),H) * (1 + (1-hf)/(H-1)); % horizontal origins of subplots
end

if V==1
  v = v0 + vw*floor((0:HV-1)/H);                      % vertical   origins of subplots
else
  v = v0 + vw*floor((0:HV-1)/H) * (1 + (1-vf)/(V-1)); % vertical   origins of subplots
end

% note: for hf=vf=1 the horizontal and vertical origins of the subplots are:
%
% h = h0 + hw*[0 1 2 .. H  ....  0 1 2   .. H   ....   0   1   2  ..  H  ....  0 1 2 .. H];
% v = v0 + vw*[0 0 0 .. 0  ....  1 1 1   .. 1   ....  V-1 V-1 V-1 .. V-1 ....  V V V .. V];
%
% The factor for H>1 (V>1) optimizes the inter-subplot distances for hf<1 (vf<1).
% Without that factor the space 1-hf (1-vf) would be lost at 
% the right (top) side of the rightmost (top) plot. 
% Instead this space is divided by the number of inter-subplot distances H-1 (V-1) 
% and distributed equally to augment the inter-subplot distances.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. generate figure with subplots at positions determined by vectors h and v %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if find(mode=='f')
  fh = figure;
else                             
  fh = gcf;                       % last figure opened or new figure, if none exists
end
set(gcf,'units','normalized','position',[0.32 0.05 0.67 0.85]); % ca. port

for hv=1:HV
  sh(hv)=subplot('position',[h(hv) v(hv) hf*hw vf*vw]);
  zoom on; hold on;               % zoom and hold on generally for every subplot
  if find(mode=='g') grid on; end % grid dependent on mode
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. permute the handle-vector sh and ticklabel-cell-arrays according to mode %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sh = reshape(sh,[H V])';          % change sh to matrix for permutations
                                  % transpose because matlab counts quicker vertically
if find(mode=='t')                % top desired instead of default bottom ?
  sh = flipud(sh);                % flip handle matrix up-down
  vtl = fliplr(vtl);              % flip yticklabels-row left-right, i.e. up-down
end %'b'

if find(mode=='r')                % right desired instead of default left ?
  sh = fliplr(sh);                % flip handle matrix left-right
  htl = fliplr(htl);              % flip xticklabels-row left-right
end %'r'

if find(mode=='h')                % quicker counting horizontal desired ?
  sh = sh';                       % transpose sh (after tb- and lr-flips !)
end %'h'

if isempty(find(mode=='m'))       % vector possibly better, caution with matrix
  sh = sh(:);                     % because of different orders in h- and v-case
end % ~'m'                        % normal matrix order in v-mode, reversed for h

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. number subplots in testmode and clear inner ticklabels if desired %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if find(mode=='n')                % number mode, write numbers of subplots into them
  fontsize = floor(min(80,80*15/HV)); % smaller fonts if more subplots
  for hv=1:HV
    subplot(sh(hv)); 
    text(0.5,0.5,num2str(hv),'HorizontalAlignment','center','Fontsize',fontsize);
  end %hv
end %'n'

if find(mode=='c')                % clear inner ticklabels, if desired
  clear_inner_ticklabels(sh);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. outer coordinate system, caution, this inhibits zooming into the subplots %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if find(mode=='o')
  axes('position',[0 0 1 1]); axis([0 1 0 1]); axis('off')
  text(0.98,0.12,hin ,'HorizontalAlignment','right');               % inner xlabel
  text(0.12,0.98,vin ,'HorizontalAlignment','right','rotation',90); % inner ylabel

% text(1.00 ,0.015,hout,'HorizontalAlignment','right','fontsize',16); % outer xlabel, was sometimes out
  text(0.995,0.022,hout,'HorizontalAlignment','right','fontsize',16); % outer xlabel
  text(0.03 ,1.000,vout,'HorizontalAlignment','right','rotation',90,'fontsize',16) % o y

  hl = 0.08; vl=0.08;
  line([hl 0.98],[vl vl  ],'color','k');    % outer horizontal line
  line([hl hl  ],[vl 0.98],'color','k');    % outer vertical   line

  if length(htl)==1                         % only one outer hticklabel provided ?
    htl = repmat(htl,[H 1]);                % replicate it H times
  end
  if length(vtl)==1                         % only one outer vticklabel provided ?
    vtl = repmat(vtl,[V 1]);                % replicate it V times
  end

  for hh=1:min(H,length(htl))
    line([h(hh)+hf*hw/2 h(hh)+hf*hw/2],[vl-0.01 vl+0.01],'color','k');% outer xticks
    text(h(hh)+hf*hw/2,vl-0.04,htl{hh},'HorizontalAlignment','center','fontsize',16);
  end %hh                                                       % outer xticklabels
  for vv=1:min(V,length(vtl))
    line([hl-0.01 hl+0.01],[v(H*vv)+vf*vw/2 v(H*vv)+vf*vw/2],'color','k');%  outer yticks
    text(hl-0.02,v(H*vv)+0.50*vf*vw,vtl{vv},'HorizontalAlignment','right','fontsize',16);
  % text(hl-0.02,v(H*vv)+0.48*vf*vw,vtl{vv},'HorizontalAlignment','right','fontsize',16);
  % note: correct height in matlab-plot with 0.48*vf*vw but:
  %       correct height in eps-file    with 0.50*vf*vw
  end %vv                                                       % outer yticklabels
end %'o'

return
