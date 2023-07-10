5552204function varargout = ida_data_browser(varargin)
% Database breower of IDA 
%
% Part of the IDA Toolbox 
% Copyright (C) Torsten Marquardt
% Terms of the GNU General Public License apply
% (www.http://www.fsf.org/licensing/licenses/gpl.html).
% v01 06Jan23 TM (EI)
%
% TO DO:
%=========================================================================
if nargin == 0
    initialise;
elseif (nargout)
    [varargout{1:nargout}] = feval(varargin{1},varargin{2:end});
else
    feval(varargin{1}, varargin{2:end}); % FEVAL switchyard
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initialise(animalID, config)

if ~exist('animalID','var'), animalID = [];end
if ~exist('config','var'), config = [];end

eval(['cd ' fileparts(which(mfilename))])
h_fig = uifigure;
h_fig.Parent = gcbf;
data_path = getappdata(h_fig.Parent,'data_path');
h_fig.Name = ['IDA Browser ' data_path]; 
h_fig.HandleVisibility = 'on';
h_fig.Tag = 'ida_browser';

files = my_dir(data_path,'mat');
if ~isempty(animalID) % filter files with animalID
    n=0;
    for q=1+size(files,1)
        tmp = strsplit(deblank(files(q,:)),'-');
        if tmp{2} == animalID
            n=n+1;
            files_selected(n,:) = deblank(files(q,:));
        end
    end
    files = files_selected;
end
if ~isempty(config) % filter files with config type (type of experiment)
    n=0;
    for q=1+size(files,1)
        tmp = strsplit(deblank(files(q,:)),'-');
        if tmp{3} == config
            n=n+1;
            files_selected(n,:) = deblank(files(q,:));
        end
    end
    files = files_selected;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CB_show_all

function CB_plot_indvdl
eval(['cd USERS\' userID '\' results.header.config_filename])
eval(['plot_' results.header.config_filename '(results)'])