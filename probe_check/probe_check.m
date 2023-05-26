function probe_check(varargin)
% measures and displays acoustic transfer function of soud delivery system.
%
% Part of the IDA_asio Toolbox
% Copyright (C) Torsten Marquardt
% Terms of the GNU General Public License apply
% (www.http://www.fsf.org/licensing/licenses/gpl.html).
% v0.1 19-07-10 TM (Oldenburg)
% v0.2 03-10-12 TM (EI)
%
% TO DO:
% - include absolute level measuremnt of knowles mic
% - playrecord(): return last the appropriate amountToAllocateSecs bytes of
%   the 1.5*amountToAllocateSecs (compensatinf for dealy)
% - replace reference file list with text box that has callback calling file
%   dialog window
% - unite with calibration software since many routines/functions are
%   similar
% - incorporate additional probe_check into RSP routine

%=========================================================================

cd(fileparts(which(mfilename))), % make sure it runs local

% FEVAL switchyard
if nargin == 0
    initialise(10); % parameter: default number of averages 
elseif (nargout),
    [varargout{1:nargout}] = feval(varargin{1},varargin{2:end});
else
    if isnumeric(varargin{1})
        initialise(varargin{1})
    else
        feval(varargin{1}, varargin{2:end}); % FEVAL switchyard
    end
end,

%=========================================================================
%=========================================================================
function initialise(ext_avg)

if ~exist(['probe_check_data'],'dir')
    eval(['!mkdir ', 'probe_check_data'])
end

h_fig= findobj('Tag','probe_check_fig');
if(isempty(h_fig)),
    open probe_check.fig;
    h_fig= gcf;
end,
figure(h_fig)
set(h_fig,'Toolbar','none')
h_ref_list = findobj(h_fig,'Tag','references');
set(h_ref_list,'Value',1);
probe_check('fill',h_fig),
set(h_fig,'Name',['Probe_check, averages: ' num2str(ext_avg)])
setappdata(h_fig, 'ext_avg', ext_avg)


%=========================================================================
function fill(h_fig)    % fill listbox with references

h_ref_list = findobj(h_fig,'Tag','references');
tmp = what('probe_check_data');
if ~isempty(tmp.mat)
    str = getfield(tmp,'mat');
else
    set(h_ref_list,'String',''),
    return;
end,
files(length(str),30)=' ';
for(n=1:length(str)),
    files(n,1:length(str{n})) = str{n};
end,
str = cellstr(sortrows(files));
set(h_ref_list,'String',str),

% fill comment boxes
val = get(h_ref_list,'Value');
path_str = ['probe_check_data\',str{val}];
try,
    load (path_str,'comment'),
catch,
    h = errordlg(['Unable to open file: ',path_str],...
        'ERROR - Probe_check: fill() listbox');
    return,
end,
if (~exist('comment')),
    comment = {''};
end,
set(findobj(h_fig,'Tag','ref_comment'),'String',comment),
if strcmp(str{val},'00000_NONE__.mat')
    set(findobj(h_fig,'Tag','comment'),'String',comment),
end,
%% ==========================================================================

function callback_retry(h_fig)

if exist('probe_check_data\00000_NONE__.mat','file')
    load('probe_check_data\00000_NONE__.mat','calib_left','calib_right',...
        'H_lsp_left','H_lsp_right','comment','time','fs','min_freq','ext_avg','n_avg'),
    save(['probe_check_data\00000_last_trial.mat'],'calib_left','calib_right',...
        'H_lsp_left','H_lsp_right','comment','time','fs','min_freq','ext_avg','n_avg'),
end

ext_avg = getappdata(gcbf,'ext_avg');
n_avg = 5;

% load current calibration data
if exist('..\calibration\calibration_data\00000_last_calib_left.mat','file')
    load('..\calibration\calibration_data\00000_last_calib_left.mat'),
    calib_left.H_lsp = H_lsp;
    calib_left.H_mic = H_mic;
    calib_left.H_abs_ = H_abs;
    calib_left.max_level = max_level;
    calib_left.min_freq = min_freq;
    calib_left.fs = fs;
    calib_left.comment = comment;
    calib_left.time = time;
    if exist('H_lsp_equalisation','var')
        calib_left.H_eq = H_lsp_equalisation;
    else
        uiwait(errordlg('Equalisation for left speaker not found!'))
        calib_left.H_eq = [];
    end
else,
    errordlg('No calibation file found (calibration\calibration_data\00000_last_calib_left.mat)!')
    return
end
if exist('..\calibration\calibration_data\00000_last_calib_right.mat','file')
    load('..\calibration\calibration_data\00000_last_calib_right.mat'),
    calib_right.H_lsp = H_lsp;
    calib_right.H_mic = H_mic;
    calib_right.H_abs_ = H_abs;
    calib_right.max_level = max_level;
    calib_right.min_freq = min_freq;
    calib_right.fs = fs;
    calib_right.comment = comment;
    calib_right.time = time;
    if exist('H_lsp_equalisation','var')
        calib_right.H_eq = H_lsp_equalisation;
    else
        uiwait(errordlg('Equalisation for right speaker not found!'))
        calib_right.H_eq = [];
    end
else,
    errordlg('No calibation file found (calibration\calibration_data\00000_last_calib_right.mat)!')
    return
end


if (calib_left.fs == calib_right.fs) && ...
        (calib_left.min_freq == calib_right.min_freq)
    fs = calib_left.fs;
    min_freq = calib_left.min_freq;
else
    errordlg('sample_rates or min_freqs between left and right calibration differ)!')
    return
end

% Link PsychPortAudio to the portaudio DLL
olddir = pwd;
drv_path = fileparts(which('portaudio_x86.dll'));
cd(drv_path);
d = PsychPortAudio('GetDevices', 3);
cd(olddir);

l = fs/min_freq; % length of a single sweep
if isempty(calib_left.H_eq),  calib_left.H_eq = ones(l,1); end
if isempty(calib_right.H_eq),  calib_right.H_eq = ones(l,1); end
H_lsp_left = 0;
H_lsp_right = 0;

for q = 1:ext_avg
    set(gcbf,'Name',['Probe_check in progress!  averages: ' ...
        num2str(q) ' (of ' num2str(ext_avg) ')'])
    drawnow
    % flat spectrum noise generation
    [re im] = pol2cart([2*pi*rand(1,l/2) 0],ones(1,l/2+1));
    tmp = (re +im*1i)';
    tmp(1)= 1e-16;
    waveOut = [tmp tmp]./[calib_left.H_eq calib_right.H_eq];
    waveOut = real(ifft([waveOut; conj(waveOut(l/2:-1:2,:))]));
    
    % data aquisition
    h_sound_device = sound_io('open_device',fs, 4, 2);
    waveIn = sound_io('playrecord',h_sound_device,...
        [zeros((n_avg+1)*l,2), repmat(waveOut,(n_avg+1),1)],...
        (n_avg+1)*l/fs);
    sound_io('close_device',h_sound_device);
    if sum(sum(waveIn)) == 0,
        msgbox('Sound output failed! See Matlab window for details!',...
            'ERROR calibration.m:', 'error'),
        return,
    end,

    waveIn = [sum(reshape(waveIn(end-l*n_avg+1:end,1),l,n_avg),2), ...
        sum(reshape(waveIn(end-l*n_avg+1:end,2),l,n_avg),2)];
    H_in_out = fft([waveIn waveOut]);
    H_lsp_left = H_lsp_left + H_in_out(:,1)./H_in_out(:,3);
    H_lsp_right = H_lsp_right + H_in_out(:,2)./H_in_out(:,4);
end
% H_lsp like it would be recorded by the (BK + mic_level_correction):
H_lsp_left = H_lsp_left./calib_left.H_mic/n_avg/ext_avg;
H_lsp_right = H_lsp_right./calib_right.H_mic/n_avg/ext_avg;

% save probe_check data
time = now;
time = now; comment = {''};
save('probe_check_data\00000_NONE__.mat','calib_left','calib_right',...
    'H_lsp_left','H_lsp_right','comment','time','fs','min_freq','ext_avg','n_avg'),

set(findobj(gcbf,'Tag','retry'),'String','Retry'),
set(h_fig,'Name','Probe_check done. OK or RETRY?')
probe_check('fill',h_fig),
probe_check('plot_all_TFs', h_fig),

%% ==========================================================================

function plot_all_TFs(h_fig)

h_plots = findobj('Tag','probe_check_plots');
if(isempty(h_plots)),
    open probe_check_plots.fig;
    h_plots = gcf;
end
figure(h_plots)

set(h_plots,'Name', ... % includes originally ref. magn
    [' Upper - Magnitude; Lower - L/R phase diff.[deg];' ...
    ' current/reference - solid/dotted.']),

% get reference filename
h_references = findobj(h_fig,'Tag','references');
val = get(h_references, 'Value');
str = get(h_references, 'String');
if isempty(str),
    reference = 'probe_check_data\00000_NONE__.mat';
else,
    reference = ['probe_check_data\',str{val}];
end,
set(findobj(h_plots,'Tag','txt_ref_equalisation_file_dates'),'String','')
set(findobj(h_plots,'Tag','txt_current_equalisation_file_dates'),'String','')


if get(findobj(gcbf,'Tag','cb_show_NONE'),'Value')
    plot_TF(h_plots, 'probe_check_data\00000_NONE__.mat', 0)
    plot_TF(h_plots, reference, 1)
else,
    plot_TF(h_plots, reference, 0)
end
%% ==========================================================================

function plot_TF(h_plots, file_name, ref),

Xmin = 10;
Xmax = 20000;
x = linspace(2,Xmax, Xmax - 2);

try
    load(file_name),
catch,
    if ~strcmp(file_name, 'probe_check_data\00000_NONE__.mat')
        msgbox (['File ', file_name, ' not found! Press ''Start'' first'],...
            'ERROR probe_check.m:','error'),
    end
    return,
end,

if ref
    if strfind(file_name, 'NONE__')
        return
    end,
    style = {'b:','r:','k:'};
    set(findobj(h_plots,'Tag','txt_ref_calib_dates_left'),...
        'String',datestr(calib_left.time))
    set(findobj(h_plots,'Tag','txt_ref_calib_dates_right'),...
        'String',datestr(calib_right.time))
else
    style = {'b','r','k'};
    set(findobj(h_plots,'Tag','txt_curr_calib_dates_left'),...
        'String',datestr(calib_left.time))
    set(findobj(h_plots,'Tag','txt_curr_calib_dates_right'),...
        'String',datestr(calib_right.time))
end

% plot magnitude: the equalised loudspeaker repsonse
axes(findobj(h_plots,'Tag','axes_magnitude')),
if ref
    hold on
else
    hold off
end
semilogx(x,20*log10(abs(H_lsp_left(3:Xmax)./calib_left.H_eq(3:Xmax))),style{1},...
    x, 20*log10(abs(H_lsp_right(3:Xmax)./calib_right.H_eq(3:Xmax))),style{2})
axis([Xmin Xmax -20 10])
set(gca,'Tag','axes_magnitude'),
grid on, zoom on,

% plot interaural level difference (ILD) and the diversion from the calibration condition
axes(findobj(h_plots,'Tag','axes_ILD')),
if ref
    hold on
else
    hold off
end
semilogx(x, 20*log10(abs(H_lsp_left(3:Xmax)./calib_left.H_lsp(3:Xmax)./ ...
    H_lsp_right(3:Xmax).*calib_right.H_lsp(3:Xmax))),style{3},...
    x, 20*log10(abs(H_lsp_left(3:Xmax)./calib_left.H_lsp(3:Xmax))),style{1}, ...
    x, 20*log10(abs(H_lsp_right(3:Xmax)./calib_right.H_lsp(3:Xmax))),style{2})
axis([Xmin Xmax -20 10])
set(gca,'Tag','axes_ILD'),
grid on, zoom on,
% plot interaural phase difference (IPD)
axes(findobj(h_plots,'Tag','axes_IPD')),
if ref
    hold on
else
    hold off
end

semilogx(x(48:end), (unwrap(angle(H_lsp_left(50:Xmax))) -  ...
    unwrap(angle(H_lsp_right(50:Xmax))))/pi*180,style{3})
axis([Xmin Xmax -200 200])
set(gca,'Tag','axes_IPD'),
grid on, zoom on,

%% ==========================================================================

function callback_ok(h_fig)

try,
    load(['probe_check_data\00000_NONE__.mat']),
catch,
    msgbox ('Do first calibrate!',...
        'ERROR probe_check.m:','error'),
    return,
end,
comment = char(get(findobj(h_fig,'Tag','comment'),'String'));
[m n] = size(comment);
comment = cellstr([[datestr(now),zeros(1,n-length(datestr(now)))];...
    [comment,zeros(m,length(datestr(now))-n)]]);

str = ['probe_check_data\00000_last_probe_check'];
save(str,'H_lsp_left','H_lsp_right','comment','time','fs','min_freq', ...
    'ext_avg','n_avg','calib_left','calib_right','H_lsp_left','H_lsp_right', ...
    'comment','time','fs','min_freq','ext_avg','n_avg'),

new_reference = get(findobj(h_fig, 'Tag','new_name'),'String');
if ~isempty(new_reference)
    str = ['probe_check_data\00_',new_reference,'.mat'];
    if exist(str,'file'),
        if (~strcmp(questdlg(strvcat('File: ', str, ...
                ' already exists. Save anyway?'),...
                'WARNING - Probe_check:','warn'),'Yes')),
            return,
        end,
    end,
    try,
        save(str,'H_lsp_left','H_lsp_right','comment','time','fs','min_freq', ...
            'ext_avg','n_avg','calib_left','calib_right','H_lsp_left','H_lsp_right',...
            'comment','time','fs','min_freq','ext_avg','n_avg'),
    catch,
        msgbox ('Reference name must be valid File name !',...
            'ERROR probe_check.m:','error'),
        return,
    end,
end,

str = ['probe_check_data\probe_check_' ,datestr(time,'yyyymmddTHHMMSS')];
save(str,'H_lsp_left','H_lsp_right','comment','time','fs','min_freq', ...
    'ext_avg','n_avg','calib_left','calib_right','H_lsp_left','H_lsp_right',...
    'comment','time','fs','min_freq','ext_avg','n_avg'),

probe_check('fill',h_fig),
%% ==========================================================================

function callback_ref_comment(h_fig) % Callback of 'comment'-editbox of reference
h = findobj(h_fig,'Tag','references');
str = get(h, 'String');
if (isempty(str)), return, end,
val = get(h, 'Value');
comment = get(findobj(h_fig,'Tag','ref_comment'),'String');
if ~iscell(comment),
    comment = comment;
end,
save(['probe_check_data\',str{val}],...
    'comment','-append'),
%% ==========================================================================

function callback_delete_reference(h_fig)
h = findobj(h_fig,'Tag','references');
str = get(h, 'String');
if (isempty(str)), return, end,
val = get(h, 'Value');
item = str{val};
if (strcmp(questdlg(strvcat('Are you sure you want delete : ',...
        item),'WARNING - Probe_check:','warn'),'Yes')),
    delete (['probe_check_data\',item])
    set(h,'Value',1);
    probe_check('fill',h_fig),
end
%% ==========================================================================

function deleteFcn
h = findobj('Tag','probe_check_plots');
if (~isempty(h)), close(h), end,
% if exist('probe_check_data\00000_NONE__.mat','file')
%     delete ('probe_check_data\00000_NONE__.mat')
% end
cd('..')

%% ==========================================================================

function help
msgbox (strvcat('ERROR probe_check.m:',...
    'Not implemented yet. Sorry!')),


A