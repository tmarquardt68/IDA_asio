function varargout = search(varargin)
% Search stimulus generation for the IDA_asio
%
% Part of the IDA_asio Toolbox
% Copyright (C) Torsten Marquardt
% Terms of the GNU General Public License apply
% (www.http://www.fsf.org/licensing/licenses/gpl.html).
% v0.1 19-07-10 TM (Oldenburg)
% v1.1 18-10-12 TM (EI)
%
% TO DO:
% - save into data directory with filename convention
% - SAMPLE_RATE issue
% - include narrow band noise as search stimulus
% - check for levels, calibrate! (new lsp_transfer_fcn.mat!)
% - reverse order of electrodes (nimbers in legend inverse to selection in
%   drop-down menu);
%   define better color map by hand  to include more colors and avoid 
%   repetition of colors. currently "lines" . where? - search for "lines"
% - make changes to edit background color when clicking radio buttons (not
%   just in the while loop)
% - get max_tone_level from stimulus_parameter and set YTickLabel
%   automatically
%=========================================================================

if nargin == 0
    initialise;
elseif (nargout),
    [varargout{1:nargout}] = feval(varargin{1},varargin{2:end});
else
    feval(varargin{1}, varargin{2:end}); % FEVAL switchyard
end,
%% ========================================================================
%% ========================================================================

function ButtonPress, % mouse click within the axes area

h_mute_left = getappdata(gcbf,'h_mute_left');
h_mute_right = getappdata(gcbf,'h_mute_right');
h_attn_left = getappdata(gcbf,'h_attn_left');
h_attn_right = getappdata(gcbf,'h_attn_right');

selectionType = get(gcbf,'SelectionType');

% toggle muting of left channel
if strcmp(selectionType, 'normal'), % left mouse click contolls muting of Ch 1
    if get(h_mute_left,'Value'),
        set(h_mute_left,'Value',0),
        set(h_attn_left,'BackgroundColor',[1 0 0])
    else,
        set(h_mute_left,'Value',1),
        set(h_attn_left,'BackgroundColor',[1 1 1])
    end,
end,

% toggle muting of right channel
if strcmp(selectionType, 'alt'), % right mouse click
    if get(h_mute_right,'Value'),
        set(h_mute_right,'Value',0),
        set(h_attn_right,'BackgroundColor',[1 0 0])
    else,
        set(h_mute_right,'Value',1),
        set(h_attn_right,'BackgroundColor',[1 1 1])
    end,
end,

% set marker
if strcmp(selectionType, 'open'), % double click
    
    % undo the triggered 'normal' mouse click (double click causes 'normal'
    % followed 'open' button event
    if get(h_mute_left,'Value'),
        set(h_mute_left,'Value',0),
        set(h_attn_left,'BackgroundColor',[1 0 0])
    else,
        set(h_mute_left,'Value',1),
        set(h_attn_left,'BackgroundColor',[1 1 1])
    end,
    % end undo
    
    electrode_no = get(findobj(gcbf,'Tag','pop_electrode_no'),'String');
    cmap = colormap;
    marker_color = cmap(17-str2num(electrode_no{ ...
        get(findobj(gcbf,'Tag','pop_electrode_no'),'Value')}),:);
    
    tmp = round(get(gca,'CurrentPoint')); % set a circle at mouse postion
    bf = tmp(1,1);
    attn = tmp(1,2);
    
    % take face off previous marker
    h_markers = getappdata(gcbf,'h_markers');
    if ~isempty(h_markers) && ishandle(h_markers(end)) && ...
            ~sum(marker_color - get(h_markers(end),'FaceColor')) && ...
             ~sum([1 1] - get(h_markers(end),'Curvature'))
        set(h_markers(end),'FaceColor','none','EdgeColor',marker_color)
    end
    
    % plot new marker
    h_markers(end+1)=rectangle('Position',[bf-bf/30, ...
    attn-1, bf/15, 1.5],...
        'Curvature', [1,1],...
        'LineWidth', 3,...
        'FaceColor',marker_color, ...
        'EdgeColor',marker_color, ...
        'ButtonDownFcn','search ButtonPress');
    set(gcbf, 'Name', ['SEARCH for IDA,  BF: ' num2str(bf) ...
        ' Hz,   Attn: ' num2str(attn) ' dB'])
    setappdata(gcbf,'h_markers',h_markers)
    setappdata(gcbf,'bf',bf)
    setappdata(gcbf,'attn',attn)
end,

% delete marker
if strcmp(selectionType, 'extend'), % mouse click + SHIFT , or both buttons
    h = gco;
    if strcmp(get(h,'Type'),'rectangle')
        h_markers = getappdata(gcbf,'h_markers');
        delete(gco)
        [a index] = min(abs(h_markers - h));
        h_markers = [h_markers(1:index-1) h_markers(index+1:end)];
        setappdata(gcf,'h_markers',h_markers)
    end
end,

% =========================================================================
function CB_clear 

h_markers = getappdata(gcbf,'h_markers');
delete(h_markers)
h_markers = [];
setappdata(gcbf,'h_markers',h_markers)

% =========================================================================
function CB_load 

h_markers = getappdata(gcbf,'h_markers'),

[filename,pathname] = uigetfile('Load markers:');

if ~isempty(filename) & filename,
    try
    load([pathname,filename],'marker_positions','marker_faceColors',...
        'marker_edgeColors')
    end,
    
    if ~exist('marker_positions', 'var')
        uiwailt(errordlg('File does not contain marker data','search.m'))
        return
    end

    for q=1:length(marker_positions)
        h_new_markers(q)=rectangle('Position',marker_positions{q},...
            'Curvature', [0,0],...
            'LineWidth', 3,...
            'EdgeColor',marker_edgeColors{q}, ...
            'FaceColor',marker_faceColors{q}, ...
            'ButtonDownFcn','search ButtonPress');
    end
    
    setappdata(gcbf,'h_markers',[h_markers h_new_markers])
end

% =========================================================================
function CB_mute_left

h_attn_left = getappdata(gcbf,'h_attn_left');

if get(gcbo,'Value'),
    set(h_attn_left,'BackgroundColor',[1 1 1])
else,
    set(h_attn_left,'BackgroundColor',[1 0 0])
end,

% =========================================================================
function CB_mute_right

h_attn_right = getappdata(gcbf,'h_attn_right');

if get(gcbo,'Value'),
    set(h_attn_right,'BackgroundColor',[1 1 1])
else,
    set(h_attn_right,'BackgroundColor',[1 0 0])
end,

% =========================================================================
function CB_save 

h_markers = getappdata(gcbf,'h_markers');
filename = getappdata(gcbf,'filename');

if isempty(filename)
    filename = fileparts(which(mfilename));
end
if strcmp(filename(end-3:end), '.mat')
    filename = filename(1:end-4);
end

[filename,pathname] = uiputfile(filename(1:end-1),'Save markers as:');

marker_positions = get(h_markers,'Position');
marker_faceColors = get(h_markers,'FaceColor');
marker_edgeColors = get(h_markers,'EdgeColor');

if ~isempty(filename) & filename,
    save([pathname,filename],'marker_positions','marker_faceColors',...
        'marker_edgeColors')
    setappdata(gcbf,'filename',[pathname,filename])
end


% =========================================================================
function CB_start_stop

h_sound_device = getappdata(gcbf,'h_sound_device');

if strcmp(get(gcbo, 'String'),'start')
    set(gcbo, 'String','stop')
    % future: return to standby (PsychPortAusdio('RunMode', ...)
    play_loop(gcbf,gcbo)
else
    set(gcbo, 'String','start')
    % future: go to standby
end
% =========================================================================
function [bf attn] = closeRequestFcn % when figure gets deleted

h_fig = findobj('Tag','search_fig');
figure(h_fig)
bf = getappdata(h_fig, 'bf');
attn = getappdata(h_fig, 'attn');
h_sound_device = getappdata(h_fig,'h_sound_device');
h_stim = getappdata(h_fig, 'h_stimuli_gui');

closereq,
pause(0.1), % Let the while loop stop before closing
if ~isempty(h_sound_device),
    sound_io('close_device', h_sound_device);
end,
if ishandle(h_stim)
    close(h_stim)
end,

%set BF and THR in IDA figure
h_ida_fig = findobj('Tag','ida_fig');
if isempty(h_ida_fig),
    ida;
    h_ida_fig = findobj('Tag','ida_fig');
else,
    figure(h_ida_fig),
end,
if ~isempty(bf) && ~isempty(h_ida_fig),
    set(findobj(h_ida_fig,'Tag','bf'),'String',num2str(5*round(bf/5))),
    set(findobj(h_ida_fig,'Tag','attn'),'String',num2str(-5*round(attn/5))),
end,
cd('..')


% =========================================================================
function initialise(username)

global SAMPLE_RATE

SAMPLE_RATE = 48000;

p = fileparts(which(mfilename));
[tmp p] = strtok(p(end:-1:1),'\');
p = p(end:-1:1);
% addpath(p, [p 'tm\RSP_10'], ...
%     [p 'probe_check'], ...
%     [p 'tm\NDFs_10'], ...
%     [p 'tm\RIF_10'], ...
%     [p 'tm\IAC_JND_10'], ...
%     [p 'probe_check'], ...
%     [p 'stimuli'])
    
% == Default dimensions of axis ===========================================
% Overwrite these by loading a configurartion MAT-file 'search_config'
config_filename = 'configs\search_config.mat';
freq_scale = round(logspace(log10(100), log10(4*6400), 9));
attn_min = -90;
attn_max = 0;
defaults.stim_duration = 150;
defaults.isi = 180;
defaults.max_level_tone = 95;
defaults.max_level_noise = 63;
% end of  defaults definitions
yTickLabels = {'05' '15' '25' '35' '45' '55' '65' '75' '85' '95'};



cd(fileparts(which(mfilename))), % make sure it runs local
if exist('username','var'),
    try % load individual configuration of search
        load([username '\search_config.mat']); % load configutaion MAT file. See general comment.
        config_filename = [username '\search_config.mat'];
    catch
        uiwait(warndlg(['Could not load config file: ' username ...
            '\search_config.mat'], ...
            'IDA WARNING','modal'));
    end,
end,
freq_min = freq_scale(1);
freq_max = min(23000,freq_scale(end));

set(0,'Units', 'normalized'),

% ==== set up search figure ======
h_fig = findobj('Tag','search_fig');
if isempty(h_fig)
    h_fig = openfig('search.fig');
    set(h_fig,'Position',[0 0 0.4 0.4]),
else
    figure(h_fig),
end
set (gca,'Xtick',[freq_scale(1:end-1) freq_max])
set(gca,'YTickLabel',yTickLabels)
axis([freq_min, freq_max, attn_min, attn_max]), grid on,
colormap(lines(16))
colorbar('Ytick',1:16)

h_mute_left = findobj(h_fig,'Tag','RB_mute_left');
h_mute_right = findobj(h_fig,'Tag','RB_mute_right');
h_attn_left = findobj(h_fig,'Tag','txt_attn_left');
h_attn_right = findobj(h_fig,'Tag','txt_attn_right');

set([h_mute_right h_mute_left],'Value',1)
set([h_attn_right h_attn_left],'BackgroundColor',[1 1 1])

% ==== Initialise sound output device
h_sound_device = sound_io('open_device', 48000, 4, 0);

% ==== initialise stimuli_gui
set(h_fig,'Units', 'pixel'),
position_fig = get(h_fig, 'Position');
set(h_fig,'Units', 'normalized'),

h_stim = getappdata(h_fig,'h_stimuli_gui');
if isempty(h_stim) || ~ishandle(h_stim) || ...
    ~strcmp(get(h_stim,'Tag'),'stimuli_gui_fig'),
    h_stim = stimuli_gui(h_fig, ...
        [position_fig(1), position_fig(4)+50],config_filename);
end

if stimuli_gui('Change_parameter', h_stim, [], ...
        'duration [ms]', defaults.stim_duration)
    error('Paramter change failed.')
end
if stimuli_gui('Change_parameter', h_stim, 0, ...
        'ISI [ms]', defaults.isi)
    error('Paramter change failed.')
end

if ishandle(h_stim),
    setappdata(h_fig, 'h_stimuli_gui', h_stim);
else,
    uiwait(warndlg(['The stimuli_gui could not be opened!'], ...
        'IDA SEARCH WARNING','modal'));
    closeRequestFcn;
    return
end,

setappdata(h_fig, 'h_sound_device',h_sound_device)
setappdata(h_fig,'min_max_axis',[freq_min freq_max attn_min attn_max])
setappdata(h_fig,'defaults',defaults)
setappdata(h_fig,'h_stim',h_stim)
setappdata(h_fig,'h_markers',[])
setappdata(h_fig,'h_mute_left',h_mute_left)
setappdata(h_fig,'h_mute_right',h_mute_right)
setappdata(h_fig,'h_attn_left',h_attn_left)
setappdata(h_fig,'h_attn_right',h_attn_right)

h_start_stop = findobj(h_fig,'Tag','PB_start_stop');
set(h_start_stop, 'String','stop')
play_loop(h_fig, h_start_stop)

% =========================================================================
function play_loop(h_fig, h_start_stop)

global SAMPLE_RATE

min_max_axis = getappdata(h_fig,'min_max_axis');
defaults = getappdata(h_fig,'defaults');
h_stim = getappdata(h_fig,'h_stim');
h_mute_left = getappdata(h_fig,'h_mute_left');
h_mute_right = getappdata(h_fig,'h_mute_right');
h_attn_left = getappdata(h_fig,'h_attn_left');
h_attn_right = getappdata(h_fig,'h_attn_right');
h_sound_device = getappdata(h_fig,'h_sound_device');

h_noise = findobj(h_fig,'Tag','RB_noise');
h_binbeat = findobj(h_fig,'Tag','RB_binbeat');
h_search_axis = findobj(h_fig,'Tag','axis_search');
h_freq = findobj(h_fig,'Tag','txt_freq');

tic
% === Contineous update of stimuli_gui
while ishandle(h_start_stop) && strcmp(get(h_start_stop, 'String'),'stop')
    if ~ishandle(h_fig)
        sound_io('close_device',h_sound_device)
        return,
    end
    stimuli_parameter = stimuli_gui('Get_stimuli_parameter', h_stim);
    
    % get position of mouse pointer
    ptr_location = get(0,'PointerLocation');
    fig_pos = get(h_fig,'Position');
    noise = get(h_noise,'Value'); %  radio button noise?
    binbeat = get(h_binbeat,'Value'); %  radio button binaural beat?
    ax_pos = get(h_search_axis,'Position');
    x = ((ptr_location(1) - fig_pos(1))/fig_pos(3)- ax_pos(1))/ax_pos(3);
    y = ((ptr_location(2) - fig_pos(2))/fig_pos(4) - ax_pos(2))/ax_pos(4);
    % translate position to frequency and attentuation
    freq1 = round(min_max_axis(1)*(min_max_axis(2)/min_max_axis(1))^x);
    freq2 = freq1;
    attn = -round(min_max_axis(3) + y*(min_max_axis(4) - min_max_axis(3)));

    if x<0 || x>1 ||y<0 || y>1, % if mouse outside axis
        attn = 120;
        set(h_freq, 'String', '----'),
        set(h_freq,'BackgroundColor',[1 0 0])
        set(h_attn_right,'BackgroundColor',[1 0 0])
        set(h_attn_left,'BackgroundColor',[1 0 0])
        outside_fig = 500/stimuli_parameter.isi;
    elseif outside_fig > 0
        attn = 120;
        outside_fig = outside_fig - 1;
    end,

    stim_duration = defaults.stim_duration;
    
    if attn < 120,
        figure(stimuli_parameter.stimuli(2).h_stimulus_gui)
        phase_tone2 = 0;
        interval = 'B';
        BackgroundColor_freq = [1 1 1];
        attn_tone2 = 999;
        pink_spctrm = 1;
        ild1 = 0;
        ild2 = 0;

        if binbeat,
            set(h_noise, 'Value', 0),
            noise = 0;
            stim_duration = 2000;
            freq1 = freq1 - 1;
            freq2 = freq2 + 1;
            BackgroundColor_freq = [1 1 0];
            figure(stimuli_parameter.stimuli(2).h_stimulus_gui)
            attn = attn+6;
            attn_tone2 = attn;
            ild1 = 99;
            ild2 = -99;
        end

        if noise,
            BackgroundColor_freq = [0 0 1];
            figure(stimuli_parameter.stimuli(1).h_stimulus_gui)
            interval = 'A';
            set(h_binbeat, 'Value', 0),
        end

        % set new stim parameter
        err = 0;
        err = err + stimuli_gui('Change_parameter', h_stim, [], ...
            'duration [ms]', stim_duration);
        err = err + stimuli_gui('Change_parameter', h_stim, ...
            1,'spct_level [dB]', max(-40,defaults.max_level_noise - attn));
        err = err + stimuli_gui('Change_parameter', h_stim, ...
            2,'level [dB]', max(-40,defaults.max_level_tone - attn));
        err = err + stimuli_gui('Change_parameter', h_stim, ...
            3,'level [dB]', max(-40,defaults.max_level_tone - attn_tone2));
        err = err + stimuli_gui('Change_parameter', h_stim, ...
            2,'freq [Hz]', freq1);
        err = err + stimuli_gui('Change_parameter', h_stim, ...
            3,'freq [Hz]', freq2);
        err = err + stimuli_gui('Change_parameter', h_stim, ...
            2,'Phase [deg]', phase_tone2);
        err = err + stimuli_gui('Change_parameter', h_stim, ...
            1,'pink_spctrm', pink_spctrm);
        err = err + stimuli_gui('Change_parameter', h_stim, ...
            2,'ILD [dB]', ild1);
        err = err + stimuli_gui('Change_parameter', h_stim, ...
            3,'ILD [dB]', ild2);

        if err
            error('Paramter change failed.')
        end

        wave = stimuli_gui('Get_waveform',interval,h_stim);
        if isempty(wave)
            error(['stimulus wave is empty. Is stimulus selected for interval ' interval '?'])
        end
        
        % update serach window
        set(h_freq,'BackgroundColor',BackgroundColor_freq),
        if sum(BackgroundColor_freq) == 1 
            set(h_freq, 'String', '----'),
        else
            set(h_freq, 'String', num2str((freq1+freq2)/2)),
        end

        % check for muting
        if get(h_mute_left,'Value'),
            set(h_attn_left, 'String', num2str(attn)),
            set(h_attn_left,'BackgroundColor',[1 1 1])

        else,
            wave(:,1) = 0;
            set(h_attn_left,'BackgroundColor',[1 0 0])
            set(h_attn_left, 'String', num2str(120)),
        end,

        if get(h_mute_right,'Value'),
            set(h_attn_right, 'String', num2str(attn)),
            set(h_attn_right,'BackgroundColor',[1 1 1])

        else,
            wave(:,2) = 0;
            set(h_attn_right, 'String', num2str(120)),
            set(h_attn_right,'BackgroundColor',[1 0 0])

        end,
    end,
    
    pause((stimuli_parameter.isi+stim_duration)/1000-toc)
    tic

    % present wave form
    if attn < 120,   
        sound_io('play',h_sound_device,[zeros(size(wave,1),2), wave]);
    end,
end, % while

