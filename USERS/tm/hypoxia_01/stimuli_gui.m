function varargout = stimuli_gui(varargin),
% Graphical user interface to combine several Stimulus GUIs
%
% External usage:
% Get_waveform_XXX() returns the waveform precomputed and stored in the
% stimulus_guis.
%
% Hint: After pressing "save, "plot" or "play_XX" button, the waveform will
% be available in variable "ans" of the current workspace (command window).

% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Version history:
% 1.1  - new functions: Set_interval_specs, Get_interval_specs,
%      - resets levelmeter before playing next sound (not red overload LED)
% 1.2  - reopens with the same figure handle, no FIG-file anymore
% 1.3  - Adaptations for PSYCHO application
%      - functions which are also intended for external calls begin with a
%        capital letter.
% 2.0  - 'User' dat of figure holds the whole stimulus_parameter structure
%      - stimuli_gui.m now in same directory like stimulus scripts (not in
%        '\General\' anymore.
% 2.1  - introduced fileparts.m to get pathname
% 2.2  - Get_stimuli_parameter() fills
%        stimuli_parameter.stimuli(q).h_stimulus_gui       (remained empty
%        before. Needed now in afc_staircase.m play_waveform())
%      - Select_edit_uicontrol() returns now fields
%        'variable.stimulus.interval' and 'variable.parameter.name'
% 3.0  - allows frozen or freshly generated noise in Get_waveform_XXX
%      - precomputed waveforms are stored in 'User' of playbuttons (as in
%        stimulus_gui.m.
%      - A{AB}A shortened to A{AB}
% 4.0  - storage of object data using setappdata() and getappdata(), i.e.
%        not in 'User' variables.
%      - major rewrite of Precompute() => refresh_waveform() and
%		 update_waveform(). Refreshed waveforms{} stored as appdata.
%	   - Get_waveform_X calls these functions if parameter has changed.
%      - Get_waveform_X assembles now trials from waveforms{} and 'isi'
%       (read out from GUI) only on demand.
% 4.1  - some adaptations to work with IDA asio
%      - Change_parameter() extended to work for stimulus GUIs parameter
%      - config folder is now one level up
%
% Part of the Stimulus Toolbox
% Copyright (C) Torsten Marquardt
% Terms of the GNU General Public License apply
% (www.http://www.fsf.org/licensing/licenses/gpl.html).

vers = '4.1';


if nargin == 0, % initialisation without parameters (usually command line)
    h_fig = initialise;
    if nargout > 0, varargout{1} = h_fig; end,
elseif ~ischar(varargin{1}), %initialisation with parameters (parent appl.)
    h_fig = initialise(varargin{1:end});
    if nargout > 0, varargout{1} = h_fig; end,
else, % INVOKE NAMED SUBFUNCTION OR CALLBACK
    if strcmp(varargin{1},'open_stimuli_gui'), % especially to pass variable 'vers'
        varargout{1} = open_stimuli_gui(vers,varargin{2:end});
    else,
        if (nargout),
            [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        else,
            feval(varargin{:}); % FEVAL switchyard
        end,
    end,
end,
%% ==========================================================================
%% ==========================================================================
%% ==========================================================================

%% ==========================================================================

function ButtonDownFnc_figure,

% windowButtonUpFnc will  bring stimulus GUIs to foreground
switch get(gcbf,'SelectionType')
    case 'normal', % left mouse click: store current position as User data of text_isi
        setappdata(gcbf,'position',get(gcbf,'CurrentPoint')),
    case 'alt', % right mouse click:

    case 'open', % double click: Reset level meter
        Level_meter(gcbf, 'full_reset');
end,
%% ==========================================================================

function ButtonDownFnc_text,
% set associated edit uicontrol to CurrentObject.
% delete dialog box if opened by fcn 'Select_edit_uicontrol'
% double click: reset level meter.

set(gcbf,'CurrentObject', findobj(gcbf,'Tag','edit_isi')),
% delete dialog box if opened by fcn 'Select_edit_uicontrol'
delete(findobj(0,'Tag','stimulus_gui_select_edit_uicontrol_fig')),
switch get(gcbf,'SelectionType')
    case 'normal', % left mouse click:
    case 'alt', % right mouse click: Close all stimuli GUI
    case 'open', % double click: Reset level meter
        Level_meter(gcbf, 'full_reset');
end,
%% ==========================================================================

function callback_checkbox_A,

stimuli_parameter = getappdata(gcbf,'stimuli_parameter');
checkbox_tag = get(gcbo,'Tag');

q = str2num(fliplr(strtok(fliplr(checkbox_tag),'_')));
if ~isempty(stimuli_parameter.stimuli(q).specs),
    stimuli_parameter.stimuli(q).interval_A = get(gcbo,'Value');
end,

setappdata(gcbf,'stimuli_parameter',stimuli_parameter),
Refresh_waveforms(gcbf,[2 3]),
%% ==========================================================================

function callback_checkbox_B,

stimuli_parameter = getappdata(gcbf,'stimuli_parameter');
checkbox_tag = get(gcbo,'Tag');

q = str2num(fliplr(strtok(fliplr(checkbox_tag),'_')));
if ~isempty(stimuli_parameter.stimuli(q).specs),
    stimuli_parameter.stimuli(q).interval_B = get(gcbo,'Value');
end,

setappdata(gcbf,'stimuli_parameter',stimuli_parameter),
Refresh_waveforms(gcbf,1),
%% ==========================================================================

function waveform = callback_plot_A,% left click: plot wave form in new figure
% Not finished! % right click: open context menu with plot options

waveform = plot_interval(gcbf,'A');
%% ==========================================================================

function waveform = callback_plot_B,
% left click: plot wave form in new figure
% Not finished! % right click: open context menu with plot options

waveform = plot_interval(gcbf,'B');
%% ==========================================================================

function waveform = callback_save_A,
% Not finished!
% left click: open file dialog box, save wave form as WAV file
% Not finished! % right click: open context menu with saving options

waveform = Save_waveform(gcbf,'A');
%% ==========================================================================

function waveform = callback_save_B,
% Not finished!
% left click: open file dialog box, save wave form as WAV file
% right click: open context menu with saving options

waveform = Save_waveform(gcbf,'B');
%% ==========================================================================

function callback_popupmenu,

popupmenu_tag = get(gcbo,'Tag');
q = str2num(fliplr(strtok(fliplr(popupmenu_tag),'_')));
selected_stimulus = get_menu_entry(gcbo);
stimuli_parameter = Get_stimuli_parameter(gcbf);

if strcmp(selected_stimulus,'Remove Stimulus'),
    q3 = 1;
    for q2 = 1:length(stimuli_parameter.stimuli),
        if q2 ~= q,
            stimuli_new(q3) = stimuli_parameter.stimuli(q2);
            q3 = q3 + 1;
        end,
    end,
    if q2 > 1, % keep at least one popupmenu
        stimuli_parameter.stimuli = stimuli_new;
    end,
else
    if strcmp(strtok(selected_stimulus),'Like'),
        if ~isempty(stimuli_parameter.stimuli(1).specs)
            stimuli_parameter.stimuli(q).specs=stimuli_parameter.stimuli(...
                str2num(fliplr(strtok(fliplr(selected_stimulus))))).specs;
        end,
    else,
        selected_stimulus = strtok(selected_stimulus,'.');
        eval(['stimuli_parameter.stimuli(q).specs = ', ...
            fliplr(strtok(fliplr(selected_stimulus),':')), ...
            '(''return_default_parameter'');']),
    end,
end,
reopen_GUIs(gcbf,stimuli_parameter),
%% ==========================================================================

function err = Change_parameter(h_fig, stimulus_no, parameter, ...
    new_value)
% !!! Use this routine to change paramter from external applications !!!
% check new value, when failure keep old value.
% $stimulus_no - is index in stimuli_parameter.stimuli struct.
% If  stimulus_no == 0: the parameter is in the stimuli_gui (e.g. 'isi')
% If  stimulus_no == []: change in all stimulu_guis the parameter named
% editUI
% TO DO: include paramter 'frozen'

stimuli_parameter = getappdata(h_fig,'stimuli_parameter');

verbose = 1;
if ~exist('new_value','var') % manual parameter change in stimuli_gui GUI
    new_value=str2num(get(gcbo,'String'));
    stimulus_no = 0;
    parameter = gcbo;
end,

if stimulus_no == 0 % The parameter is in the stimuli_gui (e.g. 'isi')
    err = change_stimuli_gui_parameter(h_fig,parameter,new_value);
else % It's a parameter of one of the stimulus_guis
    if isempty(stimulus_no) % change in all stimulus GUIs
        stimulus_no = 1:length(stimuli_parameter.stimuli);
        verbose = 0;
    end

    for q = stimulus_no
        err(q) = stimulus_gui('Change_parameter', ...
            stimuli_parameter.stimuli(q).h_stimulus_gui, ...
            parameter, new_value, verbose);
    end,
end

if verbose && sum(err),
    err = 1;
else
    err = 0;
end,



%% ==========================================================================

function h = change_stimuli_gui_parameter(h_fig,parameter,new_value)
% 'parameter' can be handle or text string of parameter in GUI

h = 0; % if h > 0 error occured
if ischar(parameter)
    % make here a case list of possible stimuli_gui parameter
    % At the moment there is only the parameter ISI
    switch parameter
        case 'ISI [ms]'
            parameter = findobj(h_fig, 'Tag', 'edit_isi');
        otherwise
            h = errordlg(['No such stimuli_gui paramter:' parameter],...
                'stimulus_gui.m - change_stimuli_gui_parameter');
            uiwait(h),
            return
    end
end

if ishandle(parameter)
    if isempty(new_value)|| length(new_value)>1,
        pause(0.2)
        h = errordlg('Parameter must be a single number!' ...
            ,'stimulus_gui.m - change_stimuli_gui_parameter');
        uiwait(h),
        new_value = get(gcbo,'User');
    end
    set(parameter,'User',new_value),
    set(parameter,'String',num2str(new_value)),
else
    h = errordlg('Invalid stimuli_gui handle!',...
        'stimulus_gui.m - change_stimuli_gui_parameter');
    uiwait(h),
end
%% ==========================================================================

function mixes = create_interval_mixes(h_fig, interval_nos)
% To be done: leave waveforms of interval B stored in stimulus_guis
% To be done: frozen == 1

waveforms = getappdata(h_fig,'waveforms');
stimuli_parameter = getappdata(h_fig,'stimuli_parameter');
frozen = get(findobj(h_fig,'Tag','checkbox_frozen'),'Value');


% get list of refreshed stimulus_guis relevant for the interval_nos
if max(interval_nos == 1),
    interval_B = 1;
else
    interval_B = 0;
end,
if max(interval_nos > 1),
    interval_A = 1;
else
    interval_A = 0;
end,

list_h_stimulus_gui_dirty_new = [];

for q = 1:length(stimuli_parameter.stimuli);
    if (interval_A & stimuli_parameter.stimuli(q).interval_A) | ...
            (interval_B & stimuli_parameter.stimuli(q).interval_B),

        if stimulus_gui('Check_and_reset_dirty_flag', ...
                stimuli_parameter.stimuli(q).h_stimulus_gui, h_fig),

            list_h_stimulus_gui_dirty_new=[list_h_stimulus_gui_dirty_new, ...
                stimuli_parameter.stimuli(q).h_stimulus_gui];
        end,
    end,
end,

% get update for all dirty waveforms
mixes = [];
q3 = 0;
for q = interval_nos;
    q3 = q3+1;  stimuli = [];
    waveforms(q).dirty_h_stimulus_list = ...
        [waveforms(q).dirty_h_stimulus_list list_h_stimulus_gui_dirty_new];
    clear('stimuli'),
    for q2 = 1:length(waveforms(q).stimuli),
        if max(waveforms(q).stimuli(q2).h_stimulus_gui == ...
                waveforms(q).dirty_h_stimulus_list),

            if ~frozen
                stimulus_gui('Refresh_waveform', ...
                    waveforms(q).stimuli(q2).h_stimulus_gui),
            end,
            stimuli(q2).h_stimulus_gui = ...
                waveforms(q).stimuli(q2).h_stimulus_gui;
            [stimuli(q2).waveform stimuli(q2).levels] = ...
                stimulus_gui('Get_waveform', ...
                waveforms(q).stimuli(q2).h_stimulus_gui);
            % reset dirty_flag in stimulus_gui also in case Get_waveform()
            % caused refreshing and has set the flag again
            stimulus_gui('Check_and_reset_dirty_flag', ...
                waveforms(q).stimuli(q2).h_stimulus_gui, h_fig);
        else,
            stimuli(q2) = waveforms(q).stimuli(q2);
        end,
        sizes(q2,:) = size(stimuli(q2).waveform);
        mixes(q3).level_struct.stimuli(q2,:,:) = stimuli(q2).levels;
    end,

    if exist('stimuli','var'),
        stimuli = equalise_Nchannels(stimuli);

        % create interval mixes
        if isstruct(stimuli),

            % adjust waveforms to interval length
            max_size = max(sizes,[],1);
            mixes(q3).waveform = zeros(max_size);
            for q2 = 1:length(waveforms(q).stimuli),
                mixes(q3).waveform=mixes(q3).waveform + [stimuli(q2).waveform;...
                    zeros(max_size(1)-sizes(q2,1), max_size(2))];
            end,

            % set n_samples and levels
            mixes(q3).n_samples = size(mixes(q3).waveform,1);
            max_levels = 20*log10(max(abs(mixes(q3).waveform)));
            for q2 = 1:length(max_levels),
                if max_levels(q2) == 0, % log(0) gives warning!
                    rms_levels(q2) = max_levels(q2);
                else,
                    rms_levels(q2) = ...
                        20*log10(sqrt(mean(mixes(q3).waveform(:,q2).^2)));
                end,
            end,
            mixes(q3).level_struct.mix = [max_levels; rms_levels];
        else,
            mixes(q3) = stimuli;
        end,
        waveforms(q).stimuli = stimuli;
    end,
end,

if ~isempty(mixes),
    mixes = equalise_Nchannels(mixes);
end,

waveforms = set_dirty_list(h_fig, waveforms, interval_nos, ...
    list_h_stimulus_gui_dirty_new);
setappdata(h_fig,'waveforms',waveforms),
%% ==========================================================================

function delete_stimuli_gui(h_fig),
% CloseRequestFcn & external call

h_fig_parent = getappdata(h_fig, 'h_fig_parent');
if isempty(h_fig_parent)| ~ishandle(h_fig_parent),
    delete_stimulus_guis(h_fig),
    figure(h_fig)
    closereq,
end,
%% ==========================================================================

function delete_stimulus_guis(h_fig),

stimuli_parameter = getappdata(h_fig,'stimuli_parameter');

for q = 1:length(stimuli_parameter.stimuli),
    if ~isempty(stimuli_parameter.stimuli(q).specs) & ...
            ishandle(stimuli_parameter.stimuli(q).h_stimulus_gui),
        stimulus_gui('delete_stimulus_gui', ...
            stimuli_parameter.stimuli(q).h_stimulus_gui, h_fig),
        stimuli_parameter.stimuli(q).h_stimulus_gui = [];
    end,
end,
setappdata(h_fig,'stimuli_parameter',stimuli_parameter),
%% ==========================================================================

function previous_state = Enable_uicontrols(h_fig, enable),

previous_state = get(findobj(h_fig,'Style','pushbutton'),'Enable');
previous_state = previous_state{1};

% disable changing of parmeter
children = get(h_fig,'Children');
for q = 1:length(children),
    if strcmp(get(children(q),'Type'),'uicontrol'),
        if ~strcmp(get(children(q),'Style'),'text'),
            set(children(q),'Enable', enable),
        end,
    end,
end,
set(findobj(h_fig,'Tag','pushbutton_play_firstA_thenB'),'Enable', 'on'),
drawnow,

% disable/enable changing Stimulus GUI's uicontrols
stimuli_parameter = getappdata(h_fig,'stimuli_parameter');
for q = 1:length(stimuli_parameter.stimuli),
    if ~isempty(stimuli_parameter.stimuli(q).specs),
        stimulus_gui('Enable_uicontrols',...
            stimuli_parameter.stimuli(q).h_stimulus_gui,enable),
    end,
end,
%% ==========================================================================

function stimuli = equalise_Nchannels(stimuli),

for q = 1:length(stimuli),
    sizes(q,:) = size(stimuli(q).waveform);
end,

if min(sizes(:,1)) == 1,
    % wave invalid (e.g. parameter changed, or wave generation failed)
    stimuli = 0;
elseif sizes(:,1) == 0 % wave generation in progress
    stimuli = [];
elseif max(sizes(:,2)) == 2,
    for q = 1:length(stimuli),
        if sizes(q,2) == 1, % duplicate mono channel
            stimuli(q).waveforms = ...
                [stimuli(q).waveforms, stimuli(q).waveforms];
        end,
    end,
end,
%% ==========================================================================

function interval_specs = Get_interval_specs(h_fig, interval),

% OBSOLETE. SEE SET_INTERVAL_SPECS

stimuli_parameter = Get_stimuli_parameter(h_fig);

q2 = 1;
for q = 1: length(stimuli_parameter.stimuli),
    if eval(['stimuli_parameter.stimuli(q).interval_' interval]),
        interval_specs(q2) = stimulus_gui('Get_specs', ...
            stimuli_parameter.stimuli(q).h_stimulus_gui);
        q2 = q2+1;
    end,
end,
%% ==========================================================================

function entry = get_menu_entry(h_menu),

menu_entries = get(h_menu,'String');
if isempty(menu_entries),
    entry = [];
    return,
else,
    entry = deblank(menu_entries(get(h_menu,'Value'),:));
    if iscell(entry)
        entry = entry{1};
    end,
end,

%% ==========================================================================

function stimuli_parameter = Get_stimuli_parameter(h_fig),

global SAMPLE_RATE,

stimuli_parameter = getappdata(h_fig,'stimuli_parameter');

for q = 1: length(stimuli_parameter.stimuli),
    if ~isempty(stimuli_parameter.stimuli(q).specs) && ...
            ~isempty(stimuli_parameter.stimuli(q).h_stimulus_gui)
        stimuli_parameter.stimuli(q).specs = stimulus_gui('Get_specs', ...
            stimuli_parameter.stimuli(q).h_stimulus_gui);
    end,
end,

stimuli_parameter.frozen = get(findobj(...
    h_fig,'Tag','checkbox_frozen'),'Value');
stimuli_parameter.isi = str2num(get(findobj(h_fig,'Tag','edit_isi'), ...
    'String'));
stimuli_parameter.sample_rate = SAMPLE_RATE;
%% ==========================================================================

function [waveform, waveform_info] = Get_waveform(button_tag, h_fig),

global SAMPLE_RATE,

previous_state = Enable_uicontrols(h_fig, 'off');

switch(button_tag),
    case 'A',
        interval_nos = 2;
    case 'B',
        interval_nos = 1;
    case {'AB','A_B'},
        interval_nos = [1 2];
    case {'ABA','A_AB'},
        interval_nos = [1 2 3];
    otherwise, error('Fcn ''Get_waveform'': unkown button_tag')
end, %switch

mixes = create_interval_mixes(h_fig, interval_nos);

if ~isstruct(mixes),
    waveform = mixes;
    waveform_info = mixes;
    Enable_uicontrols(h_fig, previous_state);
    return,
end,

isi_n =round(get(findobj(h_fig,'Tag','edit_isi'),'User')*SAMPLE_RATE/1000);
isi = zeros(isi_n, size(mixes(1).waveform,2));

switch(button_tag),
    case {'A', 'B'},
        waveform_info.button_string = button_tag;
        waveform = mixes(1).waveform;
        waveform_info.timing = mixes(1).n_samples/SAMPLE_RATE;
        waveform_info.level_struct = []; % to be done!
    case 'AB',
        if rand > 0.5,
            waveform_info.button_string = '{A B}';
            waveform = [mixes(2).waveform; isi; ...
                mixes(1).waveform];
            waveform_info.timing = [ ...
                mixes(2).n_samples isi_n mixes(1).n_samples ]/SAMPLE_RATE;
        else,
            waveform_info.button_string = '{B A}';
            waveform = [mixes(1).waveform; isi; ...
                mixes(2).waveform];
            waveform_info.timing = [...
                mixes(1).n_samples isi_n mixes(2).n_samples ]/SAMPLE_RATE;
        end,
        waveform_info.level_struct = []; % to be done!
    case 'ABA',
        no = rand;
        if no <= 1/3,
            waveform_info.button_string = '{B A A}';
            waveform = [mixes(1).waveform; isi; mixes(2).waveform; isi; ...
                mixes(3).waveform];
            waveform_info.timing = [mixes(1).n_samples isi_n ...
                mixes(2).n_samples isi_n mixes(3).n_samples ]/SAMPLE_RATE;
        elseif no > 2/3,
            waveform_info.button_string = '{A B A}';
            waveform = [mixes(2).waveform; isi; mixes(1).waveform; isi; ...
                mixes(3).waveform];
            waveform_info.timing = [mixes(2).n_samples isi_n ...
                mixes(1).n_samples isi_n mixes(3).n_samples ]/SAMPLE_RATE;
        else,
            waveform_info.button_string = '{A A B}';
            waveform = [mixes(2).waveform; isi; mixes(3).waveform; isi; ...
                mixes(1).waveform];
            waveform_info.timing = [mixes(2).n_samples isi_n ...
                mixes(3).n_samples isi_n mixes(1).n_samples ]/SAMPLE_RATE;
        end,
        waveform_info.level_struct = []; % to be done!
    case 'A_AB',
        if rand > 0.5,
            waveform_info.button_string = 'A {A B}';
            waveform = [mixes(2).waveform; isi; mixes(3).waveform; isi; ...
                mixes(1).waveform];
            waveform_info.timing = [mixes(2).n_samples isi_n ...
                mixes(3).n_samples isi_n mixes(1).n_samples ]/SAMPLE_RATE;
        else,
            waveform_info.button_string = 'A {B A}';
            waveform = [mixes(2).waveform; isi; mixes(1).waveform; isi; ...
                mixes(3).waveform];
            waveform_info.timing = [mixes(2).n_samples isi_n ...
                mixes(1).n_samples isi_n mixes(3).n_samples ]/SAMPLE_RATE;
        end,
        waveform_info.level_struct = []; % to be done!
    case 'A_AB_A',
        if rand > 0.5,
            waveform_info.button_string = 'A {A B} A';
            waveform = [mixes(2).waveform; isi; mixes(3).waveform; isi; ...
                mixes(1).waveform; isi; mixes(2).waveform];
            waveform_info.timing = [mixes(2).n_samples isi_n ...
                mixes(3).n_samples isi_n mixes(1).n_samples isi_n ...
                mixes(2).n_samples]/SAMPLE_RATE;
        else,
            waveform_info.button_string = 'A {B A} A';
            waveform = [mixes(2).waveform; isi; mixes(1).waveform; isi; ...
                mixes(3).waveform; isi; mixes(2).waveform];
            waveform_info.timing = [mixes(2).n_samples isi_n ...
                mixes(1).n_samples isi_n mixes(3).n_samplesisi_n ...
                mixes(2).n_samples]/SAMPLE_RATE;
        end,
        waveform_info.level_struct = []; % to be done!
end, %switch
Enable_uicontrols(h_fig, previous_state);
%% ==========================================================================

function h_fig = initialise(h_fig_parent, position, filename),

global SAMPLE_RATE,

% if parameter file exist overwrite defaults
if exist('filename','var'),
    try,
        load(filename,'stimuli_parameter'),
        SAMPLE_RATE = stimuli_parameter.sample_rate;
    catch,
        warning(lasterr);
    end,
else,
    % defaults (to open two stimuli GUI with two empty stimuli)
    stimuli_parameter.stimuli(1).specs = [];
    stimuli_parameter.stimuli(1).interval_A = 0;
    stimuli_parameter.stimuli(1).interval_B = 0;
    stimuli_parameter.stimuli(1).h_stimulus_gui = [];
    stimuli_parameter.isi = 100;
    stimuli_parameter.frozen = 0;
    SAMPLE_RATE = 48000;
end,

if ~exist('h_fig_parent','var'),
    h_fig_parent = [];
    rand('state',sum(100*clock)),
end,

% Open stimuli_gui figure. Call main function to pass var 'vers'
h_fig = stimuli_gui('open_stimuli_gui',stimuli_parameter);
Change_parameter(h_fig, 0, 'ISI [ms]', stimuli_parameter.isi);
set(findobj(h_fig,'Tag','checkbox_frozen'),'Value', stimuli_parameter.frozen);

% set position
current_position = get(h_fig,'Position');
if ~exist('position','var')|| isempty(position),
    position = [0,0];
elseif length(position) ~= 2, % position does not contain y-dimension of figure
    h = errordlg('Position must be a 2 element vector! No dimensions!', ...
        'stimulus_gui.m - initialise');
    uiwait(h),
    return,
end,
set(h_fig,'Position',[position current_position(3:4)]),

stimuli_parameter = open_stimulus_guis(h_fig, stimuli_parameter);

% initialise ALL appdata (list complete)
waveforms.dirty_h_stimulus_list = [];
waveforms.stimuli = [];
setappdata(h_fig,'waveforms',waveforms),
setappdata(h_fig, 'h_fig_parent', h_fig_parent)
setappdata(h_fig, 'position', []),
setappdata(h_fig, 'stimuli_parameter',stimuli_parameter),

Refresh_waveforms(h_fig),
%% ==========================================================================

function Level_meter(h_fig, action, interval, level_struct),

stimuli_parameter = getappdata(h_fig,'stimuli_parameter');

switch action,
    case 'full_reset', % resets includes red overload LEDs
        % stimuli GUI
        set(findobj(h_fig,'Type','rectangle'), 'FaceColor',[0 0 0]),

        % stimulus GUIs
        for q = 1: length(stimuli_parameter.stimuli),
            stimulus_gui('Level_meter',stimuli_parameter.stimuli(q).h_stimulus_gui, ...
                'full_reset'),
        end,

    case 'reset', % keeps the state of red overload LEDs
        % stimuli GUI
        overload = findobj(h_fig,'FaceColor', [1 0 0]);
        leds = findobj(h_fig,'Type','rectangle');
        l_leds = length(leds);
        if ~exist('interval','var') | strcmp(interval,'A'),
            set(leds(l_leds/2+1:l_leds), 'FaceColor',[0 0 0]),
        end,
        if ~exist('interval','var') | strcmp(interval,'B'),
            set(leds(1:l_leds/2), 'FaceColor',[0 0 0]),
        end,
        set(overload, 'FaceColor', [1 0 0]),

        % stimulus GUIs
        for q = 1: length(stimuli_parameter.stimuli),
            if (~exist('interval','var') | ...
                    eval(['stimuli_parameter.stimuli(q).interval_' interval])) ...
                    & ~isempty(stimuli_parameter.stimuli(q).specs),
                stimulus_gui('Level_meter',stimuli_parameter.stimuli(q).h_stimulus_gui,'reset'),
            end,
        end,

    case 'update',
        % updates the level meter. never resets the uppermost red overload LEDs.
        return,
        % Stimulus Guis
        for q = 1: length(stimuli_parameter.stimuli),
            %             if ~exist('interval','var') | eval( ...
            %                     ['stimuli_parameter.stimuli(q).interval_' interval]),
            if ~isempty(stimuli_parameter.stimuli(q).specs),
                stimulus_gui('Level_meter',...
                    stimuli_parameter.stimuli(q).h_stimulus_gui, ...
                    'update', [-5 -5; -10 -10]), % !!! to be done!!!
            end,
            %             end,
        end,

        % !!! to be done!!!

        %         % stimuli GUI
        %         if ~exist('interval','var') | strcmp(interval,'A'),
        %             if max(waveforms{2}) > max(waveforms{3}),
        %                 waveform = waveforms{2};
        %             else,
        %                 waveform = waveforms{3};
        %             end,
        %             if size(waveform,1)>1,
        %                 n_LEDs = 4*(length(stimuli_parameter.stimuli)+4);
        %                 max_levels = max(waveform);
        %                 rms_levels = sqrt(mean(waveform.^2));
        %                 set_level_meter(h_fig, 'A', n_LEDs, max_levels, rms_levels, 1),
        %             end,
        %         end,
        %         if ~exist('interval','var') | strcmp(interval,'B'),
        %             waveform = waveforms{1};
        %             if size(waveform,1)>1,
        %                 n_LEDs = 4*(length(stimuli_parameter.stimuli)+4);
        %                 max_levels = max(waveform);
        %                 rms_levels = sqrt(mean(waveform.^2));
        %                 set_level_meter(h_fig, 'B', n_LEDs, max_levels, rms_levels, 1),
        %             end,
        %         end,
    otherwise, error('Fcn ''levelmeter'': Unknown action')
end,
%% ==========================================================================

function error = load_config(filename),

global SAMPLE_RATE, %gets loaded

if ~exist('filename','var'),
    stimuli_path = fileparts(which('stimuli_gui'));

    [filename,pathname] = uigetfile([stimuli_path ...
        '*.mat'],'Load configuration:');
end,

if ~isempty(filename) & filename,
    try
        load([pathname,filename],'stimuli_parameter'),
        error = 0;
    catch,
        error = 1;
        h = errordlg(['Can''t find configuration file: ' pathname filename], ...
            'stimulus_gui.m - load_config');
        uiwait(h),
        return,
    end,
    Set_stimuli_parameter(gcbf, stimuli_parameter),
end,
%% ==========================================================================

function files = my_dir(path, suffix, begin)
% files = MY_DIR(path,begin,suffix) returns files and directories of
% directory $path. They will be sorted by name.
% If begin is given MY_DIR returns only files/directories beginnig with
% $begin.
% If $suffix is given MY_DIR returns only files with ending of $suffix
% $suffix might be '*'.

if(nargin < 1),
    path ='.';
end,
if(nargin < 2),
    suffix = '*';
end,
if(nargin < 3),
    begin = [];
end,
tmp = dir([path, '\*.', suffix]);
str = ''; m=0;

% to avoid '.' and '..'
if (strcmp(suffix,'*')),
    start = 3;
else,
    start = 1;
end,

% get the [begin *]-files only
for (n=start:length(tmp)),
    tmp2 = getfield(tmp(n),'name');
    if(length(tmp2)>=length(begin)),
        if(isempty(begin) | strcmp(tmp2(1:length(begin)),begin)) ,
            m=m+1;
            tmp3{m} = tmp2;
        end,
    end,
end,
% sort file names (only possible as matrix)(tmp3 is cell array!)
if (exist('tmp3')),
    for (n=1:length(tmp3)),
        files(n,1:length(tmp3{n})) = tmp3{n};
    end,
    files = sortrows(files);
else,
    files = [];
end,
%% ==========================================================================

function new_stimulus,

stimuli_parameter = Get_stimuli_parameter(gcbf);
q = length(stimuli_parameter.stimuli)+1;
if isempty(stimuli_parameter.stimuli(q-1).specs),
    return,
end,
stimuli_parameter.stimuli(q).specs = stimuli_parameter.stimuli(q-1).specs;
stimuli_parameter.stimuli(q).interval_A = 0;
stimuli_parameter.stimuli(q).interval_B = 0;
stimuli_parameter.stimuli(q).h_stimulus_gui = [];
reopen_GUIs(gcbf,stimuli_parameter),
%% ==========================================================================

function h_fig = open_stimuli_gui(vers,stimuli_parameter,h_fig),

if exist('h_fig'),
    figure(h_fig),
    delete(get(h_fig, 'Children')),
else
    h_fig = figure;
end,


set(h_fig,'Color',get(0,'defaultUicontrolBackgroundColor'), ...
    'Backingstore','off', ...
    'ButtonDownFcn', 'stimuli_gui ButtonDownFnc_figure', ...
    'CloseRequestFcn', ...
    ['stimuli_gui(''delete_stimuli_gui'',' num2str(h_fig.Number) ')'], ...
    'HandleVisibility','on', ...
    'MenuBar', 'none', ...
    'Name', ['Stimuli GUI ' vers], ...
    'NumberTitle','off', ...
    'Position',[0 0 212 98+length(stimuli_parameter.stimuli)*24], ...
    'Resize','off', ...
    'Tag','stimuli_gui_fig', ...
    'windowButtonMotionFcn','stimuli_gui windowButtonMotionFcn', ...
    'windowButtonUpFcn','stimuli_gui windowButtonUpFcn'),

uicontrol('Parent',h_fig, ...
    'BackgroundColor', [1 1 1], ...
    'Callback', 'stimuli_gui(''Change_parameter'',gcbf);',...
    'FontUnits','pixel',...
    'FontSize', 10, ...
    'Units', 'pixel', ...
    'Position',[119 26 50 18], ...
    'Style','edit', ...
    'Tag','edit_isi');

uicontrol('Parent',h_fig, ...
    'ButtonDownFcn','stimuli_gui ButtonDownFnc_text', ...
    'Enable','inactive', ...
    'FontUnits','pixel',...
    'FontSize', 10, ...
    'Units', 'pixel', ...
    'Position',[78 24 40 16], ...
    'String','ISI [ms]:', ...
    'Style','text', ...
    'Tag','text_isi', ...
    'Userdata', []);

uicontrol('Parent',h_fig, ...
    'callback','ans = stimuli_gui(''callback_save_A'');', ...
    'FontUnits','pixel',...
    'FontSize', 10, ...
    'Units', 'pixel', ...
    'Position',[7 50 26 18], ...
    'String','save', ...
    'Style','pushbutton', ...
    'Tag','pushbutton_save_A');

uicontrol('Parent',h_fig, ...
    'callback','stimuli_gui save_config', ...
    'FontUnits','pixel',...
    'FontSize', 10, ...
    'Units', 'pixel', ...
    'Position',[44 50 28 18], ...
    'String','Save', ...
    'Style','pushbutton', ...
    'Tag','pushbutton_save_config');

uicontrol('Parent',h_fig, ...
    'callback','stimuli_gui load_config', ...
    'FontUnits','pixel',...
    'FontSize', 10, ...
    'Units', 'pixel', ...
    'Position',[77 50 38 18], ...
    'String','Load', ...
    'Style','pushbutton', ...
    'Tag','pushbutton_load_config');

uicontrol('Parent',h_fig, ...
    'callback','stimuli_gui new_stimulus', ...
    'FontUnits','pixel',...
    'FontSize', 10, ...
    'Units', 'pixel', ...
    'Position',[119 50 50 18], ...
    'String','New Stim.', ...
    'Style','pushbutton', ...
    'Tag','pushbutton_new_stimulus');

uicontrol('Parent',h_fig, ...
    'callback','ans = stimuli_gui(''callback_save_B'');', ...
    'FontUnits','pixel',...
    'FontSize', 10, ...
    'Units', 'pixel', ...
    'Position',[181 50 26 18], ...
    'String','save', ...
    'Style','pushbutton', ...
    'Tag','pushbutton_save_B');

uicontrol('Parent',h_fig, ...
    'callback','ans = stimuli_gui(''callback_plot_A'');', ...
    'FontUnits','pixel',...
    'FontSize', 10, ...
    'Units', 'pixel', ...
    'Position',[7 26 26 18], ...
    'String','plot', ...
    'Style','pushbutton', ...
    'Tag','pushbutton_plot_A');

uicontrol('Parent',h_fig, ...
    'Callback', 'stimuli_gui(''Refresh_waveforms'',gcbf),', ...
    'FontUnits','pixel',...
    'FontSize', 9, ...
    'Units', 'pixel', ...
    'Position',[44 26 35 18], ...
    'String','frz', ...
    'Style','checkbox', ...
    'Tag','checkbox_frozen');

uicontrol('Parent',h_fig, ...
    'callback','ans = stimuli_gui(''callback_plot_B'');', ...
    'FontUnits','pixel',...
    'FontSize', 10, ...
    'Units', 'pixel', ...
    'Position',[181 26 26 18], ...
    'String','plot', ...
    'Style','pushbutton', ...
    'Tag','pushbutton_plot_B');

uicontrol('Parent',h_fig, ...
    'callback','ans = stimuli_gui(''Play'',''A'',gcbf);', ...
    'FontUnits','pixel',...
    'FontSize', 10, ...
    'Units', 'pixel', ...
    'Position',[7 2 26 18], ...
    'String','A', ...
    'Style','pushbutton', ...
    'Tag','pushbutton_play_A');

uicontrol('Parent',h_fig, ...
    'callback','ans = stimuli_gui(''Play'',''AB'',gcbf);', ...
    'FontUnits','pixel',...
    'FontSize', 9, ...
    'Units', 'pixel', ...
    'Position',[44 2 28 18], ...
    'String','{A B}', ...
    'Style','pushbutton', ...
    'Tag','pushbutton_play_AB');

uicontrol('Parent',h_fig, ...
    'callback','ans = stimuli_gui(''Play'',''ABA'',gcbf);', ...
    'FontUnits','pixel',...
    'FontSize', 9, ...
    'Units', 'pixel', ...
    'Position',[77 2 38 18], ...
    'String','{A B A}', ...
    'Style','pushbutton', ...
    'Tag','pushbutton_play_ABA');

uicontrol('Parent',h_fig, ...
    'callback','ans = stimuli_gui(''Play'',''A_AB'',gcbf);', ...
    'FontUnits','pixel',...
    'FontSize', 9, ...
    'Units', 'pixel', ...
    'Position',[119 2 50 18], ...
    'String','A {A B}', ...
    'Style','pushbutton', ...
    'Tag','pushbutton_play_A_AB');

uicontrol('Parent',h_fig, ...
    'callback','ans = stimuli_gui(''Play'',''B'',gcbf);', ...
    'FontUnits','pixel',...
    'FontSize', 10, ...
    'Units', 'pixel', ...
    'Position',[181 2 26 18], ...
    'String','B', ...
    'Style','pushbutton', ...
    'Tag','pushbutton_play_B');

for q = 1:length(stimuli_parameter.stimuli),
    uicontrol('Parent',h_fig, ...
        'callback','stimuli_gui callback_checkbox_A', ...
        'Position',[13 (length(stimuli_parameter.stimuli)-q)*24+103 18 18], ...
        'Style','checkbox', ...
        'Tag',['checkbox_A_', num2str(q)]);
    uicontrol('Parent',h_fig, ...
        'BackgroundColor',[1 1 1], ...
        'callback','stimuli_gui callback_popupmenu', ...
        'Position',[45 (length(stimuli_parameter.stimuli)-q)*24+103 126 18], ...
        'String','Can''t find any stimuli', ...
        'Style','popupmenu', ...
        'Tag',['popupmenu_', num2str(q)]);
    uicontrol('Parent',h_fig, ...
        'callback','stimuli_gui callback_checkbox_B', ...
        'Position',[187 (length(stimuli_parameter.stimuli)-q)*24+103 18 18], ...
        'Style','checkbox', ...
        'Tag',['checkbox_B_', num2str(q)]);
end,

% create Level_meter
tags = {'left_LED_A_','right_LED_A_','left_LED_B_','right_LED_B_'};
pos_LEDs = [1 34 175 208];

for q2 = 1:4;
    axes('Units','pixels', ...
        'Position', [pos_LEDs(q2), 1, 4, (length(stimuli_parameter.stimuli)+4)*24], ...
        'XTick',[],'YTick',[], ...
        'XColor',[.5 .5 .5],'YColor',[.5 .5 .5], ...
        'XLim',[-1 -.7],'YLim',[0 (length(stimuli_parameter.stimuli)+4)*24],...
        'Color',[.5 .5 .5], ...
        'Tag','Level_meter');

    for q = 0:4*(4+length(stimuli_parameter.stimuli))-1,
        led = rectangle('Position', [-1 1+6*q .3 5], ...
            'EdgeColor', [.4 .4 .4], ...
            'FaceColor', [0 0 0]);
        set(led,'Tag',[tags{q2}, num2str(q+1)])
    end,
end,


% fill popupmenus
stimuli_path = fileparts(which('stimuli_gui'));
files = my_dir(stimuli_path,'m');
[n m] = size(files);
idx = ceil(findstr(reshape(files',1,n*m),'stimuli_gui.m')/m);
files = [files(1:idx-1,:);files(idx+1:n,:)];
idx = ceil(findstr(reshape(files',1,(n-1)*m),'stimulus_gui.m')/m);
files = [files(1:idx-1,:);files(idx+1:n-1,:)];
n_files = size(files,1);
if isempty(files), stimuli_parameter.stimuli=[]; end,
for q = 1:length(stimuli_parameter.stimuli),
    popupmenu_cell = [];
    h_popupmenu = findobj(h_fig,'Tag',['popupmenu_',num2str(q)]);

    popupmenu_cell{1} = 'Remove Stimulus';
    for  q2 = 1:length(stimuli_parameter.stimuli),
        if q2 ~= q,
            popupmenu_cell{length(popupmenu_cell)+1} = ...
                ['Like Stimulus ' num2str(q2)];
        end,
    end,
    numbering = char(ones(n_files,1)*[num2str(q) ': ']);
    popupmenu_cell = [popupmenu_cell ...
        deblank(num2cell([numbering files],2))'];

    set(h_popupmenu,'String',popupmenu_cell),

    if ~isempty(stimuli_parameter.stimuli(q).specs)
        [match value] = max(strcmp(popupmenu_cell, [num2str(q) ': ' ...
            strtok(stimuli_parameter.stimuli(q).specs.name),'.m']));
        if match,
            drawnow,
            set(h_popupmenu, 'Value', value);
        else,
            stimuli_parameter.stimuli(q).specs = [];
            h = errordlg(['Can''t find stimuli: ' ...
                strtok(stimuli_parameter.stimuli(q).specs.name),'.m'], ...
                ['Error: Stimuli GUI ' vers]);
            uiwait(h),
        end,
        set(findobj(h_fig,'Tag',['checkbox_A_' num2str(q)]), ...
            'Value', stimuli_parameter.stimuli(q).interval_A),
        set(findobj(h_fig,'Tag',['checkbox_B_' num2str(q)]), ...
            'Value', stimuli_parameter.stimuli(q).interval_B),
    end,
end,

setappdata(h_fig,'stimuli_parameter',stimuli_parameter),
%% ==========================================================================

function stimuli_parameter = open_stimulus_guis(h_fig, stimuli_parameter),

for q = 1:length(stimuli_parameter.stimuli),
    if ~isempty(stimuli_parameter.stimuli(q).specs) & ...
            isempty(stimuli_parameter.stimuli(q).h_stimulus_gui), % AND not already open
        eval(['stimuli_parameter.stimuli(q).h_stimulus_gui = ' ...
            strtok(stimuli_parameter.stimuli(q).specs.name) ...
            '([],stimuli_parameter.stimuli(q).specs,h_fig);']),
        set(stimuli_parameter.stimuli(q).h_stimulus_gui, ...
            'Name', [num2str(q) ': ' ...
            get(stimuli_parameter.stimuli(q).h_stimulus_gui,'Name')]),
    end,
end,
setappdata(h_fig,'stimuli_parameter',stimuli_parameter),
Update_stimulus_gui_position(h_fig),
drawnow,
%% ==========================================================================

function [waveform, waveform_info] = Play(button_tag, h_fig),

Level_meter(h_fig, 'reset'), % visually confirmes play button press!
[waveform, waveform_info] = Get_waveform(button_tag,h_fig);
if length(waveform) > 1,
    Level_meter(h_fig, 'update', waveform_info.level_struct),
    play_waveform(waveform);
    set(findobj(h_fig,'Tag',['pushbutton_play_' ...
        button_tag]),'String', waveform_info.button_string, ...
        'BackgroundColor', get(0,'defaultUicontrolBackgroundColor')),
    switch(button_tag),
        case 'A',
            interval_nos = 2;
        case 'B',
            interval_nos = 1;
        case {'AB'},
            interval_nos = [1 2];
        case {'ABA','A_AB'},
            interval_nos = [1 2 3];
        otherwise, error('Fcn ''Play'': unkown button_tag')
    end, %switch
    Refresh_waveforms(h_fig, interval_nos),
else,
    %    error('Stimuli_gui Fcn: Play(): Getwaveform returned invalid waveform!'),
end,
%% ==========================================================================

function play_waveform(waveform),

global SAMPLE_RATE,

try,
    if size(waveform,1) > 1,
        if strcmp(computer, 'PCWIN'),
            wavplay(waveform,SAMPLE_RATE,'async');
        else,
            sound(waveform,SAMPLE_RATE,16);
        end,
    end,
catch
    warning(lasterr);
end,
%% ==========================================================================

function waveform = plot_interval(h_fig, interval),

global SAMPLE_RATE,

if ~exist('h_fig','var'), h_fig = gcbf; end;

waveform = Get_waveform(interval, h_fig);
if size(waveform,1) > 1,
    [m n]=size(waveform); m=2*floor(m/2); if m < 3, return,  end,

    % if ~exist('PLOTS','dir'),
    %     % create context menu of what plots are available in directory 'PLOTS'
    %     return,
    % end,

    scrsz = get(0,'ScreenSize');
    position = [0,round(scrsz(4)*3/4)-18,round(scrsz(3)/3.6),...
        round(scrsz(4)/4-40)];
    h_plot_fig = findobj('Position',position);
    if (isempty(h_plot_fig)),
        h_plot_fig = figure('Position',position);
    end,
    figure(h_plot_fig); clf,
    set(h_plot_fig,'Name','wave'),

    [m n]=size(waveform);
    x=1000*linspace(0,m/SAMPLE_RATE,m)';
    plot(x,waveform(:,1)), hold on,
    plot(x,waveform(:,n),'r'), grid on, hold off, title(''), zoom on,
else,
    uiwait(warndlg('Waveform is empty!', 'Stimuli_GUI - plot_waveform')),
end,
if strcmp(interval,'B'),
    Refresh_waveforms(h_fig, 1)
else,
    Refresh_waveforms(h_fig, 2)
end,

%% ==========================================================================

function Refresh_waveforms(h_fig, interval_nos),

frozen = get(findobj(h_fig,'Tag','checkbox_frozen'),'Value');
waveforms = getappdata(h_fig,'waveforms');
stimuli_parameter = getappdata(h_fig,'stimuli_parameter');

previous_state = Enable_uicontrols(h_fig, 'off');

if  ~exist('interval_nos'),
    interval_nos = 1:3;
end,


if max(interval_nos == 1),
    inteval_B = 1;
else
    inteval_B = 0;
end,
if max(interval_nos > 1),
    inteval_A = 1;
else
    inteval_A = 0;
end,

list_h_stimulus_gui_dirty = [];

if frozen,
    for q = 1:length(stimuli_parameter.stimuli);
        if (interval_A & stimuli_parameter.stimuli(q).interval_A | ...
                interval_B & stimuli_parameter.stimuli(q).interval_B)...
                & ~isempty(stimuli_parameter.stimuli(q).specs),

            stimulus_gui('Refresh_waveform', ...
                stimuli_parameter.stimuli(q).h_stimulus_gui),
            % reset dirty_flag in stimulus_gui
            stimulus_gui('Check_and_reset_dirty_flag', ...
                stimuli_parameter.stimuli(q).h_stimulus_gui, h_fig);
            list_h_stimulus_gui_dirty = [list_h_stimulus_gui_dirty ...
                stimuli_parameter.stimuli(q).h_stimulus_gui];
        end,
    end,
end,


for q3 = interval_nos,
    q2 = 0; stimuli = [];

    for q = 1:length(stimuli_parameter.stimuli),
        if q3 == 1,
            interval = 'B';
        else,
            interval = 'A';
        end,
        if eval(['stimuli_parameter.stimuli(q).interval_' interval])...
                & ~isempty(stimuli_parameter.stimuli(q).specs),
            q2 = q2+1;
            stimuli(q2).h_stimulus_gui = ...
                stimuli_parameter.stimuli(q).h_stimulus_gui;
            if ~frozen, % if not refreshed in previous loop (frozen)
                stimulus_gui('Refresh_waveform', ...
                    stimuli_parameter.stimuli(q).h_stimulus_gui),
                % reset dirty_flag in stimulus_gui
                stimulus_gui('Check_and_reset_dirty_flag', ...
                    stimuli_parameter.stimuli(q).h_stimulus_gui, h_fig);
                list_h_stimulus_gui_dirty = [list_h_stimulus_gui_dirty ...
                    stimuli_parameter.stimuli(q).h_stimulus_gui];
            end,
            [stimuli(q2).waveform, stimuli(q2).levels] = ...
                stimulus_gui('Get_waveform', ...
                stimuli_parameter.stimuli(q).h_stimulus_gui);
        end,
    end,

    waveforms(q3).stimuli = stimuli;
end,

% set dirty flags in all intevals which have not been refreshed
waveforms = set_dirty_list(h_fig, waveforms, interval_nos, ...
    list_h_stimulus_gui_dirty);

setappdata(h_fig,'waveforms',waveforms),
Enable_uicontrols(h_fig, previous_state);
%% ==========================================================================

function reopen_GUIs(h_fig, stimuli_parameter),

position = get(h_fig,'Position');
delete_stimulus_guis(h_fig),
for q = 1:length(stimuli_parameter.stimuli),
    stimuli_parameter.stimuli(q).h_stimulus_gui = [];
end,
h_fig = stimuli_gui('open_stimuli_gui',stimuli_parameter,h_fig);
curr_position = get(h_fig,'Position');
set(h_fig,'Position', [position(1:2) curr_position(3:4)]),
Change_parameter(h_fig, 0, 'ISI [ms]', stimuli_parameter.isi);
set(findobj(h_fig,'Tag','checkbox_frozen'),'Value', stimuli_parameter.frozen);
stimuli_parameter = open_stimulus_guis(h_fig, stimuli_parameter);
setappdata(h_fig, 'stimuli_parameter',stimuli_parameter),
Refresh_waveforms(h_fig)
%% ==========================================================================

function save_config,

stimuli_path = fileparts(which('stimuli_gui'));

[filename,pathname] = uiputfile([stimuli_path, ...
    'config.mat'], ...
    'Save configuration as:');

if ~isempty(filename) & filename,
    stimuli_parameter = Get_stimuli_parameter(gcbf);
    for q = 1:length(stimuli_parameter.stimuli),
        stimuli_parameter.stimuli(q).h_stimulus_gui = [];
    end,
    stimuli_parameter.time = datestr(now);
    % '-append' because the file may contain parameter of parent program
    if exist([pathname,filename],'file')
        save([pathname,filename],'stimuli_parameter','-append'),
    else
        save([pathname,filename],'stimuli_parameter'),
    end,
end,
%% ==========================================================================

function waveform = Save_waveform(h_fig, interval, filename),

global SAMPLE_RATE,

if ~exist('filename','var'),
    [filename,pathname] = uiputfile(['c:\Windows\Temp\temp.wav'],'Save wave as:');
    [filename suffix]=strtok(filename,'.');
    filename =[pathname,filename];
end,

waveform = Get_waveform(interval, h_fig);
if size(waveform,1) > 1,
    if strcmp(suffix,'.wav'),
        wavwrite(waveform,SAMPLE_RATE,16,filename),
    elseif strcmp(suffix,'.mat'),
        save( filename,'waveform'),
    end,
else,
    uiwait(warndlg('Waveform is empty!', 'Stimuli_GUI - Save_waveform')),
end,

if strcmp(interval,'B'),
    Refresh_waveforms(h_fig, 1),
else,
    Refresh_waveforms(h_fig, 2),
end,

%% ==========================================================================

function variable = Select_edit_uicontrol(h_fig, variable),

stimuli_parameter = getappdata(h_fig,'stimuli_parameter');

if exist('variable'), % variable has only valid numbers. Return handles.
    if variable.stimulus.number,
        variable.stimulus.handle = ...
            stimuli_parameter.stimuli(variable.stimulus.number).h_stimulus_gui;
        variable.parameter.handle = findobj(variable.stimulus.handle,'Tag', ...
            ['edit_' num2str(variable.parameter.number)]);
    else,
        variable.stimulus.handle = h_fig;
        if variable.parameter.number == 1
            variable.parameter.handle = findobj(h_fig,'Tag','edit_isi');
        end
    end,
    return,
end,

set(h_fig,'CurrentObject',0),
for q = 1:length(stimuli_parameter.stimuli),
    set(stimuli_parameter.stimuli(q).h_stimulus_gui,'CurrentObject',0),
end,

h_msg = warndlg('Select parameter by mouse click on text label!','Stimuli GUI');
delete(findobj(h_msg,'Style','pushbutton')),
set(h_msg,'Tag','stimulus_gui_select_edit_uicontrol_fig'),
uiwait(h_msg),

variable.stimulus.handle = h_fig;
variable.stimulus.number = 0;
variable.stimulus.interval = 'isi';
variable.parameter.handle = findobj(h_fig,'Tag','edit_isi');
variable.parameter.number = 1;
variable.parameter.name = 'isi';


for q = 1:length(stimuli_parameter.stimuli),
    current_object = get(stimuli_parameter.stimuli(q).h_stimulus_gui,'CurrentObject');
    if current_object,
        variable.stimulus.handle = stimuli_parameter.stimuli(q).h_stimulus_gui;
        variable.stimulus.number = q;
        if stimuli_parameter.stimuli(q).interval_A,
            variable.stimulus.interval = 'A';
        elseif stimuli_parameter.stimuli(q).interval_B,
            variable.stimulus.interval = 'B';
        else,
            variable.stimulus.interval = '';
        end,
        variable.parameter.handle = current_object;
        tag = get(variable.parameter.handle,'Tag');
        variable.parameter.number = str2num(tag(length(tag)));
        variable.parameter.name = ...
            stimuli_parameter.stimuli(q).specs.param(variable.parameter.number).name;
        break,
    end,
end,
%% ==========================================================================

function waveforms = set_dirty_list(h_fig, waveforms, ...
    clean_interval_nos, list_h_stimulus_gui_dirty),
% set dirty list in all intevals which have not been refreshed
% set dirty list empty in all intevals which have been refreshed

for interval_no = 1:3,
    if interval_no ~= clean_interval_nos,
        for q = 1:length(list_h_stimulus_gui_dirty),
            % check for double entries
            if ~max(list_h_stimulus_gui_dirty(q) == ...
                    waveforms(interval_no).dirty_h_stimulus_list),
                waveforms(interval_no).dirty_h_stimulus_list = ...
                    [waveforms(interval_no).dirty_h_stimulus_list, ...
                    list_h_stimulus_gui_dirty(q)];
            end,
        end,
    else,
        waveforms(interval_no).dirty_h_stimulus_list = [];
    end,
end,
%% ==========================================================================

function set_level_meter(h_fig, interval, n_LEDs, max_levels, rms_levels, stepsize_dB),

if length(max_levels) == 1, % mono stimulus
    max_levels = [max_levels max_levels];
    rms_levels = [rms_levels rms_levels];
end,

leds = zeros(n_LEDs,size([max_levels max_levels; rms_levels rms_levels],2));
try,
    leds = hist(20*log10([max_levels max_levels; rms_levels rms_levels]), ...
        [-(n_LEDs-1)*stepsize_dB:stepsize_dB:0]+0.5);
end,
leds( flipud(cumsum(flipud(leds(:,1)))) == max(cumsum(leds(:,1))),1) = 1;
leds( flipud(cumsum(flipud(leds(:,2)))) == max(cumsum(leds(:,2))),2) = 1;
leds(1,1:2)  = max_levels ~= 0; % lowest level LED on only if max_levels > 0

% set LEDs
for q = 1:n_LEDs-1,

    if leds(q,1),
        color = [0 1 0]; % green
        if q > n_LEDs-7, color = [1 1 0]; end, % yellow
    else,
        color = [0 0 0];
    end,
    set(findobj(h_fig,'Tag',['left_LED_' interval '_' num2str(q)]), ...
        'FaceColor', color),

    if leds(q,2),
        color = [0 1 0]; % green
        if q > n_LEDs-7, color = [1 1 0]; end, % yellow
    else,
        color = [0 0 0];
    end,

    set(findobj(h_fig,'Tag',['right_LED_' interval '_' num2str(q)]), ...
        'FaceColor', color),
end,

% set uppermost red overload LEDs if max_level > 1. Never overrides an overload.
if leds(n_LEDs,1),
    set(findobj(h_fig,'Tag',['left_LED_' interval '_' num2str(n_LEDs)]), ...
        'FaceColor',[1 0 0]),
end,
if leds(n_LEDs,2),
    set(findobj(h_fig,'Tag',['right_LED_' interval '_' num2str(n_LEDs)]), ...
        'FaceColor',[1 0 0]),
end,

%% ==========================================================================
function Set_stimuli_parameter(h_fig, stimuli_parameter),

global SAMPLE_RATE,

SAMPLE_RATE = stimuli_parameter.sample_rate;
enable = get(findobj(h_fig,'Style','pushbutton'),'Enable');
reopen_GUIs(h_fig, stimuli_parameter),
Enable_uicontrols(h_fig,enable{1});
%% ==========================================================================

function Update_stimulus_gui_position(h_fig),

stimuli_parameter = getappdata(h_fig,'stimuli_parameter');

if isfield(stimuli_parameter,'stimuli'),
    units = get(0,'Units');
    set(0,'Units', 'pixels')
    screen_size = get(0,'ScreenSize');
    set(0,'Units', units)

    pos_fig = get(h_fig,'Position');
    for q = 1:length(stimuli_parameter.stimuli),
        if ~isempty(stimuli_parameter.stimuli(q).specs)& ...
                ~isempty(stimuli_parameter.stimuli(q).h_stimulus_gui),
            pos_curr(q,:) = get( ...
                stimuli_parameter.stimuli(q).h_stimulus_gui,'Position');
        else,
            pos_curr(q,:) = [0 0 150 0];
        end,
    end,

    % determine side of Stimulus GUIs relative to Stimuli GUI
    if pos_fig(1) > screen_size(3)/2 - pos_fig(3)/2,
        offset = -sum(pos_curr(:,3)+3);
    else,
        offset = pos_fig(3)+3;
    end,

    for q = 1:length(stimuli_parameter.stimuli),
        if ~isempty(stimuli_parameter.stimuli(q).specs)& ...
                ~isempty(stimuli_parameter.stimuli(q).h_stimulus_gui),
            set(stimuli_parameter.stimuli(q).h_stimulus_gui,'Position', ...
                [offset+pos_fig(1)+sum(pos_curr(1:q-1,3)+3), ...
                pos_fig(2), ...
                pos_curr(q,3),pos_curr(q,4)]),
            figure(stimuli_parameter.stimuli(q).h_stimulus_gui),
        end,
    end,
    figure(h_fig)
end,
%% ==========================================================================

function windowButtonMotionFcn,
% when moving the Stimuli GUI:

current_position  = getappdata(gcbf,'position');
if ~isempty(current_position),
    curr_pos = get(gcbf,'Position');
    units = get(0,'Units');
    set(0,'Units', 'pixels')
    pointer_location = get(0,'PointerLocation');
    set(0,'Units', units)

    set(gcbf,'Position',[pointer_location - current_position ,curr_pos(3:4)]),
end,
%% ==========================================================================

function windowButtonUpFcn,
% after moving the Stimuli GUI: reset storage place of current position
% and get Stimulus GUIs next to Stimuli GUI

h_fig = gcbf;
setappdata(h_fig,'position',[]),
if gcbo == h_fig,
    Update_stimulus_gui_position(h_fig),
end,
%% ==========================================================================