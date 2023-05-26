function varargout = stimulus_gui(varargin),
% Graphical user interface for a single stimulus generator
%
% External usage:
% Only functions starting with a capital letter are thought to be useful
% to call externally from an host application.
% Change stimulus parameter by calling Change_parameter() or Set_specs()
% in combination with Get_specs(). Get_waveform() returns the precomputed
% waveform. See e.g. stimuli_gui.m.
%
% After an external stimulus paramter change using:
%  >> stimulus_gui('Change_parameter', h_stimulus_gui, parameter, new_value),
% one must also call:
%  >> stimulus_gui('Refresh_waveform', h_stimulus_gui),
% to give the new parameter effect. A call of Refresh_waveform() within
% Change_parameter() would cause the redundant computations of a new
% waveforms during each call before all intended parameter changes are done.
% To generate a new stimulus (e.g. fresh noise wavefrom), call Refresh_waveform()
% after Get_waveform().
% Set_specs(),Save_waveform(), plot_waveform(), and Play() and also any
% parameter change using the GUI (edit box, up & down button) do already call
% Refresh_waveform().
%
% Hint: After pressing "save, "plot" or "play" button, the waveform will
% be available in variable "ans" of the current workspace (command window).
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% Version history:
% 2.1  - added level meter
% 2.2  - alphabetically sorted functions instead  switch action (allows to
%        use the button menu 'show functions' of the Matlab editor)
% 2.3  - use gcbf and gcbo whenever there is no parameter 'h_fig'
% 2.4  - string of edit uicontrol determines the current parameter.
%      - changing paramter string of edit uicontrol does not update slider.
%      - slider callback changed
% 2.5  - levelmeter has now actions 'full_reset' and 'reset' (not red
%        overload LED)
% 2.6  - level meter using axis and rectangles instead of uicontrol frame
%        (the later is to memory extensive)
% 2.7  - Get_waveform(): reset level meter before generating new waveform
%        to visualize activity
%      - negative specs.param.stepsizes indicate multiplicative steps!
% 2.8  - new function Refresh_waveform(), called at the end of initialise(),
%        Set_specs(), Save_waveform(), plot_waveform(), and Play(). Should be
%        called externally after using Get_waveform() and Change_parameter(),
%		 best just after sound output initialisation to generate new waveform
%		 during sound output.
%      - function intended also for external calls begin with a capital
%        letter.
% 2.9  - initialise() uses Set_specs() to check and set parameter
%		 given in 'specs'.
%      - the used parameter are taken from slider 'Value'. This contains
%        always a valid value and prevents to Refresh_waveform() to use parameter
%        that is out of range from edit_box 'String'.
% 2.10- introduced fileparts.m to get pathname
%      - delete_stimulus_gui():rmpath is called only after LAST stimulus_gui
%        is being closed
% 3.0  - made compatible with stimuli_gui 3.0.
%      - new: Refresh_waveform_intervals() updates precomputed waveforms in
%        stimulus_gui and stimuli_gui
% 4.0  - data storage using setappdata/getappdata instead of in  'User'
%	   - h_parent_fig is stored (appdata). As long h_fig_parant is valid
%		 handle, delete_stimulus_gui() prevents deletion of stimulus_gui
%	   - Refresh_waveform() repaces Precompute_trials() and precompute().
%	     Appdata waveforms in stimuli_gui not longer updated but set to
%	     h_fig, or if appdata waveforms contains already handle of other
%	     fig, it is set [].
% 5.0  - sript 'general_code_A.m' will not be used anyumore. Its
%        functionality is now included in stimulus_gui(initialisation)
%        Feval switch yard is back in specific stimulus script.
%
% Part of the Stimulus Toolbox 
% Copyright (C) Torsten Marquardt
% Terms of the GNU General Public License apply
% (www.http://www.fsf.org/licensing/licenses/gpl.html).
vers = '5.0';


if ischar(varargin{1}), % INVOKE NAMED SUBFUNCTION OR CALLBACK
	if strcmp(varargin{1}, 'initialise'),
		h_fig = initialise(vers,varargin{2:end});
		if nargout > 0, varargout{1} = h_fig; end,
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

function buttondownFnc_text,
% in connection with'Select_edit_uicontrol':
%  - set associated edit uicontrol to CurrentObject
%  - deletes text box opened by fcn 'Select_edit_uicontrol'
% left click: bring descendant stimulus GUIs to foreground
% double click: reset level meter.
% right click: Ask to close all stimuli GUI

h_fig = gcbf;
[dummy no] = strtok(get(gcbo,'Tag'),'_');
set(h_fig,'CurrentObject', findobj(h_fig,'Tag',['edit' no])),
% delete dialog box if opened by fcn 'Select_edit_uicontrol'
delete(findobj(0,'Tag','stimulus_gui_select_edit_uicontrol_fig')),
switch get(h_fig,'SelectionType')
	case 'normal', % left mouse click: bring descendant stimulus GUIs to 
		% foreground (opened by 'duplicate')
		h_fig = getappdata(h_fig, 'h_child_duplication');
		while ishandle(h_fig),
			figure(h_fig),
			h_fig = getappdata(h_fig, 'h_child_duplication');
		end,
		figure(gcbf),
	case 'alt', % right mouse click: Close all stimuli GUI
		if strcmp(questdlg('Close all stimulus GUIs?', 'Stimuli GUI'),'Yes'),
			close(findobj(0,'Tag','stimulus_gui_fig')),
		end,
	case 'open', % double click: Reset level meter
		Level_meter(h_fig, 'full_reset');
end,
%% ==========================================================================  

function callback_slider,
% update edit window with new slider value

[dummy no] = strtok(get(gcbo,'Tag'),'_');
h_editUI = findobj(gcbf,'Tag',['edit',no]);
edit_value = str2num(get(h_editUI,'String'));
slider_value = get(gcbo,'Value');
slider_stepsize = get(gcbo,'User');
max_value = get(gcbo,'Max');
min_value = get(gcbo,'Min');

if slider_stepsize < 0,
	value = edit_value * abs(slider_stepsize).^sign(slider_value - edit_value);
else,
	value = edit_value + abs(slider_stepsize)*sign(slider_value - edit_value);
end,

if value > max_value,
	value = max_value;
end,

if value < min_value,
	value = min_value;
end,

Change_parameter(gcf, h_editUI, value);
%% ==========================================================================  

function err = Change_parameter(h_fig, parameter, new_value, verbose)
% !!! Use this routine to change paramter from external applications !!!
% check new string, check with slider range (min - max), when failure keep old value.
% 'parameter' can be handle or text string of parameter in GUI
%
% TO DO: include paramter 'frozen'

specs = getappdata(h_fig,'specs');

if ~exist('verbose','var')
    verbose = 1;
end,
err = 1; 
if ~exist('new_value','var') % manual parameter change in stimulus_gui GUI
    new_value=str2num(get(gcbo,'String'));
    h_fig = gcbf;
    parameter = gcbo;
end,

if ischar(parameter)
    for q = 1: length(specs.param),
        if strcmp(specs.param(q).name, parameter),
            parameter = findobj(h_fig, 'Tag', ['edit_' num2str(q)]);
        end,
    end,
    if ischar(parameter) && verbose % parameter string not found
        h = errordlg(['Parameter string ''' parameter ...
            ''' does not exist!'],'stimulus_gui.m - Change_parameter()');
        uiwait(h),
        return,
    end
end,

if ~ishandle(parameter)
    if verbose
        h = errordlg('Invalidstimuli_gui handle!', ...
            'stimulus_gui.m - Change_parameter()');
        uiwait(h),
    end
    return
end

[dummy no] = strtok(get(parameter,'Tag'),'_');
h_slider = findobj(h_fig,'Tag',['slider',no]);
max_value = get(h_slider,'Max');
min_value = get(h_slider,'Min');

if isempty(new_value)| length(new_value)>1,
	h = errordlg('Entry must be a single number!','stimulus_gui.m');
elseif new_value > max_value | new_value < min_value,
	h = errordlg([get(findobj(h_fig,'Tag',['text',no]),'String'),...
		' should be between ',num2str(min_value), ' and ', ...
		num2str(max_value),' and is ',num2str(new_value),'!'], ...
		'stimulus_gui.m - Change_parameter()');
end,

if exist('h', 'var'),
	new_value = get(h_slider,'Value');
	uiwait(h),
else
    err = 0;
end,

set(h_slider,'Value', new_value),
set(parameter,'String', num2str(new_value)),

Level_meter(h_fig, 'reset')

% delete precomputed wave data in stimulus_guiset dirty_flag
stimulus = getappdata(h_fig,'stimulus');
stimulus.waveform = 0;
setappdata(h_fig,'stimulus',stimulus),
% set dirty_flag flag
setappdata(h_fig,'dirty_flag',getappdata(h_fig,'h_parent_fig')),
%% ==========================================================================  

function hit = Check_and_reset_dirty_flag(h_fig,h_parent_fig),
% for parent applications, to check and reset appdata 'dirty_flag'

dirty_flag = getappdata(h_fig,'dirty_flag');

hit = 0;
if dirty_flag == h_parent_fig
    hit = 1;
end,
setappdata(h_fig,'dirty_flag',0)
%% ==========================================================================  

function close_childs(h_fig),
% close all descendant stimuli GUI opened with the 'duplicate' if chain is
% uninterupted

if ~exist('h_fig','var'), h_fig = gcbf; end;

if strcmp(get(h_fig,'SelectionType'), 'alt'),
	h_list = [];
	h = getappdata(h_fig, 'h_child_duplication');
	while ishandle(h),
		h_list = [h_list h];
		h = getappdata(h, 'h_child_duplication');
	end,
	close(h_list),
end,
%% ==========================================================================  

function delete_stimulus_gui(h_fig, h_master),
% DeleteFcn callback. Removes the folders ~\PLOTS from path
% Execute only if no parent figure exist

h_parent_fig = getappdata(h_fig, 'h_parent_fig');
if isempty(h_parent_fig) || ~ishandle(h_parent_fig) || ...
        exist('h_master','var') && h_parent_fig == h_master
	shh = get(0,'ShowHiddenHandles');
	set(0,'ShowHiddenHandles','on');
	if length(findobj(0,'Tag','stimulus_gui_fig')) == 1,
        pathname = fileparts(which('stimulus_gui'));
        if exist([pathname, '\PLOTS'],'dir')
            rmpath([pathname, '\PLOTS']),
        end
    end,
	set(0,'ShowHiddenHandles',shh);
	delete(h_fig),
end,
%% ==========================================================================  
	
function h_new = duplicate(h_fig),
% Make a copy of the stimulus GUI (same paramter). Handles of dupications
% are daisy chained.

if ~exist('h_fig','var'), h_fig = gcbf; end;

% get parameter of original figure
specs = Get_specs(h_fig);
specs.frozen = get(findobj(h_fig,'Tag','checkbox_frozen'),'Value');
% get position of last window
h_last = h_fig;
h_child_duplication = getappdata(h_fig, 'h_child_duplication');

%get position of last duplication
while ishandle(h_child_duplication),
	h_last = h_child_duplication;
	h_child_duplication = getappdata(h_child_duplication, ...
		'h_child_duplication');
end,
pos = get(h_last,'Position');

%open new window
h_new =stimulus_gui('initialise',[pos(1)+pos(3)+3 pos(2) pos(3)],specs,[]);
setappdata(h_last, 'h_child_duplication', h_new),
%% ==========================================================================  

function previous_state = Enable_uicontrols(h_fig, enable),
% disable or enable changing of parmeter

previous_state = get(findobj(h_fig,'Style','pushbutton'),'Enable');
previous_state = previous_state{1};

children = get(h_fig,'Children');
for q = 1:length(children),
	if strcmp(get(children(q),'Type'),'uicontrol'),
		if ~strcmp(get(children(q),'Style'),'text'),
			set(children(q),'Enable', enable),
		end,
	end,
end,
drawnow,
%% ==========================================================================  

function specs = Get_specs(h_fig),
% returns the parameter from the edit uicontrols

specs = getappdata(h_fig,'specs');
specs.frozen = get(findobj(h_fig,'Tag','checkbox_frozen'),'Value');
for q = 1: length(specs.param),
	specs.param(q).value = str2num(get(findobj(h_fig,'Tag', ...
		['edit_', num2str(q)]), 'String'));
end,
%% ==========================================================================  

function [waveform levels] = Get_waveform(h_fig),

stimulus = getappdata(h_fig,'stimulus'); % get last wave data
waveform = stimulus.waveform;
levels = stimulus.levels;

if size(waveform,1) < 2, % i.e. waveform is zero or empty
	if ~(isempty(waveform)| waveform == -1), 
        % stimulus generation (Refresh_waveform) not in progress or failed
		Refresh_waveform(h_fig),
		[waveform levels] = Get_waveform(h_fig);
    end,
    stimulus.waveform = waveform;
    stimulus.levels = levels;
    setappdata(h_fig,'stimulus', stimulus),
end,
%% ==========================================================================  

function h_fig = initialise(varargin)
% varargin = [vers, position, specs, h_parent_fig],
global SAMPLE_RATE,

vers = varargin{1};
stimulus_filename = varargin{2};

%add \General to path
pathname = fileparts(which(stimulus_filename));%get path
if exist([pathname,'\PLOTS'],'dir')
    addpath([pathname,'\PLOTS'])
end
    
if nargin<3||isempty(varargin{3}), %called without arguments
    pos = [0 0 150];
else
    pos=varargin{3}; %position given in 1st argument
    if length(pos) ~= 3, % position does not contain y-dimension of figure
        h = errordlg('Position must be a 3 element vector! No Y-dimension!', ...
            'Stimulus GUI ini');
        uiwait(h),
        return,
    end,
end
if pos(3)<150, pos(3) = 150, end,

if nargin>3 && isstruct(varargin{4}) % specs given
    specs = varargin{4};
else
    eval(['specs = ',stimulus_filename, ...
        '(''return_default_parameter'',vers);']);
end


if nargin>4  % called by parent application
    h_parent_fig = varargin{5};
elseif nargin>3 && ishandle(varargin{4}),
    h_parent_fig = varargin{4};
else
    h_parent_fig = []; % random generator initialisation only if no parant application
    rand('twister',sum(100*clock)), % DO NOT MODIFY!
end

h_fig = figure('Color',get(0,'defaultUicontrolBackgroundColor'), ...
	'Backingstore','off', ...
	'HandleVisibility','on', ...
	'MenuBar','none', ...
	'Name', specs.name, ...
	'NumberTitle','off', ...
	'Position',[pos length(specs.param)*24 + 15], ...
	'Resize','off', ...
	'Tag','stimulus_gui_fig');

set(h_fig,'CloseRequestFcn', ...
	['stimulus_gui(''delete_stimulus_gui'',' num2str(h_fig) ')']),

% initialis ALL appdata (list complete)
stimulus.waveform = 0;
stimulus.levels = [];
setappdata(h_fig,'dirty_flag', 1);
setappdata(h_fig,'stimulus',stimulus);
setappdata(h_fig,'specs',specs), %specs contains the stimulus parameter
setappdata(h_fig,'h_parent_fig',h_parent_fig),
setappdata(h_fig, 'h_child_duplication',[]),
if exist('lsp_transfer_fcn.mat','file')
    load('lsp_transfer_fcn.mat')
    if SAMPLE_RATE/2 > length(lsp_transfer_fcn)
        h = errordlg(...
            'lsp_transfer_fcn.mat does not cover stimulus frequency range',...
            'stimulus_gui.m - initialise()');
        uiwait(h),
        close(h_fig)
        return,
    end
    setappdata(h_fig, 'lsp_transfer_fcn',lsp_transfer_fcn),
end

uicontrol('Parent',h_fig, ...
	'Position',[round(pos(3)/200) 1 round(pos(3)/3) 18], ...
	'String','frz', ...
	'Style','checkbox', ...
	'Value', specs.frozen, ...
	'Tag','checkbox_frozen', ...
	'TooltipString', 'Reuse last wave form if no paramter change.' );

uicontrol('Parent',h_fig, ...
	'Callback','ans = stimulus_gui(''Save_waveform'',gcbf);', ...
	'Position',[round(pos(3)/4) 1  round(pos(3)/5) 18], ...
	'String','save', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_save', ...
	'UIContextMenu',[]);

uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback','ans = stimulus_gui(''plot_waveform'');', ...
	'Position',[round(pos(3)/2.2) 1  round(pos(3)/6.1) 18], ...
	'String','plot', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_plot', ...
	'UIContextMenu',[], ...
	'TooltipString', [specs.name ', GUI vers. ' vers]);

uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'buttondownFcn','stimulus_gui close_childs', ...
	'Callback','stimulus_gui duplicate', ...
	'Position',[round(pos(3)/1.62) 1  round(pos(3)/5) 18], ...
	'String','dupl', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_duplicate', ...
	'TooltipString', ...
	'Make a copy of the stimulus GUI (same paramter). Right click deletes them.');

uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback','ans = stimulus_gui(''Play'');', ...
	'Position',[round(pos(3)/1.22) 1 round(pos(3)/5) 18], ...
	'String','play', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_play', ...
	'UIContextMenu',[]);

for q = 1:length(specs.param),
	uicontrol('Parent',h_fig, ...
		'buttondownFcn', 'stimulus_gui buttondownFnc_text', ...
		'Enable','inactive', ...
		'HorizontalAlignment','right', ...
		'Position', ...
		[0 (q-1)*24+20 round(pos(3)-17-round(pos(3)/3)) 17], ...
		'String',[specs.param(q).name,':'], ...
		'Style','text', ...
		'Tag',['text_', num2str(q)]);

	h = uicontrol('Parent',h_fig, ...
		'BackgroundColor',[1 1 1], ...
		'Callback', 'stimulus_gui(''Change_parameter'',gcbf);', ...
		'Position',[round(pos(3)-17-round(pos(3)/3)) ...
		(q-1)*24+20 round(pos(3)/3) 20], ...
		'Style','edit', ...
		'Tag',['edit_', num2str(q)], ...
		'TooltipString', specs.param(q).ToolTip);

    % mark the existence of lsp_transfer_fcn
    if strfind(specs.param(q).name,'[dB')
        set(h,'BackgroundColor',[1 0.7 0.7])
    end
    
	uicontrol('Parent',h_fig, ...
		'Callback', 'stimulus_gui(''callback_slider''),', ...
		'Max',specs.param(q).max, ...
		'Min',specs.param(q).min, ...
		'Value',specs.param(q).min, ...
		'Position',[round(pos(3)-13) (q-1)*24+19 11 22], ...
		'Style','slider', ...
		'User',specs.param(q).stepsize, ...
		'Tag',['slider_', num2str(q)]);
end,

Set_specs(h_fig,specs); % checks specs.param(q).value before setting

% LED levelmeter
figure(h_fig),
axes('Units','pixels', ...
	'Position', [round(pos(3)-16), 21, 16, length(specs.param)*24], ...
	'XTick',[],'YTick',[], ...
	'XColor', get(0,'defaultUicontrolBackgroundColor'), ...
	'YColor', get(0,'defaultUicontrolBackgroundColor'), ...
	'XLim',[-1 1],'YLim',[0 length(specs.param)*24],...
	'Color', get(0,'defaultUicontrolBackgroundColor'), ...
	'Tag','Level_meter');

for q = 0:4*length(specs.param)-2,
	led = rectangle('Position', [-1 1+6*q .4 5], ...
		'EdgeColor', [.4 .4 .4], ...
		'FaceColor', [0 0 0]);
	set(led,'Tag',['left_LED_', num2str(q+1)])

	led = rectangle('Position', [.7 1+6*q .4 5], ...
		'EdgeColor', [.4 .4 .4], ...
		'FaceColor', [0 0 0]);
	set(led,'Tag',['right_LED_', num2str(q+1)])
end,

Level_meter(h_fig, 'full_reset');
%% ==========================================================================  

function Level_meter(h_fig, action, levels),

switch action,
	case 'full_reset', % resets includes red overload LEDs
		set(findobj(h_fig,'Type','rectangle'), 'FaceColor',[0 0 0]),

	case 'reset', % keeps the state of red overload LEDs
		overload = findobj(h_fig,'FaceColor', [1 0 0]);
		set(findobj(h_fig,'Type','rectangle'), 'FaceColor', [0 0 0]),
		set(overload, 'FaceColor', [1 0 0]),

	case 'update',
		% updates the level meter. never resets the uppermost red overload LEDs.
		n_LEDs = 4*length(getfield(getappdata(h_fig,'specs'),'param')) - 1;
		if ~isempty(levels),
			set_level_meter(h_fig, n_LEDs, levels, 1),
		else,
			Level_meter(h_fig,'reset'),
		end,
	otherwise, error('Fcn ''Level_meter'': Unknown action')
end,
drawnow,
%% ==========================================================================  

function waveform = Play(h_fig),
% left click: play wave form
% Not finished! right click: open dialogbox for play settings (fs, nBits, out_ch)

global SAMPLE_RATE,

if ~exist('h_fig','var'), h_fig = gcbf; end;

Level_meter(h_fig, 'reset'), % leave in, confirmes play button press!
[waveform levels] = Get_waveform(h_fig);
if size(waveform,1) > 1,
	try,
		if strcmp(computer, 'PCWIN'),
			wavplay(waveform,SAMPLE_RATE,'async');
		else,
			sound(waveform,SAMPLE_RATE,16);
		end,
	catch
		warning(lasterr);
	end,
end,
Level_meter(h_fig, 'update', levels),
Refresh_waveform(h_fig),
%% ==========================================================================  

function waveform = plot_waveform
% Not finished!
% left click: plot wave form in new figure
% right click: open context menu with plot options

global SAMPLE_RATE,

Level_meter(gcbf, 'reset'), % leave in, confirmes plot button press!
[waveform levels] = Get_waveform(gcbf);
if size(waveform,1) > 1,
	[m n]=size(waveform); m=2*floor(m/2); if m < 3, return,  end,

	% if ~exist('PLOTS','dir'),
	%     % create context menu of what plots are available in directory 'PLOTS'
	%     return,
	% end,

	h = figure;
	set(h,'Name','wave'),
	% plot time course
	figure(h),
	[m n]=size(waveform);
	x=1000*linspace(0,m/SAMPLE_RATE,m)';
	plot(x,waveform(:,1)), hold on,
	plot(x,waveform(:,n),'r'), grid on, hold off, title(''), zoom on,
end,
Level_meter(gcbf, 'update', levels),
Refresh_waveform(gcbf),

return

% % plot spectrum
% x=linspace(0,SAMPLE_RATE/2,m/2);
% H_waveform=fft(waveform)/sqrt(length(waveform)/SAMPLE_RATE)/200; % Use this for noise
% plot([x,max(x)+x]', 20*log10(abs([H_waveform(1:m/2,1); H_waveform(m/2:-1:1,n)])),'.'),
% % H_waveform=fft(waveform)/length(waveform); % Use this for tones
% %plot([x,max(x)+x]', 20*log10(abs([H_waveform(1:m/2,1); H_waveform(m/2:-1:1,n)])) + ...
% %    fix_specs.param(1).max+6.02,'.'),
% grid, title('Press any key to view IPD!'),zoom on,
% pause,
%
% % plot phase with same X limits
% figure(h),
% x_lim=xlim;
% plot(x,angle(H_waveform(1:m/2,1))-angle(H_waveform(1:m/2,n)),'.'),
% grid, title('Press any key to view time course!'), zoom on, xlim(x_lim),pause,
%% ==========================================================================  

function Refresh_waveform(h_fig),
% generate new waveform 

stimulus = getappdata(h_fig,'stimulus');

previous_state = Enable_uicontrols(h_fig, 'off');

if ~(get(findobj(h_fig,'Tag','checkbox_frozen'),'Value') & ...
        stimulus.waveform ~=0), % not if frozen and up-to-date waveform available

    % set dirty_flag flag for tests by parent applications
    setappdata(h_fig,'dirty_flag',getappdata(h_fig,'h_parent_fig')),
    
    %set stimmulus to zero while computing new waveform
    stimulus.waveform = 0;
    stimulus.levels = [];
    setappdata(h_fig,'stimulus',stimulus); % store the wave data
    
    % calculate new waveform
    specs = Get_specs(h_fig);
    stimulus.waveform = eval([strtok(specs.name), ...
        '(''generate_waveform'', specs, h_fig)']);
	if length(stimulus.waveform)< 2,
		stimulus.waveform = -1;
		warning('Stimulus_gui-Refresh_waveform: waveform refresh failed!')
	end,
    
    % compute levels for level meter
    if isempty(stimulus.waveform),
        stimulus.levels = [];
    else,
        max_levels = 20*log10(max(abs(stimulus.waveform)));
        for q = 1:length(max_levels),
            if max_levels(q) == 0, % log(0) gives warning!
                rms_levels(q) = max_levels(q);
            else,
                rms_levels(q) = ...
                    20*log10(sqrt(mean(stimulus.waveform(:,q).^2)));
            end,
        end,
        stimulus.levels = [max_levels; rms_levels];
    end,

    setappdata(h_fig,'stimulus',stimulus); % store the wave data
end,

Enable_uicontrols(h_fig, previous_state);

%% ==========================================================================  

function waveform = Save_waveform(h_fig, filename),
%%% Not finished!
% left click: open file dialog box, save wave form as WAV file
%%% right click: open context menu with saving options

global SAMPLE_RATE,

[waveform levels] = Get_waveform(h_fig);
if size(waveform,1) > 1,
	if ~exist('filename','var'),
		[filename,pathname] = uiputfile(['c:\Windows\Temp\temp.wav'],'Save wave as:');
		[filename suffix]=strtok(filename,'.');
		filename =[pathname,filename];
	end,

	if strcmp(suffix,'.wav'),
		wavwrite(waveform,SAMPLE_RATE,16,filename),
	elseif strcmp(suffix,'.mat'),
		save( filename,'waveform'),
	end,
end,
Level_meter(h_fig, 'update', levels),
Refresh_waveform(h_fig),
%% ==========================================================================  

function parameter = Select_edit_uicontrol(h_fig, parameter),

if exist('parameter'), % variable has only valid umbers. Return handles.
	parameter.handle = findobj(h_fig,'Tag',['edit_' num2str(parameter.number)]);
	return,
end,

set(h_fig,'CurrentObject',0),
h_msg = warndlg('Select parameter by mouse click on text label!','Stimuli GUI');
delete(findobj(h_msg,'Style','pushbutton')),
set(h_msg,'Tag','stimulus_gui_select_edit_uicontrol_fig'),
uiwait(h_msg),
parameter.handle = get(h_fig,'CurrentObject');
tag = get(parameter.handle,'Tag');
parameter.number = str2num(tag(length(tag)));
%% ==========================================================================

function set_level_meter(h_fig, n_LEDs, levels, stepsize_dB),

max_levels = levels(1,:);
rms_levels = levels(2,:);

if length(max_levels) == 1, % mono stimulus
    max_levels = [max_levels max_levels];
    rms_levels = [rms_levels rms_levels];
end,

leds = zeros(n_LEDs,size([max_levels max_levels; rms_levels rms_levels],2));
try,
	leds = hist([max_levels max_levels; rms_levels rms_levels], ...
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
	set(findobj(h_fig,'Tag',['left_LED_' num2str(q)]), ...
		'FaceColor', color),

	if leds(q,2),
		color = [0 1 0]; % green
		if q > n_LEDs-7, color = [1 1 0]; end, % yellow
	else,
		color = [0 0 0];
	end,

	set(findobj(h_fig,'Tag',['right_LED_' num2str(q)]), ...
		'FaceColor', color),
end,

% set uppermost red overload LEDs if max_level > 1. Never overrides an overload.
if leds(n_LEDs,1),
	set(findobj(h_fig,'Tag',['left_LED_' num2str(n_LEDs)]), ...
		'FaceColor',[1 0 0]),
end,
if leds(n_LEDs,2),
	set(findobj(h_fig,'Tag',['right_LED_' num2str(n_LEDs)]), ...
		'FaceColor',[1 0 0]),
end,
%% ==========================================================================

function error = Set_specs(h_fig,specs),
% sets the edit and slider uicontrols as specified in specs

error = 0;
for q = 1:length(specs.param),
	h_edit = findobj(h_fig,'Tag',['edit_' num2str(q)]);
	set(h_fig, 'CurrentObject', h_edit);
	if Change_parameter(h_fig, h_edit, specs.param(q).value),
		error = error + 1;
	end,
end,
set(findobj(h_fig,'Tag','checkbox_frozen'),'Value',specs.frozen)
Refresh_waveform(h_fig),