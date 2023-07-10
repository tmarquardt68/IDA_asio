function varargout = ida_exp_gui(varargin)
% TO DO:
% - %%%% !!! put these into exp figure text boxes !!!!
%         disp(['Rpt: ' num2str(q), ' of ' num2str(stimulus.n_rpts)...
%             ', intvl: ' num2str(q2), '(' num2str(interval_no) ')' ...
%             ' of ' num2str(size(parameter_table,1))])


if (nargout)
    [varargout{1:nargout}] = feval(varargin{1},varargin{2:end});
else
    feval(varargin{1}, varargin{2:end}); % FEVAL switchyard
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h_fig = initialise
h_fig = uifigure;
h_fig.HandleVisibility = 'on';
h_fig.Tag = 'ida_exp_fig';
h_fig.Position = [0 60 120 90];
h_fig.Resize = 0; 
h_fig.CloseRequestFcn  = 'ida_exp_gui CloseRequestFcn';

uibutton(h_fig, ...
	'text','pause', ...
	'ButtonPushedFcn',@CB_pause_button, ...
	'Tag','PB_pause', ...
	'Position',[5 5 55 20]);

uibutton(h_fig, ...
	'text','stop', ...
	'ButtonPushedFcn',@CB_stop_button, ...
	'Tag','PB_stop', ...
	'Position',[65 5 55 20]);

uieditfield(h_fig, ...
	'Value','counter', ...
    'FontSize',9, ...
	'Editable','off', ...
	'Tag','txt_counter', ...
	'Position',[0 30 120 15]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CB_pause_button(callingObj,~)
if strcmp(get(callingObj,'text'),'Pause')
    set(callingObj,'text','Continue')
else
    set(callingObj,'text','Pause')
end

%% =========================================================================
function CB_stop_button(callingObj,~)
set(callingObj,'text','Stopped')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CloseRequestFcn
h_exp_fig =findobj('Tag','ida_exp_fig');
h_stimuli_gui = getappdata(h_exp_fig(end), 'h_stimuli_gui');

delete(h_exp_fig(end)),

% % last sound still needs to finish!
% presentation_period = getappdata(h_exp_fig,'presentation_period');
% pause(1.1*presentation_period)

PsychPortAudio('Close')
if ishandle(h_stimuli_gui)
    close(h_stimuli_gui)
end
cd('../..') % ?

