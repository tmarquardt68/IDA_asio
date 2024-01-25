function varargout = exp_gui(varargin)

if (nargout),
    [varargout{1:nargout}] = feval(varargin{1},varargin{2:end});
else
    feval(varargin{1}, varargin{2:end}); % FEVAL switchyard
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CB_pause_button

if strcmp(get(gcbo,'String'),'Pause')
    set(gcbo,'String','Continue')
else
    set(gcbo,'String','Pause')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CB_stop_button

set(gcbo,'String','Stopped')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CloseRequestFcn

h_exp_fig =findobj('Tag','ida_exp_fig');
h_sound_device = getappdata(h_exp_fig,'h_sound_device');
h_stimuli_gui = getappdata(h_exp_fig, 'h_stimuli_gui');
presentation_period = getappdata(h_exp_fig,'presentation_period');

delete(h_exp_fig),

% last sound still needs to finish!
pause(1.1*presentation_period)
PsychPortAudio('Close')
if ishandle(h_stimuli_gui)
    close(h_stimuli_gui)
end,
