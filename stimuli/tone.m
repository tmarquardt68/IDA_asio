function output = stimulus(varargin),
%[duration,cos_ramp,delay,freq,ild,ipd,itd]
% General purpose dichotic tone.
%
% GENERAL: This functions specifies (case 'return_default_parameter') and
% generates (case 'generate_waveform') a stimulus according to the m-file name.
% The function can be used in conjunction with stimulus_gui.m (>> help stimulus_gui).
%
% Version history:
% 2.0  - made compatible with Stiumulus_GUI Vers. 2
% 2.1  - change case to function calls
% 2.2  - include  parameter "Phase"
% 2.3  - - ILD keeps the binaural average amplitude equal (by adding the
% identical waveform to one side, and subtracting it from the other.
%
% Part of the Stimulus Toolbox 
% Copyright (C) Torsten Marquardt
% Terms of the GNU General Public License apply
% (www.http://www.fsf.org/licensing/licenses/gpl.html).
% Date: 17-07-10
vers = '2.3';

%% Begin: DO NOT MODIFY!
if nargin > 0 && strcmp(varargin{1},'generate_waveform')
	output = generate_waveform(varargin{2:end});
elseif nargin > 0 && strcmp(varargin{1},'return_default_parameter')
		output = return_default_parameter(vers);
else
    h_fig = stimulus_gui('initialise',mfilename, varargin{:});
    if nargout > 0, output = h_fig; end
end
%% End: DO NOT MODIFY!
%==========================================================================

function defaults = return_default_parameter(vers), 
% specify variables and default values of sound.
% Stepsize is for the slider in the sound figure

defaults.name = [mfilename ' ' vers]; % DO NOT MODIFY!
n = 0; % DO NOT MODIFY!

global SAMPLE_RATE   MAX_LEVEL, % if not set by caller set default values
if isempty(SAMPLE_RATE), SAMPLE_RATE = 48000; end,
if isempty(MAX_LEVEL), MAX_LEVEL = 110; end, % dependant on sound card
% A tone at MAX_LEVEL has amplitude of one.

% !!! KEEP spectrum level always as the 1st parameter (for plot function)!
n = n+1;
defaults.param(n).name = 'level [dB]';     % spectrum level in
defaults.param(n).min = -40;                    % dB SPL if MAX_LEVEL
defaults.param(n).max = MAX_LEVEL;              % is set correctly (sound_card
defaults.param(n).stepsize = 3;                 % specific - needs calibration).
defaults.param(n).value = 75;
defaults.param(n).ToolTip = ...
	'Level in dB SPL if MAX_LEVEL (global) is set correctly (Soundcard specific).';

n = n+1;
defaults.param(n).name = 'ITD [us]';
defaults.param(n).min = -8000;
defaults.param(n).max = 8000;
defaults.param(n).stepsize = 100;
defaults.param(n).value = 0;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'IPD [deg]';
defaults.param(n).min = -360;
defaults.param(n).max = 360;
defaults.param(n).stepsize = 5;
defaults.param(n).value = 0;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'ILD [dB]';
defaults.param(n).min = -100;
defaults.param(n).max = 100;
defaults.param(n).stepsize = 1;
defaults.param(n).value = 0;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'freq [Hz]';
defaults.param(n).min = 0;
defaults.param(n).max = SAMPLE_RATE/2;
defaults.param(n).stepsize = 100;
defaults.param(n).value = 500;
defaults.param(n).ToolTip = 'No ToolTip';n = n+1;

defaults.param(n).name = 'Phase [deg]';
defaults.param(n).min = 0;
defaults.param(n).max = 360;
defaults.param(n).stepsize = 45;
defaults.param(n).value = 0;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'delay [ms]';
defaults.param(n).min = 0;
defaults.param(n).max = 20000;
defaults.param(n).stepsize = 5;
defaults.param(n).value = 0;
defaults.param(n).ToolTip = 'Silence (prepend zeros) before stimulus starts.';

n = n+1;
defaults.param(n).name = 'cos_ramp [ms]';
defaults.param(n).min = 0;
defaults.param(n).max = 500;
defaults.param(n).stepsize = 1;
defaults.param(n).value = 10;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'duration [ms]';
defaults.param(n).min = 0;
defaults.param(n).max = 20000;
defaults.param(n).stepsize = -2; % negative values indicate multiplicative steps!
defaults.param(n).value = 1000;
defaults.param(n).ToolTip = 'No ToolTip';

defaults.frozen = 1;    % If non-zero: Reuse last wave form if no paramter change.

%==========================================================================
function wave = generate_waveform(specs, h_fig),

%% Begin: DO NOT MODIFY ! Usually the same for all stimuli.
global SAMPLE_RATE   MAX_LEVEL, %Set by caller, otherwise by 'return_default_parameter'

wave = 0;
for q = 1: length(specs.param),
	if exist(strtok(specs.param(q).name),'var'),
		h = errordlg('Parameter name dublication!',...
			mfilename); uiwait(h), return,
	end,
	eval([strtok(specs.param(q).name) '=' ...
		num2str(specs.param(q).value) ';']);
end,

if exist('level'),
	amplitude(2) = 10^((level - specs.param(1).max)/20); % That's why "spct_level is
	d = tanh(ILD./(20*log10(60)/.5/log(60)));       % always 1st control !
    amplitude(1) = amplitude(2)*(1-d);    
	amplitude(2) = amplitude(2)*(1+d);
	%%%%% END: DO NOT MODIFY !
else,
	% Don't call the level "level" and define your own amplitude here.
end,
%% End: DO NOT MODIFY ! Usually the same for all stimuli.
%//////////////////////////////////////////////////////////////////////////
% HERE WE GO: Specify here the code to generate the sound with the values in
% 'specs', a struct like returned by 'return_default_values' specified above
%//////////////////////////////////////////////////////////////////////////
% First check and if necessary transform the parameter
%//////////////////////////////////////////////////////////////////////////

% compensate amplitude for transfer function
if exist('h_fig','var')
    lsp_transfer_fcn = getappdata(h_fig, 'lsp_transfer_fcn');
    if ~isempty(lsp_transfer_fcn)
        amplitude = amplitude./abs(lsp_transfer_fcn(freq,:));
        Phase = Phase - angle(lsp_transfer_fcn(freq,:));
    end
    if amplitude > 1
        warning('Compensation by lsp_transfer_fcn.mat: signal amplitude >1!')
        return
    end
end

ITD = ITD/1000000;
IPD = IPD/180*pi;
Phase = Phase/180*pi;

if (freq > SAMPLE_RATE/2) | (freq <= 0),
	h = errordlg('Wrong parameter setting: freq > fs/2!',...
		mfilename); uiwait(h),  return,
end,

delay = round(delay/1000 * SAMPLE_RATE);
cos_ramp = round(cos_ramp/1000 * SAMPLE_RATE);
if cos_ramp == 1, cos_ramp = 0; end,
duration = round(duration/1000 * SAMPLE_RATE);
if (2*cos_ramp + delay) > duration,
	h = errordlg(...
		'Wrong parameter settings: (2*cos_ramp + delay) > duration!', ...
		mfilename); uiwait(h),  return,
end,
clear wave

% generate wave form %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
wave(2,:) = cos(linspace(0,freq/SAMPLE_RATE*duration*2*pi*(1-1/duration), ...
	duration) + Phase(2) + IPD/2+ITD*freq*pi);

wave(1,:) = cos(linspace(0,freq/SAMPLE_RATE*duration*2*pi*(1-1/duration), ...
	duration) + Phase(1) -(IPD/2+ITD*freq*pi));

% apply envelope
wave = wave' .* (amplitude' *...
	[zeros(1,delay),cos(pi/2:(pi-pi/2)/(cos_ramp-1):pi).^2,...
	ones(1,duration-delay-2*cos_ramp),cos(0:pi/2/(cos_ramp-1):pi/2).^2])';

% Power: p = sum(wave.^2)/length(wave)