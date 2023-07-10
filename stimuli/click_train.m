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
defaults.param(n).name = 'duration [ms]';
defaults.param(n).min = 400;
defaults.param(n).max = 1000;
defaults.param(n).stepsize = 50;
defaults.param(n).value = 500;
defaults.param(n).ToolTip = 'length of the period to record';

n = n+1;
defaults.param(n).name = 'num_clicks';
defaults.param(n).min = 1;
defaults.param(n).max = 50;
defaults.param(n).stepsize = 1;
defaults.param(n).value = 10;
defaults.param(n).ToolTip = 'Number of clicks in the train';

n = n+1;
defaults.param(n).name = 'click_rate [Hz]';
defaults.param(n).min = 1;
defaults.param(n).max = 500;
defaults.param(n).stepsize = 1;
defaults.param(n).value = 20;
defaults.param(n).ToolTip = 'Rate of click train';

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

wave=zeros(1,round(duration./1000*SAMPLE_RATE));

clickSamples=round(50e-6*SAMPLE_RATE);
iciSamples=round(SAMPLE_RATE/click_rate);
for n=1:num_clicks
    wave((1:clickSamples)+(n-1)*iciSamples)=.9;
end

wave(2,:)=wave(1,:);
wave=wave';