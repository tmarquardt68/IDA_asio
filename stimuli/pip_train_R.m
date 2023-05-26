function output = stimulus(varargin),
%[duration,cos_ramp,delay,freq,ild,ipd,itd]
% Generates tone bursts for Zwicker's masking period pattern stimulus on left channel .
%
% GENERAL: This functions specifies (case 'return_default_parameter') and
% generates (case 'generate_waveform') a stimulus according to the m-file name.
% The function can be used in conjunction with stimulus_gui.m (>> help stimulus_gui).
%
% Version history:
%
% Part of the Stimulus Toolbox 
% Copyright (C) Torsten Marquardt
% Terms of the GNU General Public License apply
% (www.http://www.fsf.org/licensing/licenses/gpl.html).
% Date: 14-10-13
vers = '0.1';

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
defaults.param(n).max = MAX_LEVEL(2);              % is set correctly (sound_card
defaults.param(n).stepsize = 3;                 % specific - needs calibration).
defaults.param(n).value = 82;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'freq [Hz]';
defaults.param(n).min = 125;
defaults.param(n).max = SAMPLE_RATE/2;
defaults.param(n).stepsize = -2^(1/3);
defaults.param(n).value = 1000;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'rate [Hz]';
defaults.param(n).min = 2;
defaults.param(n).max = 70;
defaults.param(n).stepsize = .5;
defaults.param(n).value = 16;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'phase_pos [deg]';
defaults.param(n).min = -360;
defaults.param(n).max = 360;
defaults.param(n).stepsize = 22.5;
defaults.param(n).value = 0;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'mod_idx [%]';
defaults.param(n).min = 0;
defaults.param(n).max = 100;
defaults.param(n).stepsize = 5;
defaults.param(n).value = 0;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'mod_rate [Hz]';
defaults.param(n).min = .25;
defaults.param(n).max = 8;
defaults.param(n).stepsize = 5;
defaults.param(n).value = 2;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'duration [s]';
defaults.param(n).min = 0;
defaults.param(n).max = 60;
defaults.param(n).stepsize = -2; % negative values indicate multiplicative steps!
defaults.param(n).value = 3;
defaults.param(n).ToolTip = 'No ToolTip';

defaults.frozen = 0;    % If non-zero: Reuse last wave form if no paramter change.

%==========================================================================
function wave = generate_waveform(specs, h_fig),

%% Begin: DO NOT MODIFY ! Usually the same for all stimuli.
global SAMPLE_RATE MAX_LEVEL %Set by caller, otherwise by 'return_default_parameter'

duration = [];
wave = 0;
for q = 1: length(specs.param),

	eval([strtok(specs.param(q).name) '=' ...
		num2str(specs.param(q).value) ';']);
end,


%% End: DO NOT MODIFY ! Usually the same for all stimuli.
%//////////////////////////////////////////////////////////////////////////
% HERE WE GO: Specify here the code to generate the sound with the values in
% 'specs', a struct like returned by 'return_default_values' specified above
%//////////////////////////////////////////////////////////////////////////
% First check and if necessary transform the parameter
%//////////////////////////////////////////////////////////////////////////

amplitude(1) = 0;
amplitude(2) = 10^((level - MAX_LEVEL(2))/20);

% % compensate amplitude for transfer function
% if exist('h_fig','var')
%     lsp_transfer_fcn = getappdata(h_fig, 'lsp_transfer_fcn');
%     if ~isempty(lsp_transfer_fcn)
%         amplitude = amplitude./lsp_transfer_fcn(probe_freq,:);
%     end
%     if amplitude > 1
%         amplitude = 1;
%         warning('Compensation by lsp_transfer_fcn.mat: LF signal amplitude >1!')
%         return
%     end
% end

if (freq > SAMPLE_RATE/2) | (freq <= 0),
	h = errordlg('Wrong parameter setting: freq > fs/2!',...
		mfilename); uiwait(h),  return,
end,

duration = round(duration * SAMPLE_RATE);
period = round(SAMPLE_RATE/rate);
% generate tone pip %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ 
pip_length = 3/freq * SAMPLE_RATE;
tone_pip = sin(linspace(0,freq/SAMPLE_RATE*pip_length*2*pi*(1-1/pip_length), ...
	pip_length));
tone_pip = tone_pip .* cos(pi/2:pi/(pip_length-1):1.5*pi).^2;

wave = zeros(duration,2);
phase_pos = round(period*mod(phase_pos,360)/360);

ptr = 1;
while ptr + length(tone_pip) < duration 
    wave(ptr + phase_pos : ptr + phase_pos+length(tone_pip)-1,:) = [tone_pip; tone_pip]';
    ptr = ptr + period;
end

% modulation envelope
mod_env = mod_idx/200*sin(linspace(0,mod_rate/SAMPLE_RATE*duration*2*pi*(1-1/duration), ...
	duration) + pi*rand) + 1 - mod_idx/200;
wave(:,1) = wave(:,1) .* mod_env';
wave = wave .* (amplitude' * mod_env)';
%wave(:,2) = .1* mod_env';

% silence_len = floor(SAMPLE_RATE/rate-length(tone_pip));
% wave = zeros(1,floor(SAMPLE_RATE/rate));
% phase_pos = round(length(wave)*mod(phase_pos,360)/360)+1;
% wave(phase_pos : phase_pos+length(tone_pip)-1) = tone_pip;
% 
% % % Gaussian window
% % period = floor(duration/rate);
% % x = linspace(0,freq/SAMPLE_RATE*period*2*pi*(1-1/period),period);
% % tone_pip = amplitude * sin(x);
% % tone_pip = tone_pip.*gausswin(length(x),22000*freq/SAMPLE_RATE)';
% 
% wave = [zeros(1,2*length(wave)) ...
%     repmat(wave,1,floor(duration*rate/SAMPLE_RATE)-3) ...
%     zeros(1,length(wave))]; 
% wave = [wave; zeros(1,length(wave))]';
% Power: p = sum(wave.^2)/length(wave)
%% avg. spectrum level
% 10*log10(max(abs(fft(wave).^2)))