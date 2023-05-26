function output = stimulus(varargin),
%
% pink_noise.m produces binaural bandpass noise with optionally pink
% spectrum.
%[duration,cos_ramp,bandwitdh,freq,ild,ipd,itd,coherence, pink_spctrm]
%
%
% GENERAL: This m-file specifies default paramerter
% (fcn 'return_default_parameter()') and generates a stimulus 
% (fcn 'generate_waveform()') according to the name of thhe m-file.
% The function can be used in conjunction with stimulus_gui.m 
% (>> help stimulus_gui).
%
% Version history:
% 1.0  - derived from bp_noise 3.0
% 1.1  - ILD keeps the binaural average amplitude equal (by adding the
% identical waveform to one side, and subtracting it from the other.
%
% Part of Stimulus Toolbox by Torsten Marquardt
% Date: 17-07-10
vers = '1.1';

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
defaults.param(n).name = 'spct_level [dB]';     % if pink = 1, 
defaults.param(n).min = -40;                    % the reference is arbitrary
defaults.param(n).max = MAX_LEVEL;  
defaults.param(n).stepsize = 5;                 
defaults.param(n).value = 50;
defaults.param(n).ToolTip = ...
	'If pink = 1, the level in dB has arbitrary refenrence';

n = n+1;
defaults.param(n).name = 'pink_spctrm';
defaults.param(n).min = 0;
defaults.param(n).max = 1;
defaults.param(n).stepsize = 1;
defaults.param(n).value = 1;
defaults.param(n).ToolTip = 'pink = 1: spectrum slope of -3dB/octave ';

n = n+1;
defaults.param(n).name = 'coherence';   % Amount of dispersion of IPD & ILD between
defaults.param(n).min = 0;              % the spectral components. Coherence is:
defaults.param(n).max = 1;              %               sum(R.*conj(L))
defaults.param(n).stepsize = .1;        % ----------------------------------------
defaults.param(n).value = 1;            % sqrt(sum(R.*conj(R)) .* sum(L.*conj(L)))
defaults.param(n).ToolTip = ...
	'Amount of dispersion of IPD between the spectral components.';

n = n+1;
defaults.param(n).name = 'ITD [us]';
defaults.param(n).min = -9999;
defaults.param(n).max = 9999;
defaults.param(n).stepsize = 100;
defaults.param(n).value = 0;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'IPD [deg]';
defaults.param(n).min = -360;
defaults.param(n).max = 360;
defaults.param(n).stepsize = 45;
defaults.param(n).value = 0;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'ILD [dB]';
defaults.param(n).min = -99;
defaults.param(n).max = 99;
defaults.param(n).stepsize = 3;
defaults.param(n).value = 0;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'freq centre [Hz]';
defaults.param(n).min = 0;
defaults.param(n).max = SAMPLE_RATE/2-1;
defaults.param(n).stepsize = 5;
defaults.param(n).value = 10400;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'bandwidth [Hz]';
defaults.param(n).min = 10;
defaults.param(n).max = SAMPLE_RATE/2-1;
defaults.param(n).stepsize = 25;
defaults.param(n).value = 20000;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'cos_ramp [ms]';
defaults.param(n).min = 0;
defaults.param(n).max = 500;
defaults.param(n).stepsize = 5;
defaults.param(n).value = 5;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'duration [ms]';
defaults.param(n).min = 0;
defaults.param(n).max = 20000;
defaults.param(n).stepsize = 100;
defaults.param(n).value = 1500;
defaults.param(n).ToolTip = 'No ToolTip';

defaults.frozen = 0;    % If non-zero: Reuse last wave form if no paramter change.


%==========================================================================
function wave = generate_waveform(specs, h_fig),
%%%%% Begin: DO NOT MODIFY ! Usually the same for all stimuli.
global SAMPLE_RATE   MAX_LEVEL, %Set by caller, otherwise by 'return_default_parameter'

for q = 1: length(specs.param),
	if exist(strtok(specs.param(q).name)),
		h = errordlg('Parametername dublication!',...
			mfilename); uiwait(h), wave = 0; return,
	end,
	eval([strtok(specs.param(q).name) '=' ...
		num2str(specs.param(q).value) ';']);
end,

if exist('spct_level'),
	amplitude(2) = 10^((spct_level - specs.param(1).max)/20); % That's why "spct_level is
	d = tanh(ILD./(20*log10(60)/.5/log(60)));       % always 1st control !
    amplitude(1) = amplitude(2)*(1-d);    
	amplitude(2) = amplitude(2)*(1+d);
	%%%%% END: DO NOT MODIFY !
else,
	% Don't name the level "spct_level" and define your own amplitude here.
end,

%//////////////////////////////////////////////////////////////////////////
% HERE WE GO: Specify here the code to generate the sound with the values in
% 'specs', a struct like returned by 'return_default_values' specified above
%//////////////////////////////////////////////////////////////////////////
% First check and if necessary transform the parameter
%//////////////////////////////////////////////////////////////////////////

ITD = ITD/1000000;
IPD = IPD/180*pi;

cos_ramp = round(cos_ramp/1000 * SAMPLE_RATE);
if cos_ramp == 1, cos_ramp = 0; end,
duration = round(duration/1000 * SAMPLE_RATE);
if (2*cos_ramp) > duration,
	h = errordlg('Wrong parameter settings: 2*cos_ramp > duration!', mfilename);
	uiwait(h), wave = 0; return,
end,

len = 2.^ceil(log2(duration)); % get the neares "2.^n" number of samples
if len < SAMPLE_RATE,
    len = 2.^ceil(log2(SAMPLE_RATE));
end,

f_hi = floor((freq - bandwidth/2)*len/SAMPLE_RATE)+1;
if f_hi < 0 |f_hi > len/2,
	h = errordlg('Wrong parameter setting: f_hi < 0 |f_hi > f2/2!',...
		mfilename); uiwait(h), wave = 0; return,
end,

f_lo = ceil((freq + bandwidth/2)*len/SAMPLE_RATE)+1;
if f_lo < f_hi || f_lo > len/2,
	h = errordlg('Wrong parameter setting: f_lo < f_hi |f_lo > fs/2!',...
		mfilename); uiwait(h), wave = 0; return,
end,

%//////////////////////////////////////////////////////////////////////////
% generate wave form %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%//////////////////////////////////////////////////////////////////////////

% create noise with coherence == 1
wave = [zeros(f_hi,1); 0.35355*sqrt(SAMPLE_RATE*len) ...
    * sqrt(-2*log(1-rand(f_lo - f_hi + 1,1))); zeros(len/2-f_lo,1)];
wave = [[0; 2*pi*rand(len/2-1,1); 0], wave]; % phase

% interaural phase difference
phase = wave(:,1) + (0:len/2)'/len*SAMPLE_RATE*2*pi*ITD + IPD;
[re im] = pol2cart(wave(:,1), wave(:,2));
wave(:,1) = re + i*im;
[re im] = pol2cart(phase, wave(:,2));
wave(:,2) = re + i*im;

% generate one channel of independent noise
[re im] = pol2cart([0 ; 2*pi*rand(len/2-1,1); 0], ...
    [zeros(f_hi,1); 0.35355*sqrt(SAMPLE_RATE*len) ...
    * sqrt(-2*log(1-rand(f_lo - f_hi + 1,1))); zeros(len/2-f_lo,1)]);

% mix left channel with independent noise to get required coherence
wave(:,1) = coherence*wave(:,1) + sqrt(1-coherence^2)*(re + i*im);

if pink_spctrm % shape to pink spectrum
    wave = 50* wave./repmat([1; sqrt([1:len/2]'/len*SAMPLE_RATE)],1,2);
end

% compensate amplitude for transfer function
if exist('h_fig','var')
    lsp_transfer_fcn = getappdata(h_fig, 'lsp_transfer_fcn');
    if ~isempty(lsp_transfer_fcn)
        f = linspace(0, SAMPLE_RATE/2,size(wave,1))';
        tf = interp1q([0:length(lsp_transfer_fcn)-1]',lsp_transfer_fcn,f);
        wave = wave./tf;
    end
end

% IFFT, then cut to the required length (duration) and apply envelope
wave = real(ifft([wave; conj(wave(len/2:-1:2,:))]));
wave = wave(1:duration,:) .* (amplitude' * ...
	[cos(pi/2:(pi-pi/2)/(cos_ramp-1):pi).^2, ...
	ones(1,duration - 2*cos_ramp),cos(0:pi/2/(cos_ramp-1):pi/2).^2])';

% IAC = corrcoef(wave);
% IAC = IAC(1,2),
% % Power: 
% p = sum(wave.^2)/length(wave)
