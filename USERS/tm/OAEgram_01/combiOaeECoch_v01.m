function output = stimulus(varargin)
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
% Date: 03-03-22
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

function defaults = return_default_parameter(vers)
% specify variables and default values of sound.
% Stepsize is for the slider in the sound figure

defaults.name = [mfilename ' ' vers]; % DO NOT MODIFY!
n = 0; % DO NOT MODIFY!

global SAMPLE_RATE   MAX_LEVEL, % if not set by caller set default values
if isempty(SAMPLE_RATE), SAMPLE_RATE = 48000; end
if isempty(MAX_LEVEL), MAX_LEVEL = 110; end % dependant on sound card
% A tone at MAX_LEVEL has amplitude of one.

% !!! KEEP spectrum level always as the 1st parameter (for plot function)!

n = n+1;
defaults.param(n).name = 'l_BT [dB]';     
defaults.param(n).min = -40;                    
defaults.param(n).max = 105;             
defaults.param(n).stepsize = 5;                 
defaults.param(n).value = 80;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'l_suppr [+dB]';     
defaults.param(n).min = 0;                    
defaults.param(n).max = 30;             
defaults.param(n).stepsize = 5;                 
defaults.param(n).value = 10;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'l2 [dB]';     
defaults.param(n).min = -40;                    
defaults.param(n).max = 111;             
defaults.param(n).stepsize = 5;                 
defaults.param(n).value = 45;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'l1 [+dB]';     
defaults.param(n).min = -40;                    
defaults.param(n).max = 40;             
defaults.param(n).stepsize = 5;                 
defaults.param(n).value = 10;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'f_BT [Hz]';
defaults.param(n).min = 20;
defaults.param(n).max = SAMPLE_RATE/2;
defaults.param(n).stepsize = -2^(1/3);
defaults.param(n).value = 60;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'f2 [Hz]';
defaults.param(n).min = 30;
defaults.param(n).max = SAMPLE_RATE/2;
defaults.param(n).stepsize = -2^(1/3);
defaults.param(n).value = 6010;
defaults.param(n).ToolTip = 'No ToolTip';

n = n+1;
defaults.param(n).name = 'f1 [Hz]';
defaults.param(n).min = 30;
defaults.param(n).max = SAMPLE_RATE/2;
defaults.param(n).stepsize = -2^(1/3);
defaults.param(n).value = 5010;
defaults.param(n).ToolTip = 'No ToolTip';

defaults.frozen = 1;    % If non-zero: Reuse last wave form if no paramter change.

%==========================================================================
function wave = generate_waveform(specs, h_fig)

%% Begin: DO NOT MODIFY ! Usually the same for all stimuli.
global SAMPLE_RATE %Set by caller, otherwise by 'return_default_parameter'

fs = SAMPLE_RATE;
cos_ramp = 3; % cycles
for q = 1: length(specs.param)
	eval([strtok(specs.param(q).name) '=' ...
		num2str(specs.param(q).value) ';']);
end
%% End: DO NOT MODIFY ! Usually the same for all stimuli.
%//////////////////////////////////////////////////////////////////////////
% HERE WE GO: Specify here the code to generate the sound with the values in
% 'specs', a struct like returned by 'return_default_values' specified above
%//////////////////////////////////////////////////////////////////////////
% First check and if necessary transform the parameter
%//////////////////////////////////////////////////////////////////////////
if any(f2 > fs/2) || any(f2 <= 0)
	h = errordlg('Wrong parameter setting: a frequency > fs/2, or <= 0!',...
		mfilename); uiwait(h),  return,
end

f_suppr1 = f1 - 40;
f_suppr2 = f2 - 40;
f_BTdivFs = f_BT/fs;
l1 = l2 + l1;
l_suppr = l2 + l_suppr;

CURRENT_EAR = getappdata(h_fig, 'CURRENT_EAR');
%  CURRENT_EAR(:,1) - rec1 transfer fct. (complex)
%  CURRENT_EAR(:,2) - rec2 transfer fct. (complex)
%  CURRENT_EAR(:,3) - TT   transfer fct. (complex)
%  CURRENT_EAR(:,4) - mic transfer fct. (complex)
%  CURRENT_EAR(1,5) - max. soundpressure level of rec1 [dB SPL]
%  CURRENT_EAR(2,5) - max. soundpressure level of rec2 [dB SPL]
%  CURRENT_EAR(3,5) - max. soundpressure level of TT [dB SPL]

% phase correction
phase_dev_1 = angle(CURRENT_EAR(f1/10+1,1));
phase_dev_2 = angle(CURRENT_EAR(f2/10+1,2));
phase_dev_BT = angle(CURRENT_EAR(f_BT/10+1,3));
Hdb_1 = 20*log10(abs(CURRENT_EAR(f1/10+1,1)'));
Hdb_2 = 20*log10(abs(CURRENT_EAR(f2/10+1,2)'));
Hdb_BT = 20*log10(abs(CURRENT_EAR(f_BT/10+1,3)'));
Hdb_suppr1 = 20*log10(abs(CURRENT_EAR(f_suppr1/10+1,1)'));
Hdb_suppr2 = 20*log10(abs(CURRENT_EAR(f_suppr2/10+1,1)'));

if max(l1 - Hdb_1) > CURRENT_EAR(1,5)
    msgbox('The corrected Level L1 is higher than max_level!',...
        'ERROR - DP_TT Measurement:','error');
    return,
end
if max(l_BT - Hdb_BT) > CURRENT_EAR(1,5)
    msgbox('The corrected Level L_TT is higher than max_level!',...
        'ERROR - DP_TT Measurement:','error');
    return,
end
if max(l_suppr - Hdb_suppr2) > CURRENT_EAR(1,5)
    msgbox('The corrected Level L_suppr is higher than max_level!',...
        'ERROR - DP_TT Measurement:','error');
    return,
end
if max(l2 - Hdb_2) > CURRENT_EAR(2,5)
    msgbox('The corrected Level L2 is higher than max_level!',...
        'ERROR - DP_TT Measurement:','error');
    return,
end

ampl_1 = 10^((l1-Hdb_1-CURRENT_EAR(1,5))/20);
ampl_suppr = 10^((l_suppr-Hdb_suppr2-CURRENT_EAR(1,5))/20);
ampl_BT = 10^((l_BT-Hdb_BT-CURRENT_EAR(1,5))/20);
ampl_2 = 10^((l2-Hdb_2-CURRENT_EAR(2,5))/20);

ramp_on = cos(pi/2:(pi-pi/2)/(cos_ramp/f2*fs-1):pi).^2;
ramp_off = cos(0:pi/2/(cos_ramp/f2*fs-1):pi/2).^2;
ramp_len = length(ramp_on);
BT_ramp_on = cos(pi/2:(pi-pi/2)/(9/f_BTdivFs-1):pi).^2;
BT_ramp_off = cos(0:pi/2/(9/f_BTdivFs-1):pi/2).^2;

duration = ones(1,6.5/f_BTdivFs - 2*ramp_len);
duration_2 = ones(1,7/f_BTdivFs - 2*ramp_len);

duration_total = (6+1+32+1)*9/f_BTdivFs;
tone_1 = ampl_1*sin(linspace(0,f1/fs*duration_total*2*pi*(1-1/duration_total),duration_total) - phase_dev_1);
tone_suppr = ampl_suppr*sin(linspace(0,f_suppr2/fs*duration_total*2*pi*(1-1/duration_total),duration_total));
tone_BT = ampl_BT*sin(linspace(0,f_BT/fs*duration_total*2*pi*(1-1/duration_total),duration_total) - phase_dev_BT);
tone_2 = ampl_2*sin(linspace(0,f2/fs*duration_total*2*pi*(1-1/duration_total),duration_total) - phase_dev_2);

env_1 = [zeros(1,2.5/f_BTdivFs) zeros(1,6.5/f_BTdivFs) -zeros(1,2.5/f_BTdivFs) -zeros(1,6.5/f_BTdivFs) ...
    zeros(1,2.5/f_BTdivFs) zeros(1,6.5/f_BTdivFs) -zeros(1,2.5/f_BTdivFs) -zeros(1,6.5/f_BTdivFs) ...
    zeros(1,2.5/f_BTdivFs) ramp_on duration ramp_off -zeros(1,2.5/f_BTdivFs) -ramp_on -duration -ramp_off ...
    ...
    zeros(1,9/f_BTdivFs) ...
    repmat([zeros(1,2.5/f_BTdivFs) zeros(1,6.5/f_BTdivFs) -zeros(1,2.5/f_BTdivFs) -zeros(1,6.5/f_BTdivFs) ...
    zeros(1,2.5/f_BTdivFs) ramp_on duration ramp_off -zeros(1,2.5/f_BTdivFs) -ramp_on -duration -ramp_off],1,8)...
    zeros(1,9/f_BTdivFs)];

env_2 = [repmat([zeros(1,2/f_BTdivFs) ramp_on duration_2 ramp_off -zeros(1,2/f_BTdivFs) -ramp_on -duration_2 -ramp_off],1,3) ...
    ...
    zeros(1,9/f_BTdivFs) ...
    repmat([...
    zeros(1,1/f_BTdivFs) zeros(1,1/8/f_BTdivFs) ramp_on ones(1,7/8/f_BTdivFs) duration_2 ramp_off ...
    -zeros(1,1/f_BTdivFs) -zeros(1,1/8/f_BTdivFs) -ramp_on -ones(1,7/8/f_BTdivFs) -duration_2 -ramp_off ...
    ...
    zeros(1,1/f_BTdivFs) zeros(1,2/8/f_BTdivFs) ramp_on ones(1,6/8/f_BTdivFs) duration_2 ramp_off ...
    -zeros(1,1/f_BTdivFs) -zeros(1,2/8/f_BTdivFs) -ramp_on -ones(1,6/8/f_BTdivFs) -duration_2 -ramp_off ...
    ...
    zeros(1,1/f_BTdivFs) zeros(1,3/8/f_BTdivFs) ramp_on ones(1,5/8/f_BTdivFs) duration_2 ramp_off ...
    -zeros(1,1/f_BTdivFs) -zeros(1,3/8/f_BTdivFs) -ramp_on -ones(1,5/8/f_BTdivFs) -duration_2 -ramp_off ...
    ...
    zeros(1,1/f_BTdivFs) zeros(1,4/8/f_BTdivFs) ramp_on ones(1,4/8/f_BTdivFs) duration_2 ramp_off ...
    -zeros(1,1/f_BTdivFs) -zeros(1,4/8/f_BTdivFs) -ramp_on -ones(1,4/8/f_BTdivFs) -duration_2 -ramp_off ...
    ...
    zeros(1,1/f_BTdivFs) zeros(1,5/8/f_BTdivFs) ramp_on ones(1,3/8/f_BTdivFs) duration_2 ramp_off ...
    -zeros(1,1/f_BTdivFs) -zeros(1,5/8/f_BTdivFs) -ramp_on -ones(1,3/8/f_BTdivFs) -duration_2 -ramp_off ...
    ...
    zeros(1,1/f_BTdivFs) zeros(1,6/8/f_BTdivFs) ramp_on ones(1,2/8/f_BTdivFs) duration_2 ramp_off ...
    -zeros(1,1/f_BTdivFs) -zeros(1,6/8/f_BTdivFs) -ramp_on -ones(1,2/8/f_BTdivFs) -duration_2 -ramp_off ...
    ...
    zeros(1,1/f_BTdivFs) zeros(1,7/8/f_BTdivFs) ramp_on ones(1,1/8/f_BTdivFs) duration_2 ramp_off ...
    -zeros(1,1/f_BTdivFs) -zeros(1,7/8/f_BTdivFs) -ramp_on -ones(1,1/8/f_BTdivFs) -duration_2 -ramp_off ...
    ...
    zeros(1,2/f_BTdivFs) ramp_on duration_2 ramp_off ...
    -zeros(1,2/f_BTdivFs) -ramp_on  -duration_2 -ramp_off ...
    ],1,2) ...
    zeros(1,9/f_BTdivFs)];

env_suppr = [zeros(1,2.5/f_BTdivFs) zeros(1,6.5/f_BTdivFs) -zeros(1,2.5/f_BTdivFs) -zeros(1,6.5/f_BTdivFs) ...
    zeros(1,2.5/f_BTdivFs) ramp_on duration ramp_off -zeros(1,2.5/f_BTdivFs) -ramp_on -duration -ramp_off ...
    zeros(1,2.5/f_BTdivFs) zeros(1,6.5/f_BTdivFs) -zeros(1,2.5/f_BTdivFs) -zeros(1,6.5/f_BTdivFs) ...
    ...
    zeros(1,round((1+32+1)*9/f_BTdivFs))];

env_BT = [zeros(1,(6*9)/f_BTdivFs) ...
    BT_ramp_on  ones(1,32*9/f_BTdivFs) BT_ramp_off];

wave(:,1) = [env_1.*tone_1 + env_suppr.*tone_suppr + env_BT.*tone_BT]';
wave(:,2) = [env_2.*tone_2]';
% wave(:,3) = [env_BT.*tone_BT]';

% plot(linspace(0,6,288000), env_1), hold on
% plot(linspace(0,6,288000), env_2*.9), 
% plot(linspace(0,6,288000), env_suppr*.8), 
% plot(linspace(0,6,288000), env_BT*1.1), hold off
% set(gca,'XTick',.05+[0:.15:6]), grid on
% 
% plot(env_1), hold on
% plot(env_2*.9), 
% plot(env_suppr*.8), 
% plot(env_BT*1.1), hold off

