% [buffer_in, n_used] = DATA_ACQUISITION ( ...
%							wave_out, out_chs, in_chs, n_avg, use_crit)
%
% The $wave data are played $n_avg+2 times and first and the last loop 
% is discarded and therefore not averaged.
%
% Part of the OAE toolbox 
% Copyright (C) 2008 Torsten Marquardt
% Terms of the GNU General Public License apply
% (www.http://www.fsf.org/licensing/licenses/gpl.html).

function [buffer_in, perc_used] = data_acquisition( ...
    wave_out, out_chs, in_chs, n_avg, use_crit)
%%% Recommended settings, AR with Crit = 2
global SAMPLE_RATE,

if exist('use_crit','var')
    Crit = use_crit;
else
    Crit = 0; % 1: Mean + 2SD; 2: isoutlier; Values more than three scaled median absolute deviations (MAD) away from the median
end
if Crit <= 0
    AR = 0;
else
    AR = 1;   % 1: Apply artifact rejection; 0: Do not apply. Not using extra arguments.
end

if isempty(SAMPLE_RATE),
	SAMPLE_RATE = 48000;
end,

% wave_out = circshift(wave_out,1,2); % LF tone (old ch3) goes to ch 1! Rec1 = ch2, Rec2 = 3) !!! to be ompatible with Psycho
le = length(wave_out);
h_sound_device = sound_io('open_device',SAMPLE_RATE, out_chs, in_chs);
wave_in = sound_io('playrecord',h_sound_device,...
    repmat(wave_out,(n_avg+3),1),(n_avg+3)*le/SAMPLE_RATE);
sound_io('close_device',h_sound_device);

if AR == 0
    % wave_in = sound_io(SAMPLE_RATE, wave_out, out_chs, in_chs, n_avg+3,round(0.05 * SAMPLE_RATE));
    wave_in = reshape(wave_in(2*le+1:(n_avg+2)*le,:), le, n_avg,[]);
    buffer_in = double(squeeze(sum(wave_in, 2)./n_avg));
    perc_used = 100; % Percentage of non-rejected segments
elseif AR == 1 % Highpass, then get power, discard based on stats., and do weigthed average of remaining bocks
    fc = 250;                % Hz, highpass filter cutoff 
    Wn = fc/(SAMPLE_RATE/2); % Normalized cutoff
    b = fir1(100,Wn,'high');   % Xth-order (even nunber) highpass     
    wave_inHP = conv(wave_in,b,'same');                                   % Filter signal (while keeping same length). Avoid effect of L_TT power on calculations    
    wave_in = reshape(wave_in(2*le+1:(n_avg+2)*le,:), le, n_avg,[]);      % The original, unfiltered version    
    wave_in = wave_in';
    wave_inHP = reshape(wave_inHP(2*le+1:(n_avg+2)*le,:), le, n_avg,[]);  % Then reshape as above
    wave_inHP = wave_inHP';                                               % For ease transpose, then for "buffer_in" transpose again to match.    
    for i=1:size(wave_inHP,1)                                             % Get powers
        Powers(i) = sum(wave_inHP(i,:).^2)/le;
        W(i) = 1/Powers(i);                                               % Weigths        
    end
    if Crit == 1 % Mean + 2 SD
        LimitP = (Powers > mean(Powers) + 2*std(Powers));
    elseif Crit == 2 % Isoutlier. Ref.: "an outlier is a value that is more than three scaled median absolute deviations (MAD) away from the median."            
        MAD = 1.4826*median(abs(Powers-median(Powers))); % Median absolute deviations
        lowerbound = median(Powers) - 3*MAD; upperbound = median(Powers) + 3*MAD;
        LimitP = (Powers > upperbound | Powers < lowerbound); % Sane as isoutlier, median method
    end
    W = W((find(LimitP == 0)));
    W = W';                                                               % To apply weigthts on remaining segments
    wave_in = wave_in(find(LimitP == 0),:);                               % Power exclusion of segments. Exclusion on unfiltered version      
    for i = 1:length(W)
    Aux(i,:) = wave_in(i,:)*W(i);
    end
    buffer_in = (sum(Aux)/sum(W))';                                       % Weithed average. On kept segments. 
    % buffer_in = mean(wave_in)'; % For simple average
    wave_in = wave_in';
    % Matlab version issue with mult!, does not allow buffer_in =
    % (sum(wave_in.*W)/sum(W))'; !!!      
    perc_used = (length(W)/n_avg)*100; % Percentage of non-rejected segments     
   % n_avg = length(W); % update, not to give error below
end

% if nargin > 4, % artefact suppression
% 	used = mean((wave_in - repmat(buffer_in,1,n_avg)).^2)*1e6 < use_crit;
% 	n_used = sum(used,2);
% end,
