function varargout = CAP_thresholds(varargin)
% CAP_thresholds(thrshld,n_repetitions,rep_rate,start_levels,f_l,f_u,slopes)
% thrshld is given in no. of stddevs above current noise floor
% Number of $start_levels defines n_freqs, but will be rounded up to next
% multiple of 5.
% $slopes gives slope assumption at the beginning and will later be adjusted

% To do:
% - show average CAP shapes and CMs for responses above threshold
%   (normalised using estimated slope?)
% - browse/overlay with previos CAP thesholds
% - in broweser: select for average CAP thresholds

if nargin == 0
    start
elseif (nargout)
    [varargout{1:nargout}] = feval(varargin{1},varargin{2:end});
else
    feval(varargin{1}, varargin{2:end}); % FEVAL switchyard
end

%=========================================================================
function start

fs = 48000;

CAP_window = 1800:2300;
results.cos_ramp_cycles = 12; % cycles

if ~exist('n_presentations','var')|| isempty(n_repetitions)
    n_repetitions = 10;
end

if ~exist('rep_rate','var')|| isempty(rep_rate)
    rep_rate = 20; % Hz
end

if ~exist('start_levels','var')|| isempty(start_levels) % loading of standard parameter file might be overidden
    if exist('CAP_threshold_start_param','file')
        try % load standard parameter file
            load('CAP_threshold_start_param','results')
        catch
            warning('Could not load parameter file. Using default!')
        end
    else
        start_levels = 50*ones(1,32);
    end
end
if ~exist('freqs','var') || ~exist('f_l','var') || ~exist('results','var') 
    f_l = 2000;
    f_u = 11986; % 19027;
    results.freqs = calcualte_freqs(f_l,f_u,size(start_levels));% freqs are in matrix to play in staggered order
end
n_freqs = length(results.freqs);
if ~exist('results','var')||~exist('slopes','var')|| isempty(slopes)
    results.slopes = 3*ones(n_freqs,1);
end

load('C:\Users\Admin\Documents\MATLAB\IDA_asio\DP_TT\_last_probe_check')
Hdb_2 = 20*log10(abs(CURRENT_EAR(round(results.freqs/10),2)))';


tone_pips = generate_tone_pips(results.freqs,rep_rate,results.cos_ramp_cycles,fs); % generate tone pips of all frequencies

bp_filter_settings = designfilt('bandpassiir','FilterOrder',20,...
    'PassbandFrequency1',350,'PassbandFrequency2',3000,'StopbandAttenuation1',60,...
    'PassbandRipple',1,'StopbandAttenuation2',80,'SampleRate',fs);


% record_noise_floor
h_sound_device = sound_io('open_device',fs, 2, 2);
wave_in = sound_io('playrecord',h_sound_device,zeros(3*fs,2),2);
results.stddev_noise = 2*std(filtfilt(bp_filter_settings,wave_in(:,2)));
% results.stddev_noise = 2*std(wave_in(:,2)); %XXX

n_wave_out = length(tone_pips)/2;
h_sound_device = sound_io('open_device',fs, 2, 2);
PsychPortAudio('GetAudioData', h_sound_device, 4*n_wave_out/fs);
PsychPortAudio('FillBuffer',h_sound_device,zeros(2,2*n_wave_out));
PsychPortAudio('Start', h_sound_device);

[CAP_threshold_fig,h_PB_pause,h_PB_stop,h_EF_counter] = init_CAP_threshold_fig(...
    results.stddev_noise,results.freqs);

results.CAP_ampls = NaN * ones(n_repetitions+1,n_freqs);
results.levels = NaN * ones(n_repetitions+1,n_freqs);
results.levels(1,:) = start_levels;
results.valid = zeros(n_repetitions+1,n_freqs);
results.SPL_noise(1:n_freqs) = NaN;
CAP_shapes_1 = zeros(length(CAP_window),n_freqs/2);
CAP_shapes_2 = zeros(length(CAP_window),n_freqs/2);
loop_time=tic;
for q = 1:n_repetitions+1
    set(h_EF_counter,'Value',['Rpt:' num2str(q-1), '/' num2str(n_repetitions)])

    % check pause / continue
    if ishandle(h_PB_pause)
        if strcmp(get(h_PB_pause,'text'),'Continue')
            waitfor(h_PB_pause,'text')
        end
        if strcmp(get(h_PB_stop,'text'),'Stopped')
            ida_exp_gui('CloseRequestFcn')
            break
        end
    else
        CAP_thresholds('CloseRequestFcn')
    end
    
    % Play first half, analyse second half
    ampl = 10.^((results.levels(q,1:n_freqs/2)-Hdb_2(1:n_freqs/2)-CURRENT_EAR(2,5))./20);
    ampl_vector = ones(n_wave_out/(n_freqs/2),1)*ampl;
    wave_out =  ampl_vector(1:end).*tone_pips(1:end/2);
    if abs(wave_out)>1
        warndlg('Max of summed wave > 1!')
        wave_out = wave_out/max(wave_out);
    end

    if ishandle(CAP_threshold_fig) % if CAP_threshold_fig still open
        wave_in = PsychPortAudio('GetAudioData',h_sound_device,[],2*n_wave_out/fs,2*n_wave_out/fs)';        
        wave_in = wave_in(1:2*n_wave_out,2);
        if q <= n_repetitions    % present first half of wave forms
            PsychPortAudio('RefillBuffer',h_sound_device,0,[zeros(1,2*n_wave_out);[wave_out -wave_out]]);
            PsychPortAudio('Start', h_sound_device);
            true_presentation_period = toc(loop_time);
            loop_time =tic;
            if true_presentation_period > 2*n_wave_out/fs +0.05 % allow 50 ms overun
                warning('ida:Presentation_period','Presentation period %f s larger than specified by',...
                    true_presentation_period - n_wave_out/fs)
            end
        end
    end

    if q>1 % analyse second half
        RW_signal_filtered = filtfilt(bp_filter_settings,wave_in);

%         RW_signal_filtered = wave_in; %XXX

        [CAP_ampls,noise_levels,CAP_waves] = extract_CAP_amplitude(RW_signal_filtered,CAP_window,n_freqs);
        valid_CAPshape_ampl = 3*results.stddev_noise;
        CAP_shapes_2(:,CAP_ampls>valid_CAPshape_ampl) = CAP_shapes_2(:,CAP_ampls>valid_CAPshape_ampl)+CAP_waves(:,CAP_ampls>valid_CAPshape_ampl);
        results.CAP_ampls(q-1,n_freqs/2+1:end) = CAP_ampls;
        results = calcualte_next_levels(q-1,0,results,noise_levels);
        
        update_CAP_threshold_fig(CAP_threshold_fig,q,results,[CAP_shapes_1 CAP_shapes_2])
    end

    if q <= n_repetitions    % Play second half, analyse first half
        ampl = 10.^((results.levels(q,n_freqs/2+1:end)-Hdb_2(n_freqs/2+1:end)-CURRENT_EAR(2,5))./20);
        ampl_vector = ones(n_wave_out/(n_freqs/2),1)*ampl;
        wave_out =  ampl_vector(1:end).*tone_pips(end/2+1:end);
        if abs(wave_out)>1
            warndlg('Max of summed wave > 1!')
            wave_out = wave_out/max(wave_out);
        end

        if ishandle(CAP_threshold_fig) % if CAP_threshold_fig still open
            wave_in = PsychPortAudio('GetAudioData',h_sound_device,[],2*n_wave_out/fs,2*n_wave_out/fs)';
            wave_in = wave_in(1:2*n_wave_out,2);
            % present second half of wave forms
            PsychPortAudio('RefillBuffer',h_sound_device,0,[zeros(1,2*n_wave_out);[wave_out -wave_out]]);
            PsychPortAudio('Start', h_sound_device);
            true_presentation_period = toc(loop_time);
            loop_time =tic;
            if true_presentation_period > 2*n_wave_out/fs +0.05 % allow 50 ms overun
                warning('ida:Presentation_period','Presentation period %f s larger than specified by',...
                    true_presentation_period - n_wave_out/fs)
            end
        end


        % analyse first half
        RW_signal_filtered = filtfilt(bp_filter_settings,wave_in);

        % RW_signal_filtered = wave_in; %XXX

        [CAP_ampls,noise_levels,CAP_waves] = extract_CAP_amplitude(RW_signal_filtered,CAP_window,n_freqs);
        valid_CAPshape_ampl = 3*results.stddev_noise;
        CAP_shapes_1(:,CAP_ampls>valid_CAPshape_ampl) = CAP_shapes_1(:,CAP_ampls>valid_CAPshape_ampl)+CAP_waves(:,CAP_ampls>valid_CAPshape_ampl);
        results.CAP_ampls(q,1:n_freqs/2) = CAP_ampls;
        results = calcualte_next_levels(q,1,results,noise_levels);

        if strcmp(get(h_PB_stop,'text'),'Stopped')
            break
        end
    end
end
PsychPortAudio('Close')
results.CAP_shapes = [CAP_shapes_1 CAP_shapes_2];
setappdata(CAP_threshold_fig,'results',results)
set(CAP_threshold_fig,'Color',[1 1 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function freqs = calcualte_freqs(f_l,f_u,freqs_size)
% freqs are in matrix to play in staggered order row-by-row
tmp = reshape(round(logspace(log10(f_u),log10(f_l),prod(freqs_size))),[],4)';
freqs = tmp(1:end);

%%========================================================================
function tone_pips = generate_tone_pips(freqs,rep_rate,cos_ramp_cycles,fs)
n_freqs = length(freqs);
total_duration = 1/rep_rate*fs;
tone_pips = zeros(1,n_freqs*total_duration);
for q = 0:n_freqs-1
    l_cos_ramp = round(cos_ramp_cycles/freqs(q+1)*fs);
    duration= 2*l_cos_ramp+480; 
    tone_pips(q*total_duration+1:q*total_duration+total_duration) = [...
        sin(linspace(0,freqs(q+1)/fs*duration*2*pi*(1-1/duration),duration))'... tone
        .* [cos(pi/2:(pi-pi/2)/(l_cos_ramp-1):pi).^2,ones(1,480),cos(0:pi/2/(l_cos_ramp-1):pi/2).^2]';... enveleope
        zeros(total_duration-duration,1)]; % trailing silence
end

%%========================================================================
function [CAP_ampl, noise_levels,CAP_waves] = extract_CAP_amplitude(...
    RW_signal_filtered,CAP_window,n_freqs)
n = length(RW_signal_filtered)/2;
CAP_waves = reshape(RW_signal_filtered(1:n) + RW_signal_filtered(n+1:end),[],n_freqs/2); 
% CAP_waves = reshape(RW_signal_filtered(1:n) - RW_signal_filtered(n+1:end),[],n_freqs/2); %XXX
% take out interstimulus intervals for noise estmation (stddev_noise) and remainder for CAP_wave
noise_window = ones(size(CAP_waves,1),1);
noise_window(CAP_window) = 0;
noise_levels = std(CAP_waves(noise_window>0,:)); 
CAP_waves = CAP_waves(CAP_window,:);
CAP_ampl = max(CAP_waves) - min(CAP_waves);

%%========================================================================
function results = calcualte_next_levels(q,first_half,results,noise_levels)

n_freqs = size(results.levels,2);
if first_half
    freqs = 1:n_freqs/2;
    CAP_ampls_in_stddev = results.CAP_ampls(q,1:n_freqs/2)/results.stddev_noise;
else
    freqs = n_freqs/2+1:n_freqs;
    noise_levels = [zeros(1,n_freqs/2) noise_levels];
    CAP_ampls_in_stddev = [zeros(1,n_freqs/2) results.CAP_ampls(q,n_freqs/2+1:n_freqs)/results.stddev_noise];
end

q
for freq=freqs
    if (noise_levels(freq) > 2*results.stddev_noise) && (results.CAP_ampls(q,freq) < 2*noise_levels(freq)) % invalid (too noisy)
        results.levels(q+1,freq) = results.levels(q,freq); % do again
        results.SPL_noise(freq) = results.levels(q,freq);
        results.valid(q,freq) = 0;
        1
    elseif results.CAP_ampls(q,freq) < 2* results.stddev_noise % CAP below noise floor
        results.levels(q+1,freq) = results.levels(q,freq)+10;
        results.SPL_noise(freq) = max([results.SPL_noise(freq) results.levels(q,freq)]);
        results.valid(q,freq) = 0;
        2
    elseif q<10
        results.levels(q+1,freq) = results.levels(q,freq)-results.slopes(freq)*(CAP_ampls_in_stddev(freq)-3);
%     elseif (q==1) || sum(results.valid(1:q-1,freq)) == 0 % CAP valid the first time
%         if isnan(results.SPL_noise(freq)) % no SPL estimate yet that achives CAP at noise floor level
%             results.levels(q+1,freq) = results.levels(q,freq)-results.slopes(freq)*(CAP_ampls_in_stddev(freq)-3);
%             results.valid(q,freq) = 1;
%             3
%         else % SPL estimate of noise and new valid CAP allows estimate of slope
%             results.slopes(freq) = (results.levels(q,freq)-results.SPL_noise(freq))/(CAP_ampls_in_stddev(freq)-results.stddev_noise);
%             results.levels(q+1,freq) = results.levels(q,freq) - results.slopes(freq) * (CAP_ampls_in_stddev(freq)-3);
%             results.valid(q,freq) = 1;
%             4
%         end
    else
        results.valid(q,freq) = 1;
        valid_levels = results.levels(results.valid(:,freq)==1,freq);
        valid_CAP_ampls_in_stddev = results.CAP_ampls(results.valid(:,freq)==1,freq)/results.stddev_noise;
        b = [ones(length(valid_levels),1) valid_CAP_ampls_in_stddev]\valid_levels;
        results.slopes(freq) = b(2);
        results.levels(q+1,freq) = b(1) + b(2)*3;
        results.SPL_noise(freq) = b(1) + b(2);
        5
    end
end
results.levels(q+1,:)
q;
%%========================================================================
function update_CAP_threshold_fig(h_fig,q,results,CAP_shapes)
figure(h_fig)

[freqs,f_idx] = sort(results.freqs);

axes(findobj(h_fig,'Tag','CAP_shape_axes'));
CAP_shapes = CAP_shapes(:,f_idx);
CAP_shapes = CAP_shapes./repmat(max(CAP_shapes),size(CAP_shapes,1),1);
plot(CAP_shapes(1:end))
set(gca,'Tag','CAP_shape_axes')
axis([1 numel(CAP_shapes) -1 1])

axes(findobj(h_fig,'Tag','threshold_axes'));
hold off
semilogx(freqs,results.SPL_noise(f_idx),'k-')
hold on
semilogx(freqs,results.levels(q,f_idx),'bx-')

set(gca,'Tag','threshold_axes')
axis([2000 20000 -10 60])
pause(0.1)

%%========================================================================
function [h_fig,h_PB_pause,h_PB_stop,h_EF_counter]= init_CAP_threshold_fig(stddev_noise,freqs)
h_fig = findobj('Tag','CAP_threshold_fig');
if isempty(h_fig)
    h_fig = uifigure;
else
    figure(h_fig);
    clf
end
h_fig.HandleVisibility = 'on';
h_fig.Tag = 'CAP_threshold_fig';
h_fig.Position = [0 60 900 600];
h_fig.Resize = 1;
h_fig.CloseRequestFcn  = 'CAP_thresholds CloseRequestFcn';
h_fig.NumberTitle = 'on';
set(h_fig,'Color',[.8 .8 .8])

h_ida_fig = findobj('Tag','ida_fig');
if ishandle(h_ida_fig)
    animalID = getappdata(h_ida_fig,'animalID');
else
    animalID = 'tst';
end
h_fig.Name = [datestr(now,'yyyymmddTHHMMSS'),'-',animalID,...
    '-CAPthrlds'];

h_EF_counter = uieditfield(h_fig, ...
    'Value','counter', ...
    'FontSize',9, ...
    'Editable','off', ...
    'Tag','txt_counter', ...
    'Position',[5 5 50 25]);
h_PB_pause = uibutton(h_fig, ...
    'text','pause', ...
    'ButtonPushedFcn',@CB_pause_button, ...
    'Tag','PB_pause', ...
    'Position',[65 5 55 25]);
h_PB_stop = uibutton(h_fig, ...
    'text','stop', ...
    'ButtonPushedFcn',@CB_stop_button, ...
    'Tag','PB_stop', ...
    'Position',[125 5 55 25]);
uibutton(h_fig, ...
    'text','save', ...
    'ButtonPushedFcn',@CB_save_button, ...
    'Tag','PB_save', ...
    'Position',[185 5 55 25]);

h_axes = axes;
set(gca,'Position',[.08 0.88 0.9 0.10])
semilogx([1000 30000],[.5 .5])
set(h_axes,'Tag','CAP_shape_axes')
axis([1000 30000 0 1])

h_axes = axes;
set(h_axes,'Tag','threshold_axes')
set(gca,'Position',[.08 0.11 0.9 0.765])
semilogx([freqs(end,end)/1.05 freqs(1,1)*1.05],[0 0])
set(h_axes,'Tag','threshold_axes')
axis([1000 30000 -10 60])
pause(.1)

%%========================================================================
function CB_pause_button(callingObj,~)
if strcmp(get(callingObj,'text'),'Pause')
    set(callingObj,'text','Continue')
else
    set(callingObj,'text','Pause')
end
%% =========================================================================
function CB_stop_button(callingObj,~)
set(callingObj,'text','Stopped')

%%========================================================================
function CB_save_button
% When pressing save button in CAP_threshold_fig save data and close figure
global DATA_PATH

CAP_threshold_fig = gcbf;
savefig(CAP_threshold_fig,[DATA_PATH '\' CAP_threshold_fig.Name])

%% =========================================================================
function CloseRequestFcn
h_fig =findobj('Tag','CAP_threshold_fig');
delete(h_fig(end)),
