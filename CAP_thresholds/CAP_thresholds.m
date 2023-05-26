function CAP_thresholds(thrshld,n_repetitions,rep_rate,start_levels,f_l,f_u,slopes)
% thrshld is given in no. of stddevs above current noise floor
% Number of $start_levels defines n_freqs, but will be rounded up to next
% multiple of 5.
% $slopes gives slope assumption at the beginning and will later be adjusted

% To do:
% - incorporate lsp_transfer_function.mat in generate_tone_pips() 

global SAMPLE_RATE 

if ~exist('n_presentations','var')|| isempty(n_repetitions)
    n_repetitions = 30;
end

if ~exist('rep_rate','var')|| isempty(rep_rate)
    rep_rate = 20; % Hz
end
    
if ~exist('start_levels','var')|| isempty(start_levels) % loading of standard parameter file might be overidden
    try % load standard parameter file 
        load('start_levels','start_levels','freqs','slopes')
    catch
        warning('Could not load start levels. Using default!')
        start_levels = zeros(1,40);
    end
end
if ~exist('freqs','var') || ~exist('f_l','var')
    f_l = 2000;
    f_u = 19027;
    freqs = calcualte_freqs(f_l,f_u,length(start_levels));% freqs are in matrix to play in staggered order
end
if ~exist('slopes','var')|| isempty(slopes)
    slopes = 1e-3;
end
if length(slopes) == 1 % starting slopes assumption
    slopes = slopes*ones(1,numel(start_levels));
end


load('C:\Users\Admin\Documents\MATLAB\IDA_asio\DP_TT\_last_probe_check')
tone_pips = generate_tone_pips(freqs,rep_rate,SAMPLE_RATE); % generate tone pips of all frequencies
bp_filter_settings = designfilt('bandpassiir','FilterOrder',20,...
    'PassbandFrequency1',350,'PassbandFrequency2',3000,'StopbandAttenuation1',60,...
    'PassbandRipple',1,'StopbandAttenuation2',80,'SampleRate',SAMPLE_RATE);

PsychPortAudio('Close')
h_sound_device = sound_io('open_device',SAMPLE_RATE, 2, 2);

% record_noise_floor
wave_in = sound_io('playrecord',h_sound_device,zeros(3*SAMPLE_RATE,2),2);
stdev_noise = stdev(filtfilt(bp_filter_settings,wave_in(:,2)));

if ~exist('thrshld','var')|| isempty(thrshld)
    thrshld = 3*stdev_noise;
end

CAP_ampl = NaN * ones(size(freqs,1),size(freqs,2),n_repetitions);
levels = CAP_ampl;
levels(:,1) = start_levels';

for q=1:n_repetitions
    for q2=1:numel(tone_pips)
        freq_group_idx_next = rem(q*q2,n_repetitions)+1;
        repetition_next = floor(q*q2/n_repetitions)+1;
        Hdb_2 = 20*log10(abs(CURRENT_EAR(round(fregs(q)/10)+1,2)'));
        ampl = 10^((levels(freq_group_idx_next,repetition_next,:)-Hdb_2-CURRENT_EAR(2,5))/20);
        wave_out = ampl*tone_pips(freq_group_idx_next);
        pause(1/rep_rate - toc)
        tic
        wave_in = PsychPortAudio(wave_out);
        RW_signal_filtered = filtfilt(bp_filter_settings,wave_in(:,2));
        [CAP_ampl(q2,q), noise_level] = extract_CAP_amplitude(RW_signal_filtered);
        [levels(q2,q+1),slopes(q2)] = calcualte_next_levels(...
            q,levels(q2,:),CAP_ampl(q2,:),slopes(q2),noise_level,thrshld); 

        update_CAP_threshold_fig(CAP_threshold_fig,q,q2,levels,CAP_ampl,slopes,noise_level)
    end
end

PsychPortAudio('GetAudioData', h_sound_device, 2*stimulus.presentation_period);
PsychPortAudio('FillBuffer',h_sound_device,zeros(4,stimulus.presentation_period*fs));
PsychPortAudio('Start', h_sound_device);

% play n_rpts times a random sequences of the binaural intervals
loop_time=tic;
no_presentations=(stimulus.n_rpts)*(size(parameter_table,1));
first_presentation = 1; interval_no = 0; stimGenError = 0;
for n = (q-1)*size(parameter_table,1)+q2:no_presentations+1
    if n == 2, initSeriesFig = 1; else, initSeriesFig=0; end
    prmSetToSave = interval_no;
    q = floor((n-1)/(size(parameter_table,1)))+1;
    q2 = mod(n-1,size(parameter_table,1))+1;
    save([data_path, header.title], 'q','q2','-append')

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
        ida_exp_gui('CloseRequestFcn')
    end

    if n <= no_presentations
        interval_no = stimulus.interval_order(q,q2);
        [q q2 interval_no]; % print in command window


        % generate interval ID channel (#1) waveform
        intervalID_wave = create_intervalID_wave(fs,...
            parameter_table(interval_no,1),q2);
        set(findobj(CAP_threshold_fig,'Tag','txt_counter'),'Value',...
            ['Rpt:' num2str(q), '/' num2str(stimulus.n_rpts)])
    end

    % present wave form
    if ishandle(CAP_threshold_fig) % if CAP_threshold_fig still open
        wave_in = PsychPortAudio('GetAudioData',h_sound_device, [],...
            stimulus.presentation_period,stimulus.presentation_period)';
        if n <= no_presentations && ~stimGenError
            PsychPortAudio('RefillBuffer',h_sound_device,0,[wave'; ...
                [intervalID_wave zeros(1,size(wave,1)-length(intervalID_wave))]; ...
                [spike_wave zeros(1,size(wave,1)-length(spike_wave))]]);
            PsychPortAudio('Start', h_sound_device);
            true_presentation_period = toc(loop_time);
            loop_time =tic;
            if true_presentation_period > stimulus.presentation_period + 0.15 % allow 50 ms overun
                warning('ida:Presentation_period', ...
                    'Presentation period %f s larger than specified by',...
                    true_presentation_period - stimulus.presentation_period)
            end
        end
    end

    if strcmp(get(h_PB_stop,'text'),'Stopped') || stimGenError
        ida_exp_gui('CloseRequestFcn')
        break
    end
    first_presentation = 0;
end

%%% plot CAP results (quick fix)
% for q=1:length(h_subplots)
%     subplot(h_subplots(q))
%     ylim([min(neg_peak)*1.1, max(pos_peak)*1.1])
% end

if ishandle(CAP_threshold_fig)
    close(CAP_threshold_fig)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function freqs = calcualte_freqs(f_l,f_u,n_freqs)
% freqs are in matrix to play in staggered order
n_freqs = 4*ceil(n_freqs/4);
freqs = reshape(round(logspace(log10(f_u),log10(f_l),n_freqs)),[],4);

%%========================================================================
function tone_pips = generate_tone_pips(freqs,rep_rate,fs)
cos_ramp_cycles = 3; % cycles
total_duration = 1/rep_rate*fs;
tone_pips = zeros(size(freqs,1),size(freqs,2)*total_duration);
for q1 = 1:size(freqs,1)
    for q2 = 0:size(freqs,2)-1
        l_cos_ramp = round(cos_ramp_cycles/freqs(q1,q2+1)*fs);
        duration= 2*l_cos_ramp;
        tone_pips(q1,q2*total_duration+1:q2*total_duration+total_duration) = [...
            sin(linspace(0,freqs(q1,q2+1)/fs*duration*2*pi*(1-1/duration),duration))'... tone
            .* [cos(pi/2:(pi-pi/2)/(l_cos_ramp-1):pi).^2,cos(0:pi/2/(l_cos_ramp-1):pi/2).^2]';... enveleope
            zeros(total_duration-duration,1)]; % trailing silence
    end
end

%%========================================================================
function [CAP_ampl, stdev_noise] = extract_CAP_amplitude(wave_in)

%%========================================================================
function [levels,slope] = calcualte_next_levels(...
            levels,CAP_ampl,slopes,noise_level,thrshld)

%%========================================================================
function h_fig = init_CAP_threshold_fig(stdev_noise)
h_fig = uifigure;
h_fig.HandleVisibility = 'on';
h_fig.Tag = 'CAP_threshold_fig';
h_fig.Position = [0 60 1200 900];
h_fig.Resize = 0; 
h_fig.CloseRequestFcn  = 'CAP_tresholds CloseRequestFcn';

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
uibutton(h_fig, ...
	'text','save', ...
	'ButtonPushedFcn',@CB_save_button, ...
	'Tag','PB_save', ...
	'Position',[125 5 55 20]);
uieditfield(h_fig, ...
	'Value','counter', ...
    'FontSize',9, ...
	'Editable','off', ...
	'Tag','txt_counter', ...
	'Position',[0 30 120 15]);

% plot noise floor (stdev_noise)

%%========================================================================
function update_CAP_threshold_fig(h_fig,q,q2,levels,CAP_ampl,slopes,stdev_noise)
figure(h_fig)

%%========================================================================
function CB_save_button 
% When pressing save button in CAP_threshold_fig save data and close figure
global DATA_PATH

CAP_threshold_fig = gcbf;
h_ida_fig = findobj('Tag','ida_fig');
animalID=getappdata(h_ida_fig,'animalID');
savefig(CAP_threshold_fig,[DATA_PATH '\' datestr(header.datenum,'yyyymmddTHHMMSS'),'-',animalID,...
    '-CAPthrlds'])
close(CAP_threshold_fig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CloseRequestFcn
h_fig =findobj('Tag','CAP_threshold_fig');
delete(h_fig(end)),
PsychPortAudio('Close')
