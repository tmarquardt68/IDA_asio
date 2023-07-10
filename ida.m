function varargout = ida(varargin)
% Core functions of IDA_asio 
%
% Part of the IDA Toolbox 
% Copyright (C) Torsten Marquardt
% Terms of the GNU General Public License apply
% (www.http://www.fsf.org/licensing/licenses/gpl.html).
% v0.1 19-07-10 TM (Oldenburg)
% v1.1 18-10-12 TM (EI)
% v1.2 14-10-15 TM (EI)
% v1.3 07-03-22 TM (EI)
%
% TO DO:
% - edit boxes for trackNo and PosNo
% - $protocolReps: recover after restart from existing file names 
%
% - take care of rand state to repeat the exact experiment. E.g. if exp is
% interupted, read out the state and set it again on continuation, 
% in figure UI control or figure name.
% - buttons for sending file id (or adjust delay before stimulus starts
% with automatic start on/off tick box). Receipt of file ID will later start
% analysis software.
% - when restart, read current stimuli_parameter and experiment
% setting, compare with default, and if differing, store as different
% config file.
%
% 2023:
% - better variable names (replace in old result files!): q->rpt,q2->prmSet
% interval_order->presentation_order (also better store transposed)
% =========================================================================

if nargin == 0
    initialise;
elseif (nargout)
    [varargout{1:nargout}] = feval(varargin{1},varargin{2:end});
else
    feval(varargin{1}, varargin{2:end}); % FEVAL switchyard
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initialise
  
data_path = [fileparts(which(mfilename)) '\DATA\']; 
cd(fileparts(which(mfilename)))
addpath([fileparts(which(mfilename))])
addpath([fileparts(which(mfilename)) '\utilities\'])
addpath([fileparts(which(mfilename)) '\probe_check\'])
addpath([fileparts(which(mfilename)) '\CAP_gram\'])

% open and initialise IDA main figure
h =findobj('Tag','ida_fig');
if ~isempty(h)
    figure(h)
    return
end
h_fig = uifigure;
h_fig.HandleVisibility = 'on';
h_fig.Tag = 'ida_fig';
h_fig.Position = [0 0 255 90];
h_fig.Resize = 0; 
h_fig.CloseRequestFcn  = 'ida CloseRequestFcn';

files = my_dir(data_path,'mat');
[n,~] = size(files);
last_file = deblank(files(n,:)); %% ignore animalID 'test'!
tmp = strsplit(last_file,'-');
while (strcmp(tmp{2},'test')||strcmp(tmp{2},'tst')) && n>0
    n=n-1;
    if n>0
        last_file = deblank(files(n,:));
    end
    tmp = strsplit(last_file,'-');
end

if n==0
    animal_notes='';
    userID=''; 
    animalID='';
    track_no=0;
    pos_no=0;
    last_file='';
else
    load([data_path, last_file],'results')
    header=results.header;
    userID=header.experimenter;
    animalID=header.animalID;
    track_no=header.track_no;
    pos_no=header.pos_no;
    if isfield(header,'animal_notes')
        animal_notes=header.animal_notes;
    else
        animal_notes='';
    end
end

add_ida_paths(userID);
setappdata(h_fig,'userID',userID);
setappdata(h_fig,'animalID',animalID);
setappdata(h_fig,'track_no',track_no);
setappdata(h_fig,'pos_no',pos_no);
setappdata(h_fig,'animal_notes',animal_notes)
setappdata(h_fig,'last_file',last_file);
setappdata(h_fig,'data_path',data_path);
h_fig.Name = ...
    sprintf('IDA %s,  %s, %d, %d',getappdata(h_fig,'userID'),...
    getappdata(h_fig,'animalID'),getappdata(h_fig,'track_no'),...
    getappdata(h_fig,'pos_no'));

uibutton(h_fig, ...
	'text','reset sound I/O', ...
	'ButtonPushedFcn','PsychPortAudio(''Close'');', ...
	'BusyAction','cancel', ...
	'Position',[5 5 100 20]);
uibutton(h_fig, ...
	'Text','cont. exp', ...
	'ButtonPushedFcn',@CB_cont_experiment, ...
	'BusyAction','cancel', ...
	'Position',[5 25 100 20]);
uibutton(h_fig, ...
	'Text','animal notes', ...
	'ButtonPushedFcn',@CB_animal_notes, ...
	'BusyAction','cancel', ...
	'Position',[5 45 100 20]);
uibutton(h_fig, ...
	'Text','browse data', ...
	'ButtonPushedFcn', @CB_browse, ...
	'BusyAction','cancel', ...
	'Position',[5 65 100 20]);
uibutton(h_fig, ...
	'Text','probe check', ...
	'ButtonPushedFcn',@CB_probe_check, ...
	'BusyAction','cancel', ...
	'Position',[110 5 100 20]);
uibutton(h_fig, ...
	'Text','CAP thresholds', ...
	'ButtonPushedFcn',@CB_CAP_thresholds, ...
	'BusyAction','cancel', ...
	'Position',[110 25 100 20]);
uibutton(h_fig, ...
	'Text','change subject', ...
	'ButtonPushedFcn',@CB_change_subject, ...
	'BusyAction','cancel', ...
	'Position',[110 45 100 20]);
uibutton(h_fig, ...
    'Text','change user',...
	'ButtonPushedFcn',@CB_change_user, ...
	'BusyAction','cancel', ...
	'Position',[110 65 100 20]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA AQUISITION functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function run_experiment(config_filename,stimulus,...
                        parameter_table,exp_description)
% stimulus presentation routine of ida to start from the command line

h_ida_fig = findobj('Tag','ida_fig');
data_path = getappdata(h_ida_fig,'data_path');
userID = getappdata(h_ida_fig,'userID');
animal_notes = getappdata(h_ida_fig,'animal_notes');
animalID=getappdata(h_ida_fig,'animalID');
if isnumeric(animalID), animalID = num2str(animalID); end
track_no=getappdata(h_ida_fig,'track_no');
pos_no=getappdata(h_ida_fig,'pos_no');

% create header field
header.datenum = now;
header.experimenter = userID;
header.animalID = animalID;
header.animal_notes = animal_notes;
header.track_no = track_no;
header.pos_no = pos_no;
header.comments{1}='';
header.config_filename = config_filename;
header.exp_description = exp_description; % Short string descibing experiment/configuration (typically < 50 characters)
header.title = [datestr(header.datenum,'yyyymmddTHHMMSS'),'-',animalID,...
    '-',config_filename];

results.header = header;
results.stimulus = stimulus;
results.stimulus.original_parameter_table=parameter_table;

try
    load('C:\Users\Admin\Documents\MATLAB\IDA_asio\DP_TT\_last_probe_check','CURRENT_EAR'),
    results.H_mic = CURRENT_EAR(:,4);
catch
    h = errordlg('Cannot find DP_TT\_last_probe_check.mat!');
    uiwait(h),
    return,
end

q = 1; q2 = 1;
save([data_path, header.title], 'results', 'q', 'q2')
continue_experiment(header.title)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function continue_experiment(exp_filenname)
%% Start data aquisition

h_ida_fig = findobj('Tag','ida_fig');
data_path = getappdata(h_ida_fig,'data_path');

if ~exist('exp_filenname','var')
    exp_filenname = uigetfile(data_path,'Continue experiment:');
    if ~exp_filenname, return, end
end

load([data_path, exp_filenname],'results','q','q2')
header = results.header;
stimulus = results.stimulus;
%? H_mic = results.H_mic; % for monitoring function
fs = results.stimulus.sample_rate;
parameter_table = results.stimulus.original_parameter_table;
cd([fileparts(which(mfilename)) '\USERS\' header.experimenter '\' header.config_filename])
eval(['monitor_settings = monitor_init_' header.config_filename '(results);']);

% open and customise exp_fig
h_exp_fig = ida_exp_gui('initialise');
setappdata(h_ida_fig, 'h_exp_fig', h_exp_fig)
set(h_exp_fig, 'Name', [header.config_filename,', ' header.exp_description])
ida_position_fig = get(h_ida_fig, 'Position');
position_fig = get(h_exp_fig, 'Position');
set(h_exp_fig, 'Position',...
    [ida_position_fig(1)+ida_position_fig(3), ida_position_fig(2)...
    position_fig(3) position_fig(4)]);
h_PB_pause = findobj(h_exp_fig,'Tag','PB_pause');
h_PB_stop = findobj(h_exp_fig,'Tag','PB_stop');

% initialise stimuli_GUI
position_fig = get(h_exp_fig, 'Position');
h_stimuli_gui = stimuli_gui(h_exp_fig, ...
    [position_fig(1)+position_fig(3), position_fig(2)], ...
    [header.config_filename, '.mat']);
setappdata(h_exp_fig, 'h_stimuli_gui', h_stimuli_gui);
stimuli_gui('Enable_uicontrols',h_stimuli_gui,'off');

if isempty(stimulus.presentation_period) % loop stimulus as fast as possible
    stimulus.presentation_period = ...
        length(get_stimulus_waveforms(h_stimuli_gui, stimulus.parameters(2:end), ...
            parameter_table(1,2:size(parameter_table,2))))/fs;
end
setappdata(h_exp_fig,'presentation_period',stimulus.presentation_period)

spike_wave = create_trigger_wave(fs);
% test length of intrval ID wave < presentation_period (presentation_start_marker length: 5*13 ms)
if  length(spike_wave)*(5*13 + 4*ceil(log2(size(parameter_table,1)))+40)/fs > ...
        stimulus.presentation_period
    error('intervalID_wave longer than stimulus period')
end

% % play fileID
% fileID_wave = create_fileID_wave(fs,header.datenum*1e9);
% sound_io('play',h_sound_device, [zeros(fs/10+length(fileID_wave),2) ...
%     [zeros(fs/10,1); fileID_wave'] ...
%     zeros(fs/10+length(fileID_wave),1)]);
% pause(.5) % time to start spike analysis software

% ==== Initialise sound output device
h_sound_device = sound_io('open_device',fs, 4, 2);
PsychPortAudio('GetAudioData', h_sound_device, 2*stimulus.presentation_period);
PsychPortAudio('FillBuffer',h_sound_device,zeros(4,stimulus.presentation_period*fs));
PsychPortAudio('Start', h_sound_device);

% play n_rpts times a sequences of the intervals
loop_time=tic;
no_presentations=(stimulus.n_rpts)*(size(parameter_table,1));
first_presentation = 1; prmSet = 0; stimGenError = 0;
for n = (q-1)*size(parameter_table,1)+q2:no_presentations+1
    if n == 2, initSeriesFig = 1; else, initSeriesFig=0; end
    prmSetToSave = prmSet;
    q = ceil(n/size(parameter_table,1));
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
        prmSet = stimulus.interval_order(q,q2);
        [q q2 prmSet]; % print in command window
        try
            wave = get_stimulus_waveforms(h_stimuli_gui, stimulus.parameters(2:end), ...
                parameter_table(prmSet,2:size(parameter_table,2)));
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % for electrostatic speakers only
%             wave(:,1:2) = ((200*wave(:,1:2)+212.7).^0.5-sqrt(212.7))/sqrt(200);
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        catch
            stimGenError = 1;
        end
        if abs(wave)>1
            warndlg('Max of summed wave > 1!')
            stimGenError = 1;
        end
        % generate interval ID channel (#1) waveform
        intervalID_wave = create_intervalID_wave(fs,...
            parameter_table(prmSet,1),q2);
        set(findobj(h_exp_fig,'Tag','txt_counter'),'Value',['Rpt:' num2str(q), '/' num2str(stimulus.n_rpts)...
            ', PrmSet:' num2str(q2), '(' num2str(prmSet) ')' ...
            ' of' num2str(size(parameter_table,1))])
    end
%     if n > 1
%         t = toc(processing_time);
%         pause(stimulus.presentation_period-t +0.1) % allow 100 ms ISI
%     end
%     processing_time = tic;

    % present wave form
    if ishandle(h_exp_fig) % if exp figure still open
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
                set(h_exp_fig,'Color','r')
            end
        end
    end
    if ~first_presentation % save data and call monitoring display function
        eval(['monitor_' header.config_filename '(initSeriesFig, prmSetToSave, results, wave_in, monitor_settings, 1);'])
        eval(['raw_' num2str(n-1) '.wave = wave_in;'])
        eval(['raw_' num2str(n-1) '.timeStamp = now;']) % datetime(datevec(now))
        save([data_path, header.title], ['raw_' num2str(n-1)],'-append')
    end

    if strcmp(get(h_PB_stop,'text'),'Stopped') || stimGenError
        ida_exp_gui('CloseRequestFcn')
        break
    end
    first_presentation = 0;
end

if ishandle(h_exp_fig)
    close(h_exp_fig)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions for Data acquisition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function wave = get_stimulus_waveforms(h_stimuli_gui, dynamic_parameter, ...
    parameters)
err = zeros(1,size(parameters,2));
for q = 1:size(parameters,2)
    err(q) = stimuli_gui('Change_parameter', h_stimuli_gui, ...
        dynamic_parameter(q).stimulus_no, ...
        dynamic_parameter(q).name , ...
        parameters(q));
end

if sum(err)
    error('Paramter change failed.')
end
wave = stimuli_gui('Get_waveform','B',h_stimuli_gui);
if isempty(wave)
    error('stimulus wave is empty. Is stimulus selected for interval ''B''?')
end

%% =========================================================================
function fileID_wave = create_fileID_wave(fs,fileID)
% 6ms-inteval marker sequence of 5 spikes followed by 10 silence intervals, then
% 6ms-inteval binary file ID
binary_fileID = dec2bin(fileID);

% waveform of simulated spike
spike_wave = create_trigger_wave(fs);
k = length(spike_wave);

% file start marker
fileID_wave = [repmat([spike_wave, zeros(1,5*k)],1,5), zeros(1,60*k)];

for q = length(binary_fileID):-1:1 % little endian
    if strcmp(binary_fileID(q),'1')
        fileID_wave = [fileID_wave spike_wave];
    else,
        fileID_wave = [fileID_wave zeros(1,k)];
    end
    fileID_wave  = [fileID_wave zeros(1,5*k)];
end
fileID_wave  = [fileID_wave zeros(1,fs/2)]; % half second silence

%% =========================================================================
function intervalID_wave = create_intervalID_wave(fs,intervalID,index)
% max number of intervals is 2^16!!! shortest repetition rate is 150 ms!!!
% interval start marker with 5ms raster:
if index == 1
    presentation_start_marker_binary = '1111111110000'; % start of new repetition
else
    presentation_start_marker_binary = '1110001110000';
end
% intervalID with 4ms raster
binary_intervalID = dec2bin(intervalID);

% waveform of simulated spike (Must be 1 ms duration!)
spike_wave = create_trigger_wave(fs);
k = length(spike_wave);

% create interval-start marker wave
intervalID_wave = zeros(ceil(130/fs),1);
n = 1;
for q = 1:length(presentation_start_marker_binary)
    if strcmp(presentation_start_marker_binary(q),'1')
        intervalID_wave(n:n+k-1) = spike_wave;
    else,
        intervalID_wave(n:n+k-1) = zeros(1,k);
    end
    n = n+k;
    intervalID_wave(n:n+4*k-1)  = zeros(1,4*k);
    n = n+4*k;
end
% create intervalID wave (starts after (4*5 + 4) ms silence)
n = n-k;
for q = length(binary_intervalID):-1:1 % little endian
    if strcmp(binary_intervalID(q),'1')
        intervalID_wave(n:n+k-1) = spike_wave;
    else,
        intervalID_wave(n:n+k-1) = zeros(1,k);
    end
    n = n+k;
    intervalID_wave(n:n+3*k-1)  = zeros(1,3*k);
    n = n+3*k;
end
% test: plot(intervalID_wave,'r')

%% =========================================================================
function fileID = get_fileID(event_ts)
% 6ms-inteval marker sequence of 5 spikes followed by 10 silence intervals, then
% 6ms-inteval binary file ID

bit_no = round((event_ts(event_ts < event_ts(1)+ 0.5))./0.006);
binary_str = repmat('00000',1,13);
binary_str(bit_no - bit_no(1)+1) = '1';
fileID = bin2dec(binary_str(end:-1:16))/1e9;

% datestr(fileID,'yyyymmddTHHMMSS')
% edges =  min(event_ts)-.1:0.001:min(event_ts)+.6; bar(edges,histc(event_ts, edges),'histc')

%% =========================================================================
function intvl_starts = get_intvl_starts(event_ts)
% returns matrix with 3 columns
% 1st column: interval start time
% 2nd column: 1 - within sequence; 2 - new sequence (repetition) starts
% 3rd column: interval number

% Number of rows equals total number of intervals in plx file
intvl_start_template = get_intvl_start_template;
seqnc_start_template = get_seqnc_start_template;

is_marker = zeros(length(event_ts),1);
first_sequence_started = 0;
n = 1;
while n+length(intvl_start_template)-1 < length(event_ts)
    %  if interval start
    if (sum(abs(event_ts(n:n+length(intvl_start_template)-1)- ...
            intvl_start_template - event_ts(n))) < 0.5e-3) ...
            & first_sequence_started
        is_marker(n) = 1;
        n = n + length(intvl_start_template);
        % if sequence start
    elseif n+length(seqnc_start_template)-1 < length(event_ts) && ...
            abs(sum(event_ts(n:n+length(seqnc_start_template)-1)- ...
            seqnc_start_template - event_ts(n))) < 0.5e-3
        is_marker(n) = 2;
        first_sequence_started = 1;
        n = n + length(seqnc_start_template);
    else
        n=n+1;
    end
end
intvl_starts = [event_ts(is_marker>0) is_marker(is_marker>0)];
intvl_starts = [intvl_starts, zeros(size(intvl_starts,1),1)];

% intvl_no = zeros(length(intvl_starts),1);
for q = 1:size(intvl_starts,1)
    ts1 = intvl_starts(q,1);
    if q < size(intvl_starts,1)
        ts2 = intvl_starts(q+1,1);
    else
        ts2 = ts1 + 1;
    end
    bit_no = round((event_ts(event_ts > ts1+intvl_start_template(end) ...
        +1e-3 & event_ts < ts2) - (ts1+intvl_start_template(end)))./0.004);
    bit_no = bit_no(bit_no <= 16);
    binary_str = '0000000000000000';
    binary_str(bit_no-5) = '1';
    intvl_starts(q,3) = bin2dec(binary_str(end:-1:1));
end

%% =========================================================================
function intvl_start_template = get_intvl_start_template
intvl_start_template = 0:0.005:0.04;
intvl_start_template = intvl_start_template([1:3,7:9])';

%% =========================================================================
function seqnc_start_template = get_seqnc_start_template
seqnc_start_template = (0:0.005:0.04)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CB_change_user(callingObj,~)
h_fig = callingObj.Parent;
oldUserID=getappdata(h_fig,'userID');
userList=dir([fileparts(which(mfilename)) '\USERS\']);
for q=3:length(userList), users{q-2}=userList(q).name;end
newUserID=listdlg('ListString',users,'SelectionMode','single',...
    'PromptString','Select user:',...
    'Name','User selection',...
    'ListSize',[120 100]);
if ~isempty(newUserID)
    newUserID=users{newUserID};
    if ~isempty(oldUserID)
        remove_ida_paths(oldUserID);
    end
    add_ida_paths(newUserID);
    setappdata(h_fig,'userID',newUserID);
    h_fig.Name = ...
        sprintf('IDA %s,  %s, %d, %d',getappdata(h_fig,'userID'),...
        getappdata(h_fig,'animalID'),getappdata(h_fig,'track_no'),...
        getappdata(h_fig,'pos_no'));
end
%% =========================================================================

function CB_change_subject(callingObj,~)
h_fig = callingObj.Parent;
animalID = inputdlg('Give new subject iD:','Change subject ID',1,{getappdata(h_fig,'animalID')});
if isempty(animalID)
    return,
end
setappdata(h_fig,'animalID',animalID{1});
setappdata(h_fig,'track_no',0);
setappdata(h_fig,'pos_no',0);
h_fig.Name = ...
    sprintf('IDA %s,  %s, %d, %d',getappdata(h_fig,'userID'),...
    getappdata(h_fig,'animalID'),getappdata(h_fig,'track_no'),...
    getappdata(h_fig,'pos_no'));

%% =========================================================================
function CB_CAP_thresholds(callingObj,~)
CAP_thresholds

%% =========================================================================
function CB_probe_check(callingObj,~)
dp_tt_probe_check ini

% %% =========================================================================
% function CB_incr_pos(callingObj,~)
% h_fig = callingObj.Parent;
% setappdata(h_fig,'pos_no',getappdata(h_fig,'pos_no')+1);
% h_fig.Name = ...
%     sprintf('%s,  %s, %d, %d',getappdata(h_fig,'userID'),...
%     getappdata(h_fig,'animalID'),getappdata(h_fig,'track_no'),...
%     getappdata(h_fig,'pos_no'));
% h_fig = callingObj.Parent;
% animal_notes=getappdata(h_fig,'animal_notes');
% animal_notes=inputdlg({'Add note'},'Note',[40 80],{animal_notes});
% if ~isempty(animal_notes)
%     setappdata(h_fig,'animal_notes',animal_notes);
% end

%% ========================================================================
function CB_cont_experiment(callingObj,~)
h_fig = callingObj.Parent;
continue_experiment(getappdata(h_fig,'last_file'))

%% =========================================================================
function CB_animal_notes(callingObj,~)
h_fig = callingObj.Parent;
animal_notes=getappdata(h_fig,'animal_notes');
if isempty(animal_notes{1})
    animal_notes=inputdlg({'Add note'},'Notes about the subject:',[40 80]);
else
    animal_notes=inputdlg({'Add note'},'Notes about the subject:',[40 80],{animal_notes{1}});
end
if ~isempty(animal_notes)
    setappdata(h_fig,'animal_notes',animal_notes);
end

%% =========================================================================
function CB_browse(callingObj,~)
h_fig = callingObj.Parent;
data_path = getappdata(h_fig,'data_path');
animalID = getappdata(h_fig,'animalID');

h_browser_fig = ida_data_browser('initialise',data_path,animalID);
setappdata(h_fig,'h_browser_fig',h_browser_fig)


results = []; % load data

% Establish the name of the file to save the processed data in.
[~,tempname,ext] = fileparts(results.header.title);
index=1;
done=false;
while ~done
    resultsFile = [data_path tempname '_' num2str(index) ext];
    if exist(resultsFile,'file')
        index=index+1;
    else
        done=true;
    end
end
save(resultsFile, 'results')
eval(['cd USERS\' userID '\' results.header.config_filename])
eval(['plot_' results.header.config_filename '(results)'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions for init and callbacks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function add_ida_paths(userID)
basepath=[fileparts(which(mfilename)) '\USERS\' userID];
dirs=dir(basepath);
for n=3:length(dirs)
    if ~strncmp(dirs(n).name,'.',1)&&dirs(n).isdir
        addpath([basepath '\' dirs(n).name],'-begin');
    end
end

%% =========================================================================
function remove_ida_paths(userID)
basepath=[fileparts(which(mfilename)) '\USERS\' userID];
dirs=dir(basepath);
for n=1:length(dirs)
    if ~strncmp(dirs(n).name,'.',1)&&dirs(n).isdir
        rmpath([basepath '\' dirs(n).name]);
    end
end

%% =========================================================================
function CloseRequestFcn
h_ida_fig = findobj('Tag','ida_fig');
h_exp_fig = getappdata(h_ida_fig, 'h_exp_fig');
h_browser_fig = getappdata(h_ida_fig, 'h_browser_fig');
remove_ida_paths(getappdata(h_ida_fig,'userID'));
rmpath([fileparts(which(mfilename)) '\utilities\'])

closereq,
if ishandle(h_exp_fig)
    close(h_exp_fig)
end
if ishandle(h_browser_fig)
    close(h_browser_fig)
end