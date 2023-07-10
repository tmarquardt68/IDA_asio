function [vfPulseTimes vfSpikeTimes] = get_raw_data(filename,channels,units)
% system specific routine to read physiological raw data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: SJJones
% Date: 2013/01/15
% Version: 1.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extracts spike and event timestamps from a TDT datatank.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Change Log:
% 2013/01/15 - SJJ Created File
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% STEP ONE - Check and enter user input
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vsBlock = 'IAC_3255';
% vsTankName = '20121115';
% vsTankPath = 'C:\TDT\OpenEx\MyProjects\SortExampleBckUp\DataTanks\';

% TO DO:
% - load specifiv channels/units

% assign user input
[vsTankName vsBlock] = fileparts(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% STEP TWO - Activate server
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

auTT = actxcontrol('TTank.X'); % Open Active X
auTT.ConnectServer('Local','Me'); % Contact server
auTT.OpenTank(vsTankName,'R'); % Target tank
auTT.SelectBlock(vsBlock); % Target block

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% STEP THREE - Call server and extract data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tag name needs to be accurate and max records number enough - no warning

eventCodes=auTT.GetEventCodes(0);
eventNames=cell(1,length(eventCodes));
for n=1:length(eventNames)
    eventNames{n}=auTT.CodeToString(eventCodes(n));
end

if any(strncmp('EA__',eventNames,4))
    anNumRecs = auTT.ReadEventsV(100000,'EA__',0,0,0,0,'ALL'); % number of records
elseif any(strncmp('Snip',eventNames,4))
    anNumRecs = auTT.ReadEventsV(10e6,'Snip',0,0,0,0,'ALL'); % number of records
else
    error('get_raw_data: no data store for spikes');
end
mfTrodeInfo = auTT.ParseEvInfoV(0,anNumRecs,0); % extract the events

if any(strncmp('Valu',eventNames,4))
    anNumRecs = auTT.ReadEventsV(1e6,'Valu',0,0,0,0,'ALL'); % number of records
else
    error('get_raw_data: No data store for ida tag');
end
mfEventInfo = auTT.ParseEvInfoV(0,anNumRecs,0); % extract the events

% vsRecDate = auTT.FancyTime(auTT.CurBlockStartTime,'D-O-Y'); % record date

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% STEP FOUR - Shut down tank/server
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

auTT.CloseTank; % close tank
auTT.ReleaseServer; % disconnect server


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% STEP SIX - Extract, order and adjust pulse & spike times
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pulses
if length(mfEventInfo)<=1
    error('get_raw_data: no ida pulses in data tank')
end
vfPulseTimes = mfEventInfo(6,:); % extract pulses

% spikes from pulses
if length(mfTrodeInfo)<=1
    vfSpikeTimes=[];
else
    vfSpikeTimes = mfTrodeInfo(6,:); % extract spikes
end