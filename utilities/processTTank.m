function results = processTTank(filename,varargin)
% processTTank - extract data from an OpenEx tank and put it into the
% Oldenburg data format.  The function makes strong assumptions about the
% names of the fields present -- watch for warnings!
%
% RESULTS=processTTank(FILENAME) will process the block listed in
% FILENAME. Must include the full path. RESULTS will be the Oldenburg data
% structure.
%
% RESULTS=processTTank(FILENAME,'OPTNAME',OPTVAL,...) will set the option
% with name OPTNAME to the value OPTVAL.  If you want graphical selection
% of files, set FILENAME to [].  You can specify as many option pairs as
% you want. Valid options are:
%
% 'doAnalog' - whether or not to process analog data (true(default)/false)
% 'doSave' - whether or not to save the output (true/false(default))
% 'numChannels' - exists because Bjorn is lazy and hasn't written channel
%   deciphering yet.
% 'ecogFs' - sampling rate of the ECOG data saved in the .mat file
%   (in Hz, default 1000).
%
% there are a couple more options but you probably shouldn't override them.

% v1.0 Written by GBC, 14 May 2013

%% Variable initialisation

% Specify some defaults.  Note that the default for doAnalog is true and is
% set above in the filename handling code above.
% doAnalog = true; % Process analogue data.
maxNumEvents=10e6; % TDT tank format does not allow you to say "all events"
numChannels=16; % number of analog channels -- shouldn't be hardcoded
% but have not figured out how to extract the number
% yet
ecogFs=1000;    % Sampliing rate in Hz.  Approximate.
butterN=4;      % Order of the lowpass filter for the ECoG data.
ecogFilterLeeway = 0.1; % Percentage space away from Nyquist to allow.
readWindow=10; % in seconds; for batch processing of analogue data.
samplesToRead=10000;  %  for batch processing of analogue data.
doAnalog=false;
doLFP=true;
% doSave=false;    % for some applications we prefer not to save the data.
abDoPlotting=false; % just a stupid workaround.
dataPath = '';


% Below code will set options appropriately.
warnopts(assignopts(who,varargin{:}));

% Some required globals for the analog file output structure.
if doAnalog
    lsfglobals;
end

% OpenEx returns event data with rows corresponding to types of
% information.  For the information we care about, we'll create a name to
% use rather than just the row number to make this readable later on.
% Check OpenDeveloper documentation.
epochTimeIdx=2;
chanIdx=4;
scIdx=5;
timeIdx=6;

%% TDT initialisation
% Create the ActiveX control in a new figure and then hide it because it's
% ugly.  Then initialze the TDT server.
Hf=figure;
ox = actxcontrol('TTank.X','parent',Hf); % Open Active X
set(Hf,'Visible','off');
drawnow();
ox.ConnectServer('Local','Me'); % Contact server

[tankName,blockName]=fileparts(filename);

%Load the tank and block, and index the epochs.
ox.OpenTank(tankName,'R'); % Target tank
ox.SelectBlock(blockName); % Target block
ox.CreateEpocIndexing();

% Find the names of the variables stored in the tank.
eventCodes=ox.GetEventCodes(0);
eventNames=cell(1,length(eventCodes));
for n=1:length(eventNames)
    eventNames{n}=ox.CodeToString(eventCodes(n));
end

%% Trial boundaries
% The 'Epoc' event is a single timestamp at the beginning of each trial.
% It must exist, or there's no way to process the file.  Check, then load
% it.
if ~any(strncmp('Epoc',eventNames,4))
    error('ida:TTankHandling:epoch',...
        'No Epoc field present in data tank!');
end
fprintf(1,'\t Parsing trial boundaries ...\n');
% this loop deals with the fact that there's no way to check the number of
% events before we start.
done=false;
while ~done
    trialTimes=ox.GetEpocsV('Epoc',0,0,maxNumEvents);
    if isnan(trialTimes)
        error('ida:TTankHandling:sizeExceeded',...
            'Too many events to handle!');
    elseif size(trialTimes,2)==maxNumEvents
        maxNumEvents=maxNumEvents*2;
    else
        done=true;
    end
end

% Get the timestamps, work out the number of trials.  In general the
% endpoint of a trial will be the beginning of the next trial, so add a
% stopping point for the very last trial.
trialTimes=trialTimes(epochTimeIdx,:);
numTrials=length(trialTimes);
trialTimes(end+1)=trialTimes(end)+max(diff(trialTimes));


% Now we need to parse the ida identifiers.  Again, if these are not in
% the file we will have major problems.
if ~any(strncmp('Epoc',eventNames,4))
    error('ida:TTankHandling:epoch',...
        'No Epoc field present in data tank!');
end
fprintf(1,'\t Parsing ida identifiers...\n');
idaID=nan(numTrials,1);
idaTimestamps=cell(numTrials,1);
done=false;
while ~done
    % In theory we should be able to treat this as Epoch data, but on
    % large files there seems to be a bug where the later timestamps
    % are all returned as 0.  Treating it as Events seems to prevent
    % this.
    timeStamps=ox.ReadEventsV(10e6,'Valu',0,0,0,0,'ALL');
    if isnan(timeStamps)
        error('ida:TTankHandling:sizeExceeded',...
            'Too many events to handle!');
    elseif size(timeStamps,2)==maxNumEvents
        maxNumEvents=maxNumEvents*2;
    else
        done=true;
    end
end
timeStamps=ox.ParseEvInfoV(0,timeStamps,0);
timeStamps=timeStamps(timeIdx,:);

% The following code is stolen/modified from
% ida.m/get_fileID.m, v 1.1, 8 May 2013.
% The first half second of the ida timestamps give the date identifier
% used in the filename as a unique ID.
bit_no = round((timeStamps(timeStamps < timeStamps(1)+ 0.5))./0.006);
binary_str = repmat('00000',1,13);
binary_str(bit_no - bit_no(1)+1) = '1';
fileID = datestr(bin2dec(binary_str(end:-1:16))/1e9,'HHMMSS');

% Very similar to extracting spikes, but simpler because there are
% no sort codes and only one channel.
for m=1:numTrials
    idaTimestamps{m}=timeStamps(timeStamps>=trialTimes(m) & ...
        timeStamps<trialTimes(m+1)) - (trialTimes(m));
    % The below code is stolen/modified from
    % ida.m/get_fileID.m, v1.1, 8 May 2013.
    % The selected range of ida timestamps encodes the interval number
    % presented in the trial.
    bit_no = round((timeStamps(timeStamps>trialTimes(m) + .041 & ...
        timeStamps<trialTimes(m+1)) - (trialTimes(m) + .04)) ...
        ./0.004);
    bit_no = bit_no(bit_no <= 16);
    binary_str = '0000000000000000';
    binary_str(bit_no-5) = '1';
    idaID(m)=bin2dec(binary_str(end:-1:1));
end

%% Output file initialisation
% With the ida information processed, we can load the stimuli file. Do
% error handling if we can't.
if isnan(fileID)
    error('ida:TTankHandling:noMATfile', ...
        'Unable to determine .MAT file for %s!',filename);
end
stimuliFile=dir([dataPath filesep '*' fileID '*']);
if isempty(stimuliFile)
    error('ida:TTankHandling:noMATfile', ...
        'Unable to find .MAT file for %s!',filename);
else
    stimuliFile=[dataPath filesep stimuliFile(1).name];
    load(stimuliFile,'results')
    if ~exist('results','var')
        error('ida:TTankHandling:improperMATfile', ...
            '%s does not contain results variable!',stimuliFile);
    end
end

% save mandatory fields in the header
results.header.date=datestr(results.header.datenum,'DD-mmmm-YYYY');
results.header.datasource=filename;

% Save the extracted information about trial boundaries.
results.stimulus.trialMarkers.trialTimes=trialTimes(1:end-1);
results.stimulus.trialMarkers.idaID=idaID;
results.stimulus.trialMarkers.idaTimestamps=idaTimestamps;

% Save the stimulus order -- couldn't do this until we deciphered the
% ida identifiers.
if isfield(results.stimulus,'ordering_scheme')
    results.stimulus.order = ...
        results.stimulus.ordering_scheme.block.original_table_order(idaID,:);
end

% Set up a template for the output.
physiological.measures=struct('name',{},'unit',{},'values',{},'description',{});
physiological.raw=filename;
physiological.data=[];
physiological.datenum_extraction=now;

% Because of the nature of the Oldenburg data format, we don't have a way
% to work out how many data channels we will end up with.
cnt=0;

%% Now process all the event types present in the tank.
for n=1:length(eventNames)
    switch eventNames{n}
        case {'Snip' 'EA__'}
            %% Spike times, which in the past have had a couple different names.
            done=false;
            while ~done
                spikes=ox.ReadEventsV(10e6,eventNames{n},0,0,0,0,'ALL');
                if isnan(spikes)
                    error('ida:TTankHandling:sizeExceeded',...
                        'Too many events to handle!');
                elseif size(spikes,2)==maxNumEvents
                    maxNumEvents=maxNumEvents*2;
                else
                    done=true;
                end
            end
            % The Oldenburg data format has a separate data entry for each
            % channel and each spike category.
            spikes=ox.ParseEvInfoV(0,spikes,0);
            if numel(spikes)>1
                fprintf(1,'\t Parsing spikes...\n');
                channelsPresent=unique(spikes(chanIdx,:));
                sortcodesPresent=unique(spikes(scIdx,:));
                for c=1:length(channelsPresent)
                    for sc=1:length(sortcodesPresent)
                        % Work out all spikes that match the current channel
                        % and sort code.  Might not be any, so check that.
                        events=spikes(timeIdx,spikes(chanIdx,:)==channelsPresent(c) & ...
                            spikes(scIdx,:)==sortcodesPresent(sc));
                        if ~isempty(events)
                            cnt=cnt+1;
                            physiological.measures(cnt).name='spike times';
                            physiological.measures(cnt).unit='seconds';
                            physiological.measures(cnt).values=[];
                            physiological.measures(cnt).raw=filename;
                            physiological.measures(cnt).description=...
                                sprintf('channel %d, unit %d',...
                                channelsPresent(c),sortcodesPresent(sc));
                            for m=1:numTrials
                                physiological.data{m,cnt}=events( ...
                                    events >= trialTimes(m) & ...
                                    events < trialTimes(m+1));
                            end
                        end
                    end
                end
            end
        case 'RawT'
%             %% Analogue data, field potentials.
%             if doAnalog & doLFP
%                 fprintf(1,'\t Parsing field potentials...\n');
%                 % Extract the sampling rate of the wavenform.
%                 Fs=ox.ReadEventsV(1,'RawT',0,0,0,0,'ALL');
%                 Fs=ox.ParseEvInfoV(0,Fs,0);
%                 Fs=unique(Fs(9,:));
%                 Fs=floor(Fs);
%                 ox.SetGlobals(sprintf('WaveSF=%d',Fs));
%                 for c=1:numChannels
%                     % We need to save a descriptor in the .MAT file so
%                     % that we know that the field potential data was
%                     % present and where it is stored.  However, we
%                     % don't actually save any in the .MAT file.
%                     cnt=cnt+1;
%                     physiological.measures(cnt).name='field potential';
%                     physiological.measures(cnt).unit='';
%                     physiological.measures(cnt).values=[];
%                     physiological.measures(cnt).raw=lfpFile;
%                     physiological.measures(cnt).description= ...
%                         sprintf('lfp channel %d',c);
%                     physiological.measures(cnt).samplingRate=Fs;
%                     for m=1:numTrials
%                         physiological.data{m,cnt}=[];
%                     end
%                 end
%                 
%                 % Now create the separate file containing the field
%                 % potential data.
%                 if doSave
%                     ox.SetGlobalV('Channel',0);
%                     thisWin=0;
%                     done=false;
%                     sz=0;
%                     maxval=0;
%                     
%                     % First pass: We'll read all the data out of the tank
%                     % and put it in a temporary file.  We need to calculate
%                     % the maximum absolute deflection for the rescaling
%                     % later.
%                     fid=fopen(tmpFile,'Wb');
%                     while ~done
%                         ox.SetGlobalV('T1',thisWin*readWindow);
%                         ox.SetGlobalV('T2',(thisWin+1)*readWindow);
%                         waves=double(ox.ReadWavesV('RawT'))';
%                         if isempty(waves)||(length(waves)==1&&isnan(waves))
%                             done=true;
%                         else
%                             maxval=max([maxval max(abs(waves))]);
%                             if thisWin==0
%                                 % OpenEx seems to have the first couple of values be an
%                                 % artefact.  Solve this by throwing out the first
%                                 % millisecond for the purpose of calculating the max.
%                                 maxval=max(max(abs(waves(round(1e-3*Fs)))));
%                                 fwrite(fid,waves,'float32');
%                             else
%                                 maxval = max([maxval max(max(abs(waves)))]);
%                                 fwrite(fid,waves,'float32');
%                             end
%                             sz=sz+numel(waves);
%                         end
%                         thisWin=thisWin+1;
%                     end
%                     fclose(fid);
%                     
%                     % We want to store as integer values to minimise file
%                     % size.  So we'll scale the data to go from +/- 32767
%                     % (maximum value of a 16 bit integer).
%                     scale=maxval/32767;
%                     fin=fopen(tmpFile,'r');
%                     fout=fopen(lfpFile,'Wb','ieee-be');
%                     % Create the header for the .lsf file
%                     lsfwriteheader(fout,'numchannels',numChannels,'srate',Fs,...
%                         'scale',1,'encoding',LSF_ENC_LINEAR_16,...
%                         'info',['ecog from ' filename],...
%                         'datasize',sz*2);
%                     
%                     % Now read in all the data from the temp file and add
%                     % it to the .lsf file.  If we get to a point where
%                     % we've read less than what we asked for, we've reached
%                     % the end of the file.
%                     count=samplesToRead;
%                     while count==samplesToRead
%                         [waves,count]=fread(fin,samplesToRead,'float32');
%                         waves=round(waves./scale);
%                         fwrite(fout,waves,'int16');
%                     end
%                     fclose(fin);
%                     fclose(fout);
%                     delete(tmpFile);
%                 end
%             end
        case 'ECoG'
%             %% Analogue data, ecog
%             if doAnalog
%                 fprintf(1,'\t Parsing ECoG...\n');
%                 % Extract the sampling rate of the wavenform.
%                 Fs=ox.ReadEventsV(1,'ECoG',0,0,0,0,'ALL');
%                 Fs=ox.ParseEvInfoV(0,Fs,0);
%                 Fs=unique(Fs(9,:));
% %                 Fs=floor(Fs);
% %                 ox.SetGlobals(sprintf('WaveSF=%d',Fs));
% 
%                 % Unlike field potential data, we will store some ecog
%                 % data in the .mat file.  However, we will downsample
%                 % it to reduce file size.  We will compute an integer
%                 % decimation factor to bring us as close as we can to
%                 % the requested ecog sampling rate.  Then we will
%                 % design a low-pass Butterworth filter at the Nyquist
%                 % of the target sampling rate (with the specified
%                 % leeway to make sure no Nyquist effects) to prevent
%                 % aliasing.
%                 % Some sanity checking just in case no downsampling is
%                 % required!
%                 if Fs<=ecogFs
%                     doFilter=false;
%                     ecogFsTrue=Fs;
%                 else
%                     doFilter=true;
%                     decimateFactor=round(Fs./ecogFs);
%                     ecogFsTrue=Fs./decimateFactor;
%                     [filterB,filterA]=butter(butterN,...
%                         (1-ecogFilterLeeway).*2*ecogFsTrue/Fs,'low');
%                 end
%                 for c=1:numChannels
%                     ox.SetGlobalV('Channel',c);
%                     cnt=cnt+1;
%                     physiological.measures(cnt).name='ecog';
%                     physiological.measures(cnt).unit='';
%                     physiological.measures(cnt).values=[];
%                     physiological.measures(cnt).raw=ecogFile;
%                     physiological.measures(cnt).description= ...
%                         sprintf('ecog channel %d',c);
%                     physiological.measures(cnt).samplingRate=ecogFsTrue;
%                     % So, for each trial, read in the appropriate
%                     % segment of ECoG data, filter it, resample it at
%                     % the new sampling rate, then save it.
%                     for m=1:numTrials
%                         ox.SetGlobalV('T1',trialTimes(m));
%                         ox.SetGlobalV('T2',trialTimes(m+1));
%                         tmp=double(ox.ReadWavesV('ECoG'));
%                         if doFilter
%                             tmp=filter(filterB,filterA,tmp-mean(tmp));
%                             tmp=resample(tmp,1,decimateFactor);
%                         end
%                         physiological.data{m,cnt}=tmp;
%                     end
%                 end
%                 
%                 % Now go through and save the external file, which is
%                 % not broken down by trial and not resampled.  The
%                 % handling of this is identical to the handling of the
%                 % field potentials above.
%                 if doSave
%                     ox.SetGlobalV('Channel',0);
%                     thisWin=0;
%                     done=false;
%                     sz=0;
%                     maxval=0;
%                     fid=fopen(tmpFile,'Wb');
%                     while ~done
%                         ox.SetGlobalV('T1',thisWin*readWindow);
%                         ox.SetGlobalV('T2',(thisWin+1)*readWindow);
%                         waves=double(ox.ReadWavesV('ECoG'))';
%                         if isempty(waves)||(length(waves)==1&&isnan(waves))
%                             done=true;
%                         else
%                             maxval=max([maxval max(abs(waves))]);
%                             if thisWin==0
%                                 % OpenEx seems to have the first couple of values be an
%                                 % artefact.  Solve this by throwing out the first
%                                 % millisecond for the purpose of calculating the max.
%                                 maxval=max(max(abs(waves(round(1e-3*Fs)))));
%                                 fwrite(fid,waves,'float32');
%                             else
%                                 maxval = max([maxval max(max(abs(waves)))]);
%                                 fwrite(fid,waves,'float32');
%                             end
%                             sz=sz+numel(waves);
%                         end
%                         thisWin=thisWin+1;
%                     end
%                     fclose(fid);
%                     scale=maxval/32767;
%                     fin=fopen(tmpFile,'r');
%                     fout=fopen(ecogFile,'Wb','ieee-be');
%                     lsfwriteheader(fout,'numchannels',numChannels,'srate',Fs,...
%                         'scale',1,'encoding',LSF_ENC_LINEAR_16,...
%                         'info',['ecog from ' filename],...
%                         'datasize',sz*2);
%                     count=samplesToRead;
%                     while count==samplesToRead
%                         [waves,count]=fread(fin,samplesToRead,'float32');
%                         waves=round(waves./scale);
%                         fwrite(fout,waves,'int16');
%                     end
%                     fclose(fin);
%                     fclose(fout);
%                     delete(tmpFile);
%                 end
%             end
            
        case {'Epoc' 'Mark' 'Spks' 'Valu'}
            %% We do nothing with these:
            %   Epoc and Valu have already been handled.
            %   Mark is an internal timer variable for OpenEx
            %   Spks is redundant with Snip, and only stored for OpenEx
            %    reasons
        otherwise
            %% Error handling
            warning('ida:TTankHandling:newStore',...
                'Unknown store field ''%s'' in tank, ignored',eventNames{n});
    end
end

%% Final data saving and return values
fprintf(1,'\t Merging with stimulus information...\n');

results.physiological=physiological;
% If we were only asked to process a single file, then return its
% contents.  Otherwise, return the number of files processed.
ox.CloseTank; % close tank
fprintf(1,'\t\t Done!\n');

%% Final clean-up
ox.ReleaseServer; % disconnect server
close(Hf); % close the ActiveX hosting figure.

