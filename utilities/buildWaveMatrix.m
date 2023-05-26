function [mfWaves,afFs]=buildWaveMatrix(atResults,vnMeasureMap,varargin)
asWindowType='ratio'; % valid types: window, ratio, offset, full
afWindow=1.2;
abDoSmooth=true;
anSmoothN=5;

warnopts(assignopts(who,varargin{:}));

afFs=unique([atResults.physiological.measures(vnMeasureMap).samplingRate]);
if length(afFs)>1
    error('ida:buildWaveMatrix:samplingRate','Different sampling rates?!')
end

afIOI=1000*min(diff(atResults.stimulus.trialMarkers.trialTimes));

anDurationIndex=strncmpi(...
    {atResults.stimulus.stimuli_parameters.stimuli.specs.param.name},...
    'duration',8);
afStimDur=atResults.stimulus.stimuli_parameters.stimuli.specs.param(anDurationIndex).value;

switch lower(asWindowType)
    case 'ratio'
        afWindow=afStimDur*afWindow;
    case 'offset'
        afWindow=afStimDur+afWindow;
    case 'window'
        %afWindow=afWindow;
    case 'full'
        afWindow=afIOI;
    otherwise
        error('ida:buildWaveMatrix:windowType','Unknown window type');
end

if afWindow>afIOI
    afWindow=afIOI;
    warning('ida:buildWaveMatrix:windowTooLong',...
        'Desired window exceeds minimum inter-stimulus interval, truncating');
end
anNumSamps=ceil(afWindow*1e-3*afFs);

mfWaves=nan(anNumSamps,length(vnMeasureMap),size(atResults.stimulus.order,1));

for n=1:length(vnMeasureMap)
    for m=1:size(atResults.stimulus.order,1)
        mfWaves(:,n,m)=atResults.physiological.data{m,vnMeasureMap(n)}...
            (1:anNumSamps);
        if abDoSmooth
            vfSMOOTH=ones(1,anSmoothN)./anSmoothN;
            mfWaves(:,n,m)=filter(vfSMOOTH,1,mfWaves(:,n,m));
        end
        
        
    end
end
mfWaves=mfWaves-mean(mfWaves(:));