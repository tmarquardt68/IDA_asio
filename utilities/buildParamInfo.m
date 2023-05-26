function [vfParamList,vfParamSequence]=buildParamInfo(atResults,asParamName)

anParamIdx=find(strncmpi({atResults.stimulus.parameters.name},...
    asParamName,length(asParamName)));
vfParamList=atResults.stimulus.parameters(anParamIdx).values;
if isempty(vfParamList)
    vfParamSequence=atResults.stimulus.order(:,anParamIdx);
else
    if any(atResults.stimulus.order(:,anParamIdx)==0)
        vfParamList(end+1)=NaN;
        tmp=atResults.stimulus.order(:,anParamIdx);
        tmp(tmp==0)=length(vfParamList);
        vfParamSequence=vfParamList(tmp);
        vfParamList=vfParamList(1:end-1);
    else
        vfParamSequence=vfParamList(atResults.stimulus.order(:,anParamIdx));
    end
end