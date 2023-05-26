function vnMeasureMap=buildChannelToMeasureMap(atResults,asMeasureName)

vnMeasureIdx=find(strncmpi({atResults.physiological.measures.name},...
    asMeasureName,length(asMeasureName)));

vnMeasureMap=zeros(size(vnMeasureIdx));

for n=1:length(vnMeasureMap)
    anChannelIdx=strfind(lower(...
        atResults.physiological.measures(vnMeasureIdx(n)).description),'channel');
    if strncmp(atResults.physiological.measures(vnMeasureIdx(n)).description(anChannelIdx+7),...
            ' ',1)
        anOffset=8;
    else
        anOffset=7;
    end
    vsDescriptionRemainder=...
        atResults.physiological.measures(vnMeasureIdx(n)).description(...
        anChannelIdx+anOffset:end);
    vnTokenIdx=strfind(vsDescriptionRemainder,' ');
    if ~isempty(vnTokenIdx)
        vnMeasureMap(str2double(vsDescriptionRemainder(1:vnTokenIdx(1)-1)))=...
            vnMeasureIdx(n);
    else
        vnMeasureMap(str2double(vsDescriptionRemainder))=...
            vnMeasureIdx(n);
    end
end