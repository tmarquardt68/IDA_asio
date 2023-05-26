function mfWaves=buildEcogWaveMatrix(atResults,anNumSamps,vnMeasureMap)

mfWaves=nan(anNumSamps,length(vnMeasureMap),size(results.stimulus.order,1));

for n=1:length(vnMeasureMap)
    for m=1:size(results.stimulus.order,1)
        mfWaves(:,n,m)=atResults.physiological.data{m,vnMeasureMap(n)}...
            (1:anNumSamps);
    end
end