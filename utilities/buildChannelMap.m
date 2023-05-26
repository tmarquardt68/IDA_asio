function mnChannelMap=buildChannelMap(atResults,asMeasure)

if strncmpi(asMeasure,'ecog',4)
    mnChannelMap = [3 14 11 6;
                    2 15 10 7;
                    4 13 12 5;
                    1 16 9 8];
elseif strncmpi(asMeasure,'field',5) || strncmpi(asMeasure,'spike',5)
    mnChannelMap = [6 11 3 14 1 16 2 15 5 12 4 13 7 10 8 9];
else
    error('ida:buildChannelMap:measureType','Unknown measure type');
end