function badTrials = filterTrials(trialsData, minBadChannels, minBadPercentage)

badTrials = [];
for iTrial = 1:length(trialsData)
    
    badPoints = sum(trialsData(iTrial).BadSig > 30) > minBadChannels;
    if sum(badPoints)/size(trialsData(iTrial).BadSig, 2)*100 > minBadPercentage
        
        badTrials(end+1) = iTrial;
        
    end
    
end