function [riseDomain,recoverDomain] = generateBinaryDomain(data_in,threshold)

    domainWave = data_in(1,:) > threshold;
    domainWave = bwlarge(domainWave);
    
    [~,peakIDX] = max(data_in(1,:));
    
    
    riseDomain = domainWave;
    riseDomain(peakIDX+1:end) = 0;
     
    
    recoverDomain = domainWave;
    recoverDomain(1:peakIDX) = 0;
    
end