function [wData] = fixData(wData)
    
    for e = 1:size(wData,1)
        nidx = find(~isnan(wData(e,:)));
        if ~isempty(nidx)
            wData(e,1:nidx(1)) = 0;
        end
    end
    ridx = sum(isnan(wData),2) > 12;
    wData = knnimpute(wData, 11);
    wData(ridx,:) = NaN;
end