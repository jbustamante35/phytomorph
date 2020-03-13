function [truePoints,falsePoints] = generateRegionPoints(trueMask,dilateSize,re,modV,N,DIST)
    sz = size(trueMask);    
    trueMask = imresize(trueMask,re);
    trueMask = imdilate(trueMask,strel('disk',dilateSize,0));
    trueMask = imresize(trueMask,sz);
    trueMask = trueMask > 0;
    if DIST > 0
        D = bwdist(trueMask);
        falseMask = D < DIST & ~trueMask;
        
    else
        falseMask = trueMask == 0;
    end
    
    [truePoints(:,2),truePoints(:,1)] = find(trueMask);
    [falsePoints(:,2),falsePoints(:,1)] = find(falseMask);
    
    
    
    
    idx = mod(truePoints(:,1),modV) == 0 & mod(truePoints(:,2),modV) == 0;
    truePoints = truePoints(idx,:);
    
    
    if N < 0
        N = size(truePoints,1);
    end
    
    falsePoints = falsePoints(randperm(size(falsePoints,1)),:);
    N = min(N,size(falsePoints,1));
    falsePoints = falsePoints(1:N,:);
end





















