function [probMap] = makeProbMap(tX,tY,N,sigma)

    NX = N(1);
    NY = N(2);
    %probMap = zeros(NY,NX);
    ptX = tX - min(tX(:));
    ptX = ptX / max(ptX(:));
    ptY = tY - min(tY(:));
    ptY = ptY / max(ptY(:));
    dX = discretize(ptX,linspace(0,1,NX));
    dY = discretize(ptY,linspace(0,1,NY));
    
    [py,px] = ndgrid(1:NY,1:NX);
    D = pdist2([dX dY],[px(:) py(:)]);
    
    
    probMap = sum(exp(-(D*sigma^-1).^2),1);
    probMap = reshape(probMap,N);
    
    %probMap = imfilter(probMap,fspecial('gaussian',61,31),'replicate');
    probMap = imfilter(probMap,fspecial('disk',51),'replicate');
    
    probMap = probMap / sum(probMap(:));
    %{
    didx = sub2ind(size(probMap),dY,dX);
    
    
    for p = 1:numel(didx)
        probMap(didx(p)) = probMap(didx(p)) + 1;
    end
    D = squareform(pdist([dY dX]));
    
    %probMap = imdilate(probMap,strel('disk',21,0));
    probMapDist = bwdist(probMap);
    pM = exp(-(probMapDist*sigma^-1).^2);
    %pM = imfilter(pM,fspecial('gaussian',61,31),'replicate');
    pM = pM / sum(pM(:));
    %}
end