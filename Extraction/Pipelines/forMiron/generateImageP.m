function [IMG] = generateImageP(MASK,R,N)
    idx = find(MASK);
    IMG = zeros(size(MASK));
    idx = idx(randperm(numel(idx)));
    idx = idx(1:N);
    for e = 1:N
        IMG(idx(e)) = 1;
    end
    IMG = bwdist(IMG) < R;
    
    %IMG = imdilate(IMG,strel('disk',R,0));
end