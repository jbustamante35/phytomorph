function [pMap] = pMapFromPoints(pts,template)
    pMap = zeros(size(template));
    idx = sub2ind(size(template),pts(:,1),pts(:,2));
    for e = 1:numel(idx)
        pMap(idx) = pMap(idx) + 1;
    end
end