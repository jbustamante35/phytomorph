function [tasselMask] = tasselMask2OriginalSize(tasselMask, BOX, orig)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    
    tidx = [];
    [tidx(:,2),tidx(:,1)] = find(tasselMask);

    if ~isempty(tidx)
        tidx(:,1) = tidx(:,1) + round(BOX(1));
        tidx(:,2) = tidx(:,2) + round(BOX(2));
    end
    T = sub2ind([size(orig,1) size(orig,2)],tidx(:,2),tidx(:,1));
    tasselMask = zeros(size(orig,1),size(orig,2));
    tasselMask(T) = 1;
    tasselMask = logical(tasselMask);
    %tasselM = imfill(tasselM,'holes');
    %tasselM = imopen(tasselM,strel('disk',7,0));
end

