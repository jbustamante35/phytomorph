function [tassel_mask] = makeTasselMaskBlur(raw, reSZ, hole_max)
%MAKETASSELMASKBLUR Summary of this function goes here
%   Detailed explanation goes here
 
toProcess = raw;
oldSZ = size(toProcess);
toProcess = imresize(toProcess,1/reSZ);
toProcess = imfilter(toProcess,fspecial('disk',27),'replicate');
toProcess = imdilate(toProcess,strel('line',101,0));
toProcess = imfilter(toProcess,fspecial('gaussian',[31 31],7),'replicate');
toProcess = imresize(toProcess,oldSZ);

tassel = toProcess - raw;
tassel = bindVec(tassel);
tassel_mask = tassel > graythresh(tassel); 

% Fill holes smaller than hole_max pixels and connect components by closing
filled = imfill(tassel_mask, 'holes');
holes = filled & ~tassel_mask;
bigholes = bwareaopen(holes, hole_max);
smallholes = holes & ~bigholes;
tassel_mask = tassel_mask | smallholes;

% 
tassel_mask = imclose(tassel_mask, strel('disk', 10));

end

