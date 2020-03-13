function [tasselM] = thresholdTasselImage(img, reSZ)
%
%
    Lab = rgb2lab(img);

    isColor = mean(abs(Lab(:,:,2:3)),3);
    %isColor = max(abs(Lab(:,:,2:3)),[],3);
    raw = imcomplement(isColor);
    %raw2 = imcomplement(Lab(:,:,1));
    %raw3 = mean(cat(3,raw,raw2),3);
    %raw = raw3;



    toProcess = raw;
    oldSZ = size(toProcess);
    toProcess = imresize(toProcess,1/reSZ);
    toProcess = imfilter(toProcess,fspecial('gaussian',[31 31],7),'replicate');
    toProcess = imdilate(toProcess,strel('line',101,0));
    toProcess = imfilter(toProcess,fspecial('gaussian',[31 31],7),'replicate');
    toProcess = imresize(toProcess,oldSZ);




    tassel = toProcess - raw;
    %tassel = imcomplement(raw);
    tassel = bindVec(tassel);
    tasselM = tassel > graythresh(tassel);
    tasselM = imclose(tasselM,strel('disk',5,0));

    tasselM = imclearborder(tasselM);

    tasselM = bwlarge(tasselM);
    tasselM = imclose(tasselM,strel('disk',5,0));


end
   