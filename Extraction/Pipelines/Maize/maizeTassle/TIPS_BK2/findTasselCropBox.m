function [BOX] = findTasselCropBox(I,reSZ,PER)

    Lab = rgb2lab(I);
    %Lab = rgb2gray(I);


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
    oldSZ = size(tasselM);
    tasselM = imresize(tasselM,.25);
    tasselM = imclearborder(tasselM);
    tasselM = imclose(tasselM,strel('disk',25,0));
    tasselM = imclearborder(tasselM);
    tasselM = imresize(tasselM,oldSZ);
    tasselM = logical(tasselM);

    tasselObjects = regionprops(tasselM,'Centroid','BoundingBox');
    centerPoint = [size(I,2) size(I,1)]/2;
    % fill holes
    dist = [];
    for e = 1:numel(tasselObjects)
        dist(e) = norm(centerPoint - tasselObjects(e).Centroid);
    end
    [~,midx] = min(dist);
    BOX = tasselObjects(midx).BoundingBox;


    BOXExpand = BOX(3:4)*PER/2;
    BOX(1:2) = BOX(1:2) - BOXExpand;
    BOX(3:4) = BOX(3:4) + 2*BOXExpand;


    BOX(1) = max(BOX(1),1);
    BOX(2) = max(BOX(2),1);

    MAX_RIGHT = BOX(1) + BOX(3);
    if MAX_RIGHT > size(I,2)
        BOX(3) = size(I,2) - BOX(1);
    end

    MAX_BOTTOM = BOX(2) + BOX(4);
    if MAX_BOTTOM > size(I,1)
        BOX(4) = size(I,1) - BOX(2);
    end

end