function [M,eR,wM] = wholeCotImageMask(I,sen,erValue)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make the mask for the cots

    %I = imerode(I,strel('disk',erValue,0));
    

    %I = rgb2gray(I);
    M = I > graythresh(I);
    wM = M;
    %M = imerode(M,strel('disk',erValue,0));

    %T = adaptthresh(I,sen,'ForegroundPolarity','bright');
    %M = imbinarize(I,T);
    %M = bwareaopen(M,10000);
    M = imclearborder(M);
    M = bwareaopenRange(M,[2000 100000]);
    R = regionprops(logical(M),'BoundingBox','Area','PixelIdxList');
    cidx = count([R.Area]);
    kidx = find(cidx < 2);
    M = zeros(size(M));
    for k = 1:numel(kidx)
        M(R(kidx(k)).PixelIdxList) = 1;
    end
    M = imfill(M,'holes');

    
    R = regionprops(logical(M),'PixelIdxList');
    eR = regionprops(logical(zeros(size(M))),'Centroid','BoundingBox','Area','PixelIdxList','Orientation','MajorAxisLength','MinorAxisLength');
    for e = 1:numel(R)
        tmpM = zeros(size(I));
        tmpM(R(e).PixelIdxList) = 1;
        %tmpM = imdilate(tmpM,strel('disk',erValue,0));
        eR(e) = regionprops(logical(tmpM),'Centroid','BoundingBox','Area','PixelIdxList','Orientation','MajorAxisLength','MinorAxisLength');
    end
    
end