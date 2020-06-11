function [out,colorSquare,totalCount] = measureGeneralSeeds(I,probMap)
    smallObjectThreshold = 50;
    localT = false;
    
    if ~localT
        % run ostu on value channel
        Mask = probMap > graythresh(probMap);
    else
        th = adaptthresh(probMap,'NeighborhoodSize',601);
        Mask = imbinarize(probMap,th);
    end
    
    
    % remove small objects
    Mask = bwareaopen(Mask,smallObjectThreshold);
    
    % fill in the object - added for wilson controls of coins
    %M = imfill(M,'holes');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Measure with region props.\n']);
    R = regionprops(Mask,'Area','PixelIdxList','MajorAxisLength','MinorAxisLength','Centroid','Eccentricity');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Running counting method.\n']);
    % count on: area & major axis & minor axis
    cidx = count([R.Area]);
    cidx1 = count([R.MajorAxisLength]);
    cidx2 = count([R.MinorAxisLength]);
    fidx = find(cidx==1 & cidx1==1 & cidx2==1);
    singleMask = zeros(size(Mask));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Filling single seed mask.\n']);
    for e = 1:numel(fidx)
        singleMask(R(fidx(e)).PixelIdxList) = 1;
    end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Sampling color channels over all seeds.\n']);
    meanRGB = zeros(numel(fidx),size(I,3));
    stdRGB = meanRGB;
    % for each color channel in image
    for k = 1:size(I,3)
        tmp = I(:,:,k);
        % for each seed
        for e = 1:numel(fidx)
            % sample color channel at seed
            cl = double(tmp(R(fidx(e)).PixelIdxList));
            meanRGB(e,k) = mean(cl);
            stdRGB(e,k) = std(cl);
            seedColorSamples(e).colorData(:,k) = cl;
            seedColorSamples(e).Centroid = R(fidx(e)).Centroid;
        end
    end
    
    
    allColors = [];
    for e = 1:numel(fidx)
        allColors = cat(1,allColors,seedColorSamples(e).colorData);
    end
    
    alpha = 1;
    steps = 5;
    [U,E,L] = PCA_FIT_FULLws(allColors,1);
    LS = alpha*linspace(-L.^.5,L.^.5,steps);
    colorSquare = [];
    for s = 1:numel(LS)
        tmpC = [];
        for k = 1:3
            tmpC(:,:,k) = (U(k)+E(k)*LS(s))*ones(101,101);
            eyeColor(s,k) = U(k)+E(k)*LS(s);
        end
        colorSquare = [colorSquare;tmpC];
    end
    colorSquare = colorSquare/255;
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make clump image
    clumpMask = zeros(size(singleMask));
    clidx = find(cidx~=1);
    totalCount = cidx;
    totalCount(totalCount > 100) = 0;
    totalCount = sum(totalCount);
    for e = 1:numel(clidx)
        clumpMask(R(clidx(e)).PixelIdxList) = 1;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['Make overlay.\n']);
    out = flattenMaskOverlay(I,logical(singleMask),.6,'r');
    out = flattenMaskOverlay(out,logical(clumpMask),.6,'g');
    
end