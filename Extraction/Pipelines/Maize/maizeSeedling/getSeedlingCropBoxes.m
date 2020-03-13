function [cropPackage,I] = getSeedlingCropBoxes(fileName,coneDetNetL_region,coneDetNetR_region,coneDetNetL_fine,coneDetNetR_fine,oPath,rPath,disp)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fileName :- the name of the file to process
    % coneDetNetL_region :- CNN for detecting the left of the cone region
    % coneDetNetR_region :- CNN for detecting the right of the cone region - rough shot
    % coneDetNetL_fine :- CNN for fine detection on left
    % coneDetNetR_fine :- CNN for fine detection on right
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mkdir(oPath);
    [pth,nm,ext] = fileparts(fileName);
    fileList = {};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % left and right points
    LP = [];
    RP = [];
    % amount of buffer needed in pixels from the edge of the images for
    % detection
    CROP_PAD = 200;
    % rough shot down sample
    DS = 10;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['start: generating domains for multi-scale detection.\n']);
    [regionDetectionDomain] = generateMultiResSampleDomain([800 800],[1],[20 20]);
    [fineDetectionDomain] = generateMultiResSampleDomain([800 800],[.1 .5 1],[100 100]);
    fprintf(['end: generating domains for multi-scale detection.\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['start: reading image.\n']);
    I = double(imread(fileName));
    I(:,(end-70):end,:) = [];
    % make gray scale image
    G = rgb2gray(I/255);
    fprintf(['start: reading image.\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['start: finding the meta data sheet data.\n']);
    pad = 50;
    Lab = rgb2lab(I/255);
    RED = Lab(:,:,2) > 25.5;
    clear Lab % converse memory
    RED = imclose(RED,strel('square',51));
    RED = bwlarge(RED);
    R = regionprops(RED);
    box = R(1).BoundingBox;
    box(1:2) = box(1:2) - pad;
    box(3:4) = box(3:4) + 2*pad;
    fprintf(['end: finding the meta data sheet data.\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['start: generating the domain for rough scan.\n']);
    NP = [round(size(I,1)/DS) round(size(I,2)/DS)];
    [s2,s1] = ndgrid(linspace(CROP_PAD+1,size(G,1)-CROP_PAD,NP(1)),linspace(CROP_PAD+1,size(G,2)-CROP_PAD,NP(2)));
    S = [s1(:) s2(:)];
    fprintf(['end: generating the domain  for rough scan.\n']);

    
    % init the view with the pad value removed from edges
    toView = G;
    for r = 1:4
        toView(1:CROP_PAD,:) = [];
        toView = imrotate(toView,90);
    end
    sizeToView = size(toView);
    clear toView;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['start: rough scan left.\n']);tic
    func_rough_left = @(X)coneDetNetL_region.predict(X);
    [toL] = sampleApply(G,S,regionDetectionDomain,func_rough_left);
    fprintf(['end: rough scan left.' num2str(toc) '\n']);
    fprintf(['start: rough scan right.\n']);tic
    func_rough_right = @(X)coneDetNetR_region.predict(X);
    [toR] = sampleApply(G,S,regionDetectionDomain,func_rough_right);
    fprintf(['end: rough scan right.' num2str(toc) '\n']);
    toL = squeeze(toL)';     
    toL = reshape(toL(:,2),NP);
    toR = squeeze(toR)';
    toR = reshape(toR(:,2),NP);
    toR = imresize(toR,sizeToView);
    toL = imresize(toL,sizeToView);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %RRegion = bwareaopen(toR > .90,3000);
    %LRegion = bwareaopen(toL > .90,3000);
    % filter image and make "dots" image at max intensities for each region
    RRegionF = imfilter(toR,fspecial('disk',31),'replicate');
    LRegionF = imfilter(toL,fspecial('disk',31),'replicate');
    LRegion = bindVec(LRegionF);
    RRegion = bindVec(RRegionF);
    RRegionM = bwareaopen(RRegion > .5,1000);
    LRegionM = bwareaopen(LRegion > .5,1000);
    RRegionM = imclearborder(RRegionM);
    LRegionM = imclearborder(LRegionM);
    
    RRR = regionprops(RRegionM,RRegionF,'PixelIdxList','Area','MaxIntensity');
    LLL = regionprops(LRegionM,LRegionF,'PixelIdxList','Area','MaxIntensity');
    RRegion = zeros(size(RRegion));
    LRegion = zeros(size(LRegion));
    
    [~,sidx1] = sort([LLL.Area],'descend');
    [~,sidx2] = sort([RRR.Area],'descend');
    LLL = LLL(sidx1(1:3));
    RRR = RRR(sidx2(1:3));
    
    for r = 1:numel(RRR)
        didx = find(RRegionF(RRR(r).PixelIdxList) == RRR(r).MaxIntensity);
        RRegion(RRR(r).PixelIdxList(didx)) = 1;
    end
    for r = 1:numel(LLL)
        didx = find(LRegionF(LLL(r).PixelIdxList) == LLL(r).MaxIntensity);
        LRegion(LLL(r).PixelIdxList(didx)) = 1;
    end
    RRegion = imdilate(RRegion,strel('disk',31));
    LRegion = imdilate(LRegion,strel('disk',31));
    rightPT = [];
    leftPT = [];
    [rightPT(:,2),rightPT(:,1)] = find(RRegion);
    ridx = find(RRegion==1);
    lidx = find(LRegion==1);
    [leftPT(:,2),leftPT(:,1)] = find(LRegion);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['start: fine scan left.\n']);tic
    rightPT = bsxfun(@plus,rightPT,[CROP_PAD,CROP_PAD]);
    leftPT = bsxfun(@plus,leftPT,[CROP_PAD,CROP_PAD]);
    func_find_left = @(X)coneDetNetL_fine.predict(X);
    func_find_right = @(X)coneDetNetR_fine.predict(X);
    [toL_FINE] = sampleApply(G,leftPT,fineDetectionDomain,func_find_left);
    fprintf(['end: fine scan left.' num2str(toc) '\n']);
    fprintf(['start: fine scan right.\n']);tic
    [toR_FINE] = sampleApply(G,rightPT,fineDetectionDomain,func_find_right);
    fprintf(['end: fine scan right.' num2str(toc) '\n']);
    toRP = squeeze(toR_FINE(1,2,:));
    toLP = squeeze(toL_FINE(1,2,:));
    %toRP2 = coneDetNetR_fine.predict(toR_FINE2,'MiniBatchSize',2048);
    %toLP = coneDetNetL_fine.predict(toL_FINE,'MiniBatchSize',2048);
    toR_FINE = zeros(size(toR));
    toR_FINE(ridx) = toRP(:);
    toL_FINE = zeros(size(toL));
    toL_FINE(lidx) = toLP(:);
    %out = imresize(out,.25);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [~,sidx1] = sort([LLL.Area],'descend');
    [~,sidx2] = sort([RRR.Area],'descend');
    LLL = LLL(sidx1(1:3));
    RRR = RRR(sidx2(1:3));
    for p = 1:3
        tmpLoc = [];
        MSK = zeros(sizeToView);
        MSK(LLL(p).PixelIdxList) = 1;
        MSK = MSK.*toL_FINE;
        %MSK = logical(MSK > 0);
        tmp = regionprops(logical(MSK>0),toL_FINE,'MeanIntensity','MaxIntensity');
        [~,midx] = max(tmp.MeanIntensity);
        %LP(p,:) = tmp(midx).WeightedCentroid;
        [tmpLoc(:,2),tmpLoc(:,1)] = find(MSK == tmp(midx).MaxIntensity);
        [LP(p,:)] = mean(tmpLoc,1);
    end
    for p = 1:3
        tmpLoc = [];
        MSK = zeros(sizeToView);
        MSK(RRR(p).PixelIdxList) = 1;
        MSK = MSK.*toR_FINE;
        %MSK = logical(MSK > 0);
        tmp = regionprops(logical(MSK>0),toR_FINE,'MeanIntensity','MaxIntensity');
        [~,midx] = max(tmp.MeanIntensity);
        [tmpLoc(:,2),tmpLoc(:,1)] = find(MSK == tmp(midx).MaxIntensity);
        [RP(p,:)] = mean(tmpLoc,1);
    end
    [~,kidx] = sort(LP(:,1));
    LP = LP(kidx,:);
    [~,kidx] = sort(RP(:,1));
    RP = RP(kidx,:);

    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    


    toR_FINE = padarray(toR_FINE,[CROP_PAD,CROP_PAD],0,'both');
    out = flattenMaskOverlay(I/255,toR_FINE>.7,.7,'r');
    clear toR_FINE
    
    
    toL_FINE = padarray(toL_FINE,[CROP_PAD,CROP_PAD],0,'both');
    out = flattenMaskOverlay(out,toL_FINE>.7,.7,'b');
    clear toL_FINE

    
    
    %oSTORE{e} = out;



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    




    SEEDLING_BOX = {};
    SEEDLING_BOX{1} = [1 box(2)+box(4) .5*(LP(2,1) + RP(1,1)+2*CROP_PAD) mean([LP(1,2) RP(1,2)])+CROP_PAD - (box(2)+box(4))];
    SEEDLING_BOX{2} = [SEEDLING_BOX{1}(1)+SEEDLING_BOX{1}(3) box(2)+box(4) .5*(LP(2,1) + RP(3,1)+2*CROP_PAD) - SEEDLING_BOX{1}(3) mean([LP(2,2) RP(2,2)])+CROP_PAD - (box(2)+box(4))];
    SEEDLING_BOX{3} = [SEEDLING_BOX{2}(1)+SEEDLING_BOX{2}(3) box(2)+box(4) size(I,2) - .5*(LP(3,1) + RP(2,1)+2*CROP_PAD) - 1  mean([LP(3,2) RP(3,2)])+CROP_PAD - (box(2)+box(4))];


    DELTA = (RP - LP);
    ANG = -atan2(DELTA(:,2),DELTA(:,1));
    CENTER_POINT = .5*(RP + LP + 2*CROP_PAD);
    rBOXpoints = {};
    for p = 1:numel(SEEDLING_BOX)
        rBOXpoints{p} = box2points(SEEDLING_BOX{p});
        rBOXpoints{p} = bsxfun(@minus,rBOXpoints{p},CENTER_POINT(p,:));
        rBOXpoints{p} = [[cos(ANG(p)) sin(ANG(p))];[-sin(ANG(p)) cos(ANG(p))]]*rBOXpoints{p}';
        rBOXpoints{p} = rBOXpoints{p}';
        rBOXpoints{p} = bsxfun(@plus,rBOXpoints{p},CENTER_POINT(p,:));
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if disp
        imshow(out,[]);
        drawnow
        hold on
        plot(LP(:,1)+CROP_PAD,LP(:,2)+CROP_PAD,'r*')
        plot(RP(:,1)+CROP_PAD,RP(:,2)+CROP_PAD,'b*')
        for p = 1:3
            plot([LP(p,1),RP(p,1)]+CROP_PAD,[LP(p,2),RP(p,2)]+CROP_PAD,'y');
        end
        toc
        rectangle('Position',box,'EdgeColor','r')
        for p = 1:numel(SEEDLING_BOX)
            rectangle('Position',SEEDLING_BOX{p},'EdgeColor','g');
        end



        for p = 1:numel(rBOXpoints)
            plot(rBOXpoints{p}(:,1),rBOXpoints{p}(:,2),'b')
        end
        hold off
        drawnow
    end
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        
    cropPackage.out = out;
    clear out % for memory
    cropPackage.LP = LP;
    cropPackage.RP = RP;
    cropPackage.rBOXpoints = rBOXpoints;
    cropPackage.SEEDLING_BOX = SEEDLING_BOX;
    for e = 1:numel(SEEDLING_BOX)
        cropPackage.croppedImage{e} = imcrop(I,cropPackage.SEEDLING_BOX{e})/255;  
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        
    for e = 1:numel(cropPackage.croppedImage)
        fileList{end+1} = [oPath nm '{stage1croppedImage_' num2str(e) '}.tif'];
        imwrite(cropPackage.croppedImage{e},fileList{end});
    end
        
    pushToiRods(rPath,fileList);
        
        
end

%{



cropPackage =
getSeedlingCropBoxes(FileList{9},coneDetNetL_region,coneDetNetR_region,coneDetNetL_fine,coneDetNetR_fine,true);

oPath = '/mnt/tetra/nate/tmpSeedling/';
rPath = '';
for w = 17:numel(wholeFileList)
    smartMain_new_ver1(wholeFileList{w},coneDetNetL_region,coneDetNetR_region,coneDetNetL_fine,coneDetNetR_fine,GMModel,GMModel2,oPath,rPath,true);
end

%}



