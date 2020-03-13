FilePath = '/home/nate/Vein_pattern_images/';
%FilePath = '/mnt/spaldingdata/nate/makotoReturn/';
%FilePath = '/home/nate/Vein_pattern_images/OMY183_Col/';
cFileList = {};
FileExt = {'tif'};
cFileList = fdig(FilePath,cFileList,FileExt,1);
%% look at the image sizes
for e = 1:numel(cFileList)
    I = imread(cFileList{e});
    SZ(e,:) = size(I);
end

%% make model version 2
close all
T = [];
for e = 1:numel(cFileList)
    fileName = cFileList{e};
    I = double(imread(fileName))/255;
    I = rgb2gray(I);
    
    I = imresize(I,[1024 1440]);
    
    
    th = graythresh(I);
    M = I > th;
    %M = bwareaopenRange(M,[4000 10000]);
    M = bwareaopen(M,5000);
    M = imclearborder(M);
    R = regionprops(M,'Orientation','Centroid','BoundingBox','Area','Perimeter','MajorAxisLength','MinorAxisLength');
    c1 = find(count([R.Area]) == 1);
    imshow(I,[]);
    hold on
    for r = 1:numel(R)
        R(r).fileName = fileName;
    end
    for r = 1:numel(R)
        rectangle('Position',R(r).BoundingBox,'EdgeColor','r');
    end
    for r = 1:numel(c1)
        rectangle('Position',R(c1(r)).BoundingBox,'EdgeColor','g');
    end
    if ~isempty(R)
        T = [T;R];
    end
    drawnow
    hold off
end
%% count the total number for template building
close all
fidx = find(count([T.Area])==1 & count([T.MajorAxisLength])==1);
fidx = find(count([T.Area])==1);
%fidx = find(count([T.Area])==1 & count([T.MajorAxisLength].*[T.MinorAxisLength].^-1)==1);
for e = 1:numel(fidx)
    I = double(imread(T(fidx(e)).fileName))/255;
    I = rgb2gray(I);
    I = imresize(I,[1024 1440]);
    size(I)
    imshow(I,[]);
    hold on
    rectangle('Position',T(fidx(e)).BoundingBox,'EdgeColor','r');
    drawnow
    hold off   
end
%% crop each template
close all
fidx = find(count([T.Area])==1);
cropDims = 1.5*[[T(fidx).MajorAxisLength];[T(fidx).MinorAxisLength]]';
cropLoc = reshape([T(fidx).Centroid],[2 numel(fidx)])';
cropTheta = (pi/180)*[T(fidx).Orientation];
for b = 1:size(cropDims,1)
    I = double(imread(T(fidx(b)).fileName))/255;
    I = rgb2gray(I);
    I = imresize(I,[1024 1440]);
    th = graythresh(I);
    M = double(I > th);
    tmpI = rimcrop(I,cropDims(b,2),cropDims(b,1),cropTheta(b),cropLoc(b,1:2));
    tmpM = rimcrop(M,cropDims(b,2),cropDims(b,1),cropTheta(b),cropLoc(b,1:2));
    tmpM = logical(tmpM);
    %tmpM = imclearborder(tmpM);
    tmpM = bwlarge(tmpM);
    R = regionprops(tmpM,'Orientation','Centroid','BoundingBox','Area','Perimeter','MajorAxisLength','MinorAxisLength');
    T(b).Area = R.Area;
    T(b).MajorAxisLength = R.MajorAxisLength;
    T(b).MinorAxisLength = R.MinorAxisLength;
    T(b).Perimeter = R.Perimeter;
    
    out = flattenMaskOverlay(bindVec(tmpI),tmpM);
    img{b} = out;
    imshow(out,[]);
    drawnow
end
%% select X bases for templates
baseN = 5;
featureStack = [[T.Area];[T.MajorAxisLength];[T.MinorAxisLength];[T.Perimeter]];
uF = mean(featureStack,2);
deltaFS = featureStack - uF;
[deltaFS,sidx] = sort(sum(deltaFS.*deltaFS,1));
for b = 1:baseN
    imshow(img{sidx(b)},[]);
    waitforbuttonpress
end
%%
for e = 1:numel(cFileList)
    I = imread(cFileList{e});
    imshow(I,[]);
    drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
mag = -.4;
mag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% value for spuring the skeleton
spurAmount = 3;
% value to measure the cost for spurs that connect branch to end
snipAmount = 11;
% amount to dilate the skeleton
skeletonDilateAmount = 3;
% amount to erode the mask
totalMaskErodeAmount = 3;
% n-hood size for threshold
nsz = 41;
% sensetivity for threshold
sen = .40;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SEN = linspace(.35,.45,8);
SEN = .4;
for s = 1:numel(SEN)

    alignedPair = {};
    for e = 1:numel(cFileList)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % extract the single cots only from collection
        [pth,nm,ext] = fileparts(cFileList{e});
        I = double(imread(cFileList{e}))/255;
        I = rgb2gray(I);
        M = I > graythresh(I);
        M = bwareaopen(M,10000);
        M = imclearborder(M);
        R = regionprops(logical(M),'BoundingBox','Area','PixelIdxList');
        cidx = count([R.Area]);
        kidx = find(cidx < 2);
        M = zeros(size(M));
        for k = 1:numel(kidx)
            M(R(kidx(k)).PixelIdxList) = 1;
        end
        M = imfill(M,'holes');
        R = regionprops(logical(M),'Centroid','BoundingBox','Area','PixelIdxList','Orientation','MajorAxisLength','MinorAxisLength');


        if numel(R) > 1
            try
            close all
            [alignedPairMasks{1,e},alignedPairImage{1,e}] = alignMasks_ver2(I,M,R);
            catch ME
                ME
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% construct the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelSZ = [];
cnt = 1;
TOT = 1;
for e = 1:TOT
    for i = 1:size(alignedPairMasks{e},1)
        for j = 1:size(alignedPairMasks{e},2)
            tmp = size(alignedPairMasks{e}{i,j});
            if ~all(tmp==0)
                modelSZ(cnt,:) = tmp;
                cnt = cnt + 1;
            end
        end
    end
end
mSize = max(modelSZ,[],1);
mSize = mSize(1:2);
fidx = find(mod(mSize(1:2),2) == 0);
mSize(fidx) = mSize(fidx) + 1;
cnt = 1;
template = [];
% resize the data to fit the largest
for e = 1:TOT
    for i = 1:size(alignedPairMasks{e},1)
        for j = 1:size(alignedPairMasks{e},2)
            tmp = size(alignedPairMasks{e}{i,j});
            if ~all(tmp==0)
                deltaSZ = mSize - tmp(1:2);
                fidx = find(mod(deltaSZ,2) == 1);
                deltaSZ(fidx) = deltaSZ(fidx) + 1;
                deltaSZ = deltaSZ / 2;
                
                tmp = padarray(alignedPairMasks{e}{i,j}(:,:,1),deltaSZ,0,'both');
                tmp = imresize(tmp,mSize);
                
                template(:,:,cnt) = tmp;
                cnt = cnt + 1;
                
                template(:,:,cnt) = flip(tmp,2);
                cnt = cnt + 1;
                
                tmp = padarray(alignedPairMasks{e}{i,j}(:,:,2),deltaSZ,0,'both');
                tmp = imresize(tmp,mSize);
                template(:,:,cnt) = tmp;
                cnt = cnt + 1;
                
                template(:,:,cnt) = flip(tmp,2);
                cnt = cnt + 1;
                
                size(tmp)
            end
        end
    end
end
%%
close all
modelThreshold = .65;
imshow(mean(template,3) > modelThreshold,[]);
modelTemplate = double(mean(template,3) > modelThreshold);
scale = 1.8;
szM = size(modelTemplate);
modelTemplate = imresize(modelTemplate,scale);
R = regionprops(logical(modelTemplate));
modelTemplate = imcrop(modelTemplate,R(1).BoundingBox);
PAD = round((szM - size(modelTemplate))/2);
modelTemplate = padarray(modelTemplate,[PAD(1) PAD(1)],0,'both');
imshow(modelTemplate,[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% align to model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
alignedPair = {};
sen = .8;
dataPoint = {};
oPath = '/mnt/spaldingdata/nate/makotoReturn/';
mkdir(oPath);
for e = 1:numel(cFileList)
    

    toPath = [oPath nm filesep];
    mkdir(toPath);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % process the image and make mask
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fileName = cFileList{e};
    [pth,nm,ext] = fileparts(fileName);
    I = double(imread(fileName))/255;
    I = rgb2gray(I);
    [M,R,wM] = wholeCotImageMask(I,sen,0);
    wM = imclose(wM,strel('disk',5,0));


    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % clean up the model - resize and control padding
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PAD = 10;
    initResize = .75;
    useModel = imresize(modelTemplate,initResize);
    uR = regionprops(useModel > .8,'BoundingBox');
    useModel = imcrop(useModel,uR.BoundingBox);
    useModel = padarray(useModel,[PAD PAD],0,'both');
    modelDT = bwdist(~useModel);
    modelMask = useModel;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make a border mask
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    borderPad = 101;
    borderMask = ones(size(wM));
    for rot = 1:4
        borderMask(1:borderPad,:) = 0;
        borderMask = imrotate(borderMask,90);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the cot search regions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hdome = 31;
    wM = bwareaopen(wM,1000);
    dt = bwdist(~wM);
    dt = imfilter(dt,fspecial('gaussian',[31 31],11),'replicate');
    mx = imextendedmax(dt,hdome);
    mx = mx.*borderMask;
    mx = bwareaopen(mx,200);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set the moving and the fixed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    moving = single(modelDT);
    movingMask = single(modelMask);
    fixed = double(wM);
    fixedMask = [];


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set optimization parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lb = [.6 .6 -pi -50 -50];
    ub = [2.5 2.5 pi 50 50];
    b = [lb;ub];
    N = 100;
    optiPara.bounds = b;
    optiPara.swarmSize = N;

    dispImage = flattenMaskOverlay(bindVec(I),logical(mx));
    cotSeachRegions = mx;
    numberGroups = 7;
    maxCots = Inf;
    %maxCots = 1;
    templateData.image = moving;
    templateData.mask = movingMask;
    fixedData.image = fixed;
    fixedData.mask = fixedMask;
    % locate the cots via template matching
    [curveBank] = findCots(cotSeachRegions,numberGroups,maxCots,templateData,fixedData,dispImage,optiPara);
    
    % render the cot outlines and save the image
    oTemplate = moving;
    cropDisp = true;
    renderCurveBankToImage(I,curveBank,oTemplate,toPath);

    % crop the cots and make the masks
    [dataPoint] = cropCurveBank(fileName,curveBank,oTemplate,cropDisp);


    % extract the vein networks
    delay = false;
    [dataPoint] = extractVeinNetwork(dataPoint,oTemplate,fileName,delay);


    matFile = [toPath 'data.mat'];
    save(matFile,'dataPoint');


    % save single cot images
    renderSingleCotImage(dataPoint,toPath);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % render sheets to disk
    renderSheets = {};

    xlsFile = [toPath 'data.xlsx'];
    header = {'End_Point_X','End_Point_Y','Branch_Point_X','Branch_Point_Y'};
    summarySheetHeader = {'Number_of_Holes','Number_of_Branch_Points','Number_of_End_Points'};
    for cot = 1:numel(dataPoint)

        szEP = size(dataPoint(cot).veinNetwork.ePoints);
        szBP = size(dataPoint(cot).veinNetwork.bPoints);

        summarySheet{cot,1} = dataPoint(cot).numberHoles;
        summarySheet{cot,2} = num2str(szBP(1));
        summarySheet{cot,3} = num2str(szBP(1));

        
        H = max([szEP(1) szBP(1)]);
        oSheet = cell(H,4);
        %for h = 1:numel(header)
        %    oSheet{1,h} = header{h};
        %end
       
        c = 1;
        for i = 1:size(dataPoint(cot).veinNetwork.ePoints,1)
            oSheet{i,c} = num2str(dataPoint(cot).veinNetwork.ePoints(i,2));
            oSheet{i,c+1} = num2str(dataPoint(cot).veinNetwork.ePoints(i,1));
        end
        c = 3;
        for i = 1:size(dataPoint(cot).veinNetwork.bPoints,1)
            oSheet{i,c} = num2str(dataPoint(cot).veinNetwork.bPoints(i,2));
            oSheet{i,c+1} = num2str(dataPoint(cot).veinNetwork.bPoints(i,1));
        end
        
        renderSheets{cot+1} = cell2table(oSheet,'VariableNames',header);
        
    end
    renderSheets{1} = cell2table(summarySheet,'VariableNames',summarySheetHeader);
    

    
    for r = 1:numel(renderSheets)
        writetable(renderSheets{r},xlsFile,'FileType','spreadsheet','sheet',r);
    end
    

    
end
%% post loop
FilePath = '/mnt/spaldingdata/nate/makotoReturn/';
%FilePath = '/home/nate/Vein_pattern_images/OMY183_Col/';
mFileList = {};
FileExt = {'mat'};
mFileList = fdig(FilePath,mFileList,FileExt,1);
mdataPoint = [];
for e = 1:numel(mFileList)
    d = load(mFileList{e});

    %%%%%%%%%%%%%%%%%%
    % reextract network
    % extract the vein networks
    %%%%%%%%%%%%%%%%%%
    delay = false;
    [tmpDataPoint] = extractVeinNetwork(d.dataPoint,oTemplate,fileName,delay);


    mdataPoint = [mdataPoint;tmpDataPoint'];
    fprintf(['Done stacking:' num2str(e) ':' numel(mFileList) '\n']);
end
%%
[featureTable,imageStack,structData] = projectCotFeatureVector(mdataPoint);
%% added and checking region-branches
close all
for e = 1:numel(mFileList)
    try
        d = load(mFileList{e});
        for f = 1:numel(d.dataPoint)
            q = d.dataPoint(f);
            q.decorateVeins();
            q.veinNetwork.generateRegions(31);
            q.view();
            f
        end
        e
    catch
    end
end
%%
close all
qidx = find(featureTable.holeN == 2 & featureTable.endPointsN == 2);
qidx = find(featureTable.holeN == 2 & featureTable.branchPointsN <= 3 & featureTable.endPointsN <= 2);
dataSetPercent = numel(qidx)/numel(mdataPoint);

fprintf(['Selected ' num2str(100*dataSetPercent) ' percent of data.\n']);
mainPoints = [structData(qidx).mainBranchPoints]';
uI = mean(imageStack.holesStack(:,:,qidx),3);

[pMap] = pMapFromPoints(mainPoints,oTemplate);

imshow(uI,[]);
hold on
plot(mainPoints(:,2),mainPoints(:,1),'.')

%% crop curve bank
close all
per1 = .5;
per2 = 1;
perTotal = per1*per2;
perTotal = .5;
oTemplate = moving;




[d,nsz] = generateAffineImageDomain(oTemplate,2);
oTemplate = imresize(oTemplate,2);
for e = 1:numel(curveBank)

    fileName = cFileList{e};
    [pth,nm,ext] = fileparts(fileName);
    I = double(imread(fileName))/255;
    I = rgb2gray(I);

    for c = 1:size(curveBlock{e},1)
        trans0 = curveBlock{e}(c,1:5);
        trans = curveBlock{e}(c,6:10);
        trans0([1 2 4 5]) = trans0([1 2 4 5])*perTotal^-1;
        trans([1 2 4 5]) = trans([1 2 4 5])*perTotal^-1;
        
        fixedImageInMovingFrame = transformImage(oTemplate,I,d,trans,trans0,nsz);
        th = graythresh(fixedImageInMovingFrame);
        fixedMaskInMovingFrame = fixedImageInMovingFrame > th;

        modMask = imdilate(oTemplate > .8,strel('disk',5,0));
        modMask = modMask .* fixedMaskInMovingFrame;
        fixedMaskInMovingFrame = bwlarge(fixedMaskInMovingFrame);
        out = flattenMaskOverlay(bindVec(fixedImageInMovingFrame),fixedMaskInMovingFrame);
        out = flattenMaskOverlay(out,oTemplate > .8,.3,'b');
        out = flattenMaskOverlay(bindVec(fixedImageInMovingFrame),modMask>.8);
        
        dataPoint{e}(c).alignedImage = fixedImageInMovingFrame;
        dataPoint{e}(c).alignedMask = modMask;
        
        imshow(out,[]);
        drawnow
    end
end
%% quick view of aligned data
close all
for e = 1:numel(dataPoint)
    for c = 1:numel(dataPoint{e})
        out = flattenMaskOverlay(bindVec(dataPoint{e}(c).alignedImage),...
            logical(dataPoint{e}(c).alignedMask));
        out = flattenMaskOverlay(bindVec(dataPoint{e}(c).alignedImage),moving > .8);
        imshow(out,[])
        drawnow
        waitforbuttonpress
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% process and build model from aligned gray scale data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% value for spuring the skeleton
lengthFilter = 15;
% value to measure the cost for spurs that connect branch to end
snipAmount = 31;
% amount to dilate the skeleton
skeletonDilateAmount = 3;
% amount to erode the mask
totalMaskErodeAmount = 3;
% n-hood size for threshold
nsz = 71;
% sensetivity for threshold
sen = .45;
sen = .5; % This might work better
% small hole filter
smallHoleFilter = 200;
% gap close value
% this parameter uses morphological dilatation and as a result will has
% problem with V-shapes
gapClose = 3; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all
h1 = figure;
clear pointSets
cnt = 1;

for e = 1:numel(dataPoint)
    n = size(dataPoint{e});
    if ~all(n == 0)
        for m = 1:n(2)
            try
                n
                m

                imageToUse = dataPoint{e}(m).alignedImage.*dataPoint{e}(m).alignedMask;

                pointSets = singleCotFromImage(oTemplate,h1,imageToUse,dataPoint{e}(m).alignedMask,sen,gapClose,lengthFilter,smallHoleFilter,snipAmount,skeletonDilateAmount,totalMaskErodeAmount,nsz);
                
                f = fields(dataPoint{e}(m));
                for r = 1:numel(f)
                    dP(cnt).(f{r}) = dataPoint{e}(m).(f{r});
                end
                dP(cnt).pointSets = pointSets;
                dP(cnt).fileName = cFileList{e};
                cnt = cnt + 1;
                %waitforbuttonpress
            catch ME
                getReport(ME)
            end
        end
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build out prob maps for branch points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
branchPointPmask = zeros(size(modelTemplate));
for e = 1:numel(pointSets)
    idx = sub2ind(size(branchPointPmask),pointSets{e}.bPoints(:,1),pointSets{e}.bPoints(:,2));
    branchPointPmask(idx) = 1;
end
branchPointPmask = imdilate(branchPointPmask,strel('disk',11,0));
imshow(branchPointPmask,[]);
%%
    %{    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % write out single images
        for b = 1:numel(R)
            singleFile = [pth filesep nm '_' num2str(b) '_single.tif'];
            subI = imcrop(I,R(b).BoundingBox);
            imwrite(subI,singleFile);
        end




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % single cot analysis
        for b = 1:numel(R)


            singleCot(b,cFileList{e},I,M,cropBox,sen,spurAmount,snipAmount,skeletonDilateAmount,totalMaskErodeAmount,nsz)


        end

    end
end
%}



%{



    cp = size(BWmassTemplate)/2;
    wM = double(wM);
    alignPoints = [];
    massTemplate = imresize(wM,reSZ);
    %dwM = imerode(useM,strel('disk',21,0));
    dwM = useM;
    [f1,f2] = ndgrid(1:size(useM,1),1:size(useM,2));
    searchSpots = (mod(f1,filterSam) == 0) & (mod(f2,filterSam) == 0);
    searchSpots = searchSpots .* borderMask;


    pss = find(searchSpots(:));
    szQ = [max(sum(searchSpots,2)),max(sum(searchSpots,1))];
    toSearch = find(dwM(pss));

    %searchSpots = searchSpots & dwM;
    alignPoints = [];
    [alignPoints(:,1),alignPoints(:,2)] = find(searchSpots);
    searchMask = searchSpots > .8 & dwM > .8;
    imshow(searchMask,[]);
    
    tmp = find(searchSpots);
    
    [afidx] = find(searchMask(tmp) > .8);

    szQ = [max(sum(searchSpots,2)),max(sum(searchSpots,1))];
    K = zeros(size(alignPoints,1),1);
    N = 400;
    useM = double(useM);


%}