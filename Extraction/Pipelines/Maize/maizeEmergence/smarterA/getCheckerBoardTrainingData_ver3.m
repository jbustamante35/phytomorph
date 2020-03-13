%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct file list for training
% 1: get list of all files from user = jgustin
%       - order via sdig style
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 1: gather data from CyVerse
% this is the RILs only path
dataPath = '/iplant/home/jgustin/maizeData%'; 
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
%% step 2: parse data
FileList = {};
FileExt = {'tiff'};
fprintf(['Found:' num2str(numel(r)) ' records \n']);
for e = 1:numel(r)
    [p,nm,ext] = fileparts(r(e).DATA_NAME);
    if any(strcmp(ext(2:end),FileExt))
        FileList{end+1} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
    end
end
%% step 3: make stacks 
[FileList] = orderFrom_gdig(FileList,{});
%% step 4: sort stacks
rm = [];
for e = 1:numel(FileList)
    nm = [];
    n = [];
    try
        for s = 1:numel(FileList{e})
            [~,nm] = fileparts(FileList{e}{s});
            n(s) = str2num(nm);
        end

        [~,sidx]= sort(n);
        FileList{e} = FileList{e}(sidx);
    catch
        rm = [rm e];
    end
end
FileList(rm) = [];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step 2: build up color space maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numS = 10;
numI = 20;
setM = numel(FileList);
rS = randi(setM,1,numS);
per = .5;
nd = 3;
toSAM = 5*10000;
HISTO = zeros(toSAM*numS*numI,3);
str = 1;
stp = toSAM;
for s = 1:numel(rS)
    setM = numel(FileList{rS(s)});
    rI = randi(setM,1,numI);
    %{
    [tform] = getCameraPara(FileList{e}{1});
    CB = getRectifiedImage(FileList{e}{2},tform);
    %}
    for i = 1:numel(rI)
        try
            
            I = imread(FileList{rS(s)}{rI(i)});
            if per ~= 1
                I = imresize(I,per,'nearest');
            end
            I = permute(I,[3 1 2]);
            sz = size(I);
            I = reshape(I,[sz(1) prod(sz(2:3))])';
            I = I(randperm(size(I,1)),:);
            
            
            
            HISTO(str:stp,:) = I(1:toSAM,:);
            str = stp + 1;
            stp = str - 1 + toSAM;
            %{
            I = reshape(I,sz);
            I = ipermute(I,[3 1 2]);
            %}
            %{
            I = I + 1;
            lidx = sub2ind([256 256 256], I(:,1), I(:,2),I(:,3));
            for l = 1:numel(lidx)
                HISTO(lidx(l)) = HISTO(lidx(l)) + 1;
            end
            %}
        catch
            
        end
    end
    s
end
%%% clean up HISO
HISTO(str:end,:) = [];
%% 
[hgmm] = generateHGMMC(HISTO,[2 2 2]);
%%
HISTO_HSV = squeeze(rgb2hsv(permute(HISTO/255,[1 3 2])));
%%
HISTO_TOT = [HISTO/255 HISTO_HSV(:,2:3)];
%%
pdf = HISTO * sum(HISTO)^-1;
cdf = cumsum(pdf);
numD = 2000;
rS = randi(numel(cdf),numD,1);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CROP SAMPLE IMAGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[finalScore] = cropOutSamples_ver0(FileList{10});
%%
BKFileList = FileList;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read master list and intersect with file list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read master list measurements for emergence
ml = readtext('/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/EmergenceAssay_Master_list.csv');
toScan = {};
iPlantDirectory = 29;
for e = 2:size(ml,1)
    if ~isempty(ml{e,iPlantDirectory})
        if ischar(ml{e,iPlantDirectory})
            toScan{end+1} = ml{e,iPlantDirectory};
        else
            ml{e,iPlantDirectory}
            e
        end
    end
end
toScan = unique(toScan);
%% get intersection with masterlist data
CMD = 'iget -f /iplant/home/jgustin/maizeData/coleoptileEmergence/handScores/handscoreTrainSet.csv /home/nate/Downloads';
system(CMD)
D = readtext('/home/nate/Downloads/handscoreTrainSet.csv');
scoreNUMF = 250;
import java.util.Map.*;
import java.util.HashMap;
handData_t0 = HashMap();
handData_t1 = HashMap();
handData_t2 = HashMap();
SCORED = {};
for e = 2:size(D,1)
    emerFr = D{e,5};
    gidx = strfind(D{e,1},filesep);
    key = [lower(D{e,1}((gidx(end)+1):(end))) '-' lower(D{e,2}(1)) '-' lower(D{e,3}(1)) lower(num2str(D{e,4}))];  
    SCORED{end+1} = D{e,1};

    handData_t0.put(key,emerFr);
    eF = zeros(scoreNUMF,1);
    if ~isnan(emerFr) & ~isinf(emerFr) & emerFr <= 250
        eF(emerFr) = 1;
    end
    handData_t1.put(key,eF);
    handData_t2.put(key,cumsum(eF));
end
%% get intersection with masterlist data - OLD WAY
rootP = {};
for e = 1:numel(FileList)
    [rootP{e}] = fileparts(FileList{e}{1});
end
[U,sidx,~] = intersect(rootP,toScan);
FileList = FileList(sidx);
%{
%% load hand score
scoreNUMF = 250;
em = readtext('/mnt/snapper/nate/mirror_images/maizeData/jgustin/coleoptileEmergence/Emergence_hand_score_merged.csv');
iPlantDirectory = 29;
iplantName = 27+1;
iplantName = 29;
hcName = 20+1+1;
posName = 10+1+1;
genoName = 7+1+1;
FRAME_SCORE = 34;
SCORED = {};

import java.util.Map.*;
import java.util.HashMap;
handData_t0 = HashMap();
handData_t1 = HashMap();
handData_t2 = HashMap();
for e = 2:size(em,1)
    e
    gidx = strfind(em{e,iplantName},filesep);
    key = [lower(em{e,iplantName}((gidx(end)+1):(end))) '-' lower(em{e,hcName}(1)) '-' lower(em{e,29+1}(1)) lower(num2str(em{e,30+1}))];  
    SCORED{end+1} = em{e,iplantName};
    fidx = strfind(SCORED{end},filesep);
    key = lower([SCORED{end}((fidx(end)+1):end) '-' em{e,21+1}(1) '-' em{e,29+1} num2str(em{e,30+1})]);
    
    emerFr = em{e,end};
    %{
    if strcmp(lower(['20170220_Camera1-' LABELS{22}]),key)
        emerFr = 71;
    end
    if strcmp(lower(['20170220_Camera1-' LABELS{50}]),key)
        emerFr = 69;
    end
    if strcmp(lower(['20170220_Camera1-' LABELS{12}]),key)
        emerFr = 85;
    end
    if strcmp(lower(['20170220_Camera1-' LABELS{131}]),key)
        emerFr = 201;
    end
    
    
    if strcmp(lower(['20170220_Camera2-' LABELS{111}]),key)
        emerFr = 236;
    end
    if strcmp(lower(['20170220_Camera2-' LABELS{15}]),key)
        emerFr = 117;
    end
    if strcmp(lower(['20170220_Camera2-' LABELS{110}]),key)
        emerFr = 76;
    end
    if strcmp(lower(['20170220_Camera2-' LABELS{4}]),key)
        emerFr = 98;
    end
   %}
    
    
    if ischar(emerFr)
        emerFr = str2num(emerFr);
    end
    handData_t0.put(key,emerFr);
    eF = zeros(scoreNUMF,1);
    if ~isnan(emerFr) & ~isinf(emerFr) & emerFr <= 250
        eF(emerFr) = 1;
    end
    handData_t1.put(key,eF);
    handData_t2.put(key,cumsum(eF));
    
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% intersect the file list with  the scored list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UQ = unique(SCORED);
rootP = {};
for e = 1:numel(FileList)
    [rootP{e}] = fileparts(FileList{e}{1});
end
[U,sidx,~] = intersect(rootP,UQ);
FileList = FileList(sidx);
%% look for matches in file list- SEARCH
%toMatch= {'20170516_Camera3','20170608_Camera2','20170613_Camera4'};
toMatch1 = {'20170613_Camera2','20170613_Camera3'};
toMatch1 = {'20170131_Camera3'};
toMatch1 = {'20170613_Camera2'};
toMatch1 = {'20171121_Rack1_Camera1'};
toMatch1 = {'20171020_Rack1_Camera1'};
toMatch1 = {'20171013_Rack1_Camera3'};
kidx = zeros(size(FileList,1),1);
for e = 1:numel(FileList)
    for m = 1:numel(toMatch1)
        if ~isempty(strfind(FileList{e}{1},toMatch1{m}))
            kidx(e) = 1;
        end
    end
end
%FileList = FileList(find(kidx));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% extract to disk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
framesToMeasure = 250;
oPath = '/mnt/tetra/nate/forEM2';
for e = 1:2
       
        [tform] = getCameraPara(FileList{e}{1});
        
        
        CB = getRectifiedImage(FileList{e}{2},tform);
        
        %{
        I = imread(FileList{e}{1});
        imshow(I,[]);
        hold on
        for p = 1:size(imagePoints,1)
            text(imagePoints(p,1),imagePoints(p,2),num2str(p),'BackGround','r')
        end
        drawnow
        %}
        
        [boundingBox,centerPoints,MASK,I] = getCropBoxes(FileList{e},tform,168,10,150,1);
        [boundingBox,centerPoints,LABELS] = orderCropBoxes(MASK,boundingBox,centerPoints,I,true);
        cropTimeSeriesFromRectified(FileList{e},boundingBox,[1:1:framesToMeasure],LABELS,tform,oPath,[],[],101,[],[],[],[],[],[],[],[],[],[],[],[]);
%cropTimeSeriesFromRectified(FileList,boundingBox,timeSeries,LABELS,tform,storePath,emergenceNet,BESTemergence,nsz,zU,zE,zU2,zE2,WINDOW_SZ,Zmean,Zstd,CELLMASK,NUMPC,nI,subBOX,Tscore)
end
%% 
%% BIG NEXT
FilePath = '/mnt/tetra/nate/forEM2/';
tifFileList = {};
FileExt = {'tif','TIF'};
tifFileList = gdig(FilePath,tifFileList,FileExt,1);
%% split into ordered sets by key
imageKey = {};
for e = 1:numel(tifFileList)
    [p,n,ext] = fileparts(tifFileList{e});
    fidx = strfind(n,'--');
    imageKey{e} = n(1:(fidx(1)-1));
end
uniqueImageKey = unique(imageKey);
%{
%% check if any do not have 250 frames
for e = 1:numel(uniqueImageKey)
    parfor f = 1:numel(tifFileList)
        if ~isempty(strfind(tifFileList{f},uniqueImageKey{e}))
            grp(f) = e;
        end
    end
    e
end
nUQ = unique(grp);
rm = [];
for e = 1:numel(nUQ)
    if sum(grp==nUQ(e)) ~= 249
        rm = [rm e];
    end
end
%%
for e = 1:numel(rm)
    fidx = find(grp==rm(e));
    for f = 1:numel(fidx)
        delete(tifFileList{fidx(f)})
    end
end
%}
%% loader from disk
nsz = 151;



TYPE = zeros(numel(tifFileList),1);
TYPEBEST = TYPE;
VALUE = zeros(numel(tifFileList)/249,1);
cnt = 1;
cnt2 = 1;
cnt3 = 1;
CELLMASK = zeros(nsz,nsz);
CELLMASK((nsz-1)/2,(nsz-1)/2) = 1;
CELLMASK = bwdist(CELLMASK);
CELLMASK = double(CELLMASK < .9*(nsz-1)/2);
RC = regionprops(logical(CELLMASK),'BoundingBox');
BOX = RC(1).BoundingBox;
junk = imcrop(CELLMASK,BOX);
tiffLOAD = 1;
if tiffLOAD
    rZ = zeros([size(junk,1) size(junk,2) 3 numel(tifFileList)]);
end
TYPEF = [];
newM = imerode(CELLMASK,strel('disk',27,0));
vidx = find(newM);
for e = 1:numel(uniqueImageKey)
    
    baseName = uniqueImageKey{e};
    if tiffLOAD
        parfor fr = 1:249
            fileName = [FilePath baseName '--f' num2str(fr) '.tif'];
            Z{fr} = imresize(double(imread(fileName)),[nsz nsz]);
            Z{fr} = bsxfun(@times,Z{fr},CELLMASK);
            Z{fr} = imcrop(Z{fr},BOX);
            
            %nm
        end
    end
    
    
    YV = handData_t2.get(uniqueImageKey{e});
    BB = handData_t1.get(uniqueImageKey{e});
    
    BBC = imdilate(BB,strel('disk',5));
    BBC2 = imdilate(BB,strel('disk',9));
    %{
    R = regionprops(logical(BBC2-BBC),'PixelIdxList');
    BB = zeros(size(BBC));
    if ~isempty(R)

        BB1 = zeros(size(BBC));
        BB2 = zeros(size(BBC));
        BB1(R(1).PixelIdxList) = 1;
        if numel(R) == 2
            BB2(R(2).PixelIdxList) = 1;
        end
        BB = BB1 + BBC*2 + 3*BB2;
    end
    %}
    BB = BBC2;
    Yv = handData_t0.get(uniqueImageKey{e});
    for fr = 1:249
        TYPE(cnt2) = YV(fr);
        cnt2 = cnt2 + 1;
    end
    
    TYPEF = [TYPEF;YV'];
    
    for fr = 1:249
        TYPEBEST(cnt3) = BB(fr);
        cnt3 = cnt3 + 1;
    end
    
    
    
    %TYPE = [TYPE ; YV];
    VALUE(e) = Yv;
    
    if tiffLOAD
        for fr = 1:249
            rZ(:,:,:,cnt) = Z{fr};
            cnt = cnt + 1;
        end
    end
    
    
    e
end
%% make histograms for Value
parfor e = 1:size(rZ,4)
    tmp = rgb2hsv_fast(rZ(:,:,:,e)/255,'','V');
    [H(:,e)] = hist(tmp(:),linspace(0,1,256));
    e
end
%% select histogram from first stack first cell
uH = mean(H(:,1:249),2);
%% correct images
%crZ = rZ;
CELLMASK_sub = imcrop(CELLMASK,BOX);
cidx = find(CELLMASK_sub);


toMatch = rgb2hsv_fast(mean(rZ(:,:,:,1:10:end)/255,4),'','V');
  
  
parfor e = 1:size(rZ,4)
    tmp = rgb2hsv_fast(rZ(:,:,:,e)/255);
    
    
    tmp(:,:,3) = imhistmatch(tmp(:,:,3),toMatch);
    
    
    rZ(:,:,:,e) = 255*hsv2rgb(tmp);
    e
    % imshow(cat(2,rZ(:,:,:,e)/255,tmp),[]);
end
%% subtract first 20 frames
rZtmp = reshape(rZ,[size(rZ,1) size(rZ,2) 3 249 size(rZ,4)/249]);
for e = 1:size(rZtmp,5)
    tmpU = mean(rZtmp(:,:,:,1:20,e),4);
    rZtmp(:,:,:,:,e) = bsxfun(@minus,rZtmp(:,:,:,:,e),tmpU);
    e
end
rZtmp = reshape(rZtmp,size(rZ));
%% color segment
% init the arborist with rules for building trees
options = statset('Display','iter');
% creaete extract function to give to the arborist
filter = fspecial('gaussian',[21 21],5);
extractFunc = @(X)rgbAndhsvAndtextureExtract(X,filter);
% create tree growth rules for the arborist
suggestNumberOfClusterFunction = @(data,level,treePara)treePara(1)*(level<=treePara(2)) + 1*(level>treePara(2));
% create cluster parameter generating function for arborist
clusterFunctionFunctionGenerator = @(X,K)fitgmdist(X,K,'Options',options,'Start','plus','RegularizationValue',0.0001);
% feature selection function
idxSelectorFunction = @(X,L)logical([1 1 1 1 1 0]*(L < 2) + (L >= 2)*[0 0 0 0 0 1]);
% spec the tree parameters
maxBranch = 4;
maxDepth = 2;
% build the arborist
jA = arborist(suggestNumberOfClusterFunction,clusterFunctionFunctionGenerator,extractFunc,maxDepth,maxBranch,idxSelectorFunction);
%%
zCOLOR = jA.sampleElements(rZ,[3000 8000]);
%%
forest = jA.plantTrees(zCOLOR(1:10:size(zCOLOR,1),:),20);
%% view clusters over some examples
%ridx = randi(size(rZ,4),1,20);
for e = 1:100:size(rZ,4)
    oI = rZ(:,:,:,(e));
    [k p] = forest.clusterImage(oI);
    for i = 1:size(k,3)
        rgb = label2rgb(k(:,:,i));
        sig = cat(2,oI,rgb);
        imshow(sig,[]);
        title([num2str(i) '-' num2str(numel(unique(k(:,:,i))))])
        waitforbuttonpress
    end
end
%%
%% fold eigenvector into space 
toFold = zE;
for e = 1:size(zE,2)
    sz = (size(zE,1)/3).^.5;
    tmp = reshape(zE(:,e),[sz sz 3]);
    tmp = sum(tmp,3);
    foldedV(:,e) = tmp(:);
end
waveFunction2 = foldedV;
%% apply quantum framework
[X Y] = ndgrid(1:sz,1:sz);
OB{1} = sparse(diag(Y(:)));
OB{2} = sparse(diag(X(:)));
%%
[X Y] = ndgrid(1:sz(1),1:sz(1));
OB = cat(3,sparse(diag(X(:))),sparse(diag(Y(:))));
%%
zE = zE_0;
zU = zU_0;
%rZtmp = reshape(rZ,[size(rZ,1) size(rZ,2) 3 249 size(rZ,4)/249]);
close all
SIGGY = {};
h1 = figure;
h2 = figure;
h3 = figure;
[f1 f2] = ndgrid(linspace(-136/2,136/2,136),linspace(-136/2,136/2,136));
F = (f1.^2 + f2.^2).^.5;
MSK = F < 136/2;
FF = normpdf(F,0,136/6);
FF = MSK.*(FF);

CELLMASK_sub = imcrop(CELLMASK,BOX);
fidx = find(CELLMASK_sub);
str = 1;
stp = size(rZtmp,4) + str - 1;
L1 = zeros(136*136,size(rZ,4));
L2 = zeros(136*136,size(rZ,4));
for e = 1:size(rZtmp,5)
    fprintf(['Start: ' num2str(e) ':' num2str(size(rZtmp,5)) '\n'])
    tempy = [];

    %for tree = 2
    tree = 6;
    measureTrace = [];
    oTI = rZtmp(:,:,:,:,e);

    tmpL = zeros(136*136,size(rZtmp,4));
    
    
    parfor t = 1:size(rZtmp,4)
        tmpI = oTI(:,:,:,t);
        oI = tmpI;

        tmpI = tmpI(:);
        tmpC = PCA_REPROJ_T(tmpI,zE(:,1:3),zU);
        tmpI = PCA_BKPROJ_T(tmpC,zE(:,1:3),zU);
        modelImage = reshape(tmpI,size(oI));

        %[k,p] = forest.getTree(tree).clusterImage(oI);

        [fusionImage] = generateFusionImage(oI,modelImage,forest.getTree(tree),[2]);

        %ooI = oI;

        [k,p] = forest.getTree(tree).clusterImage(fusionImage);
        
        [ok,op] = forest.getTree(tree).clusterImage(oI);
        ok = ok;
        op1 = op(:,:,2);
        op2 = op(:,:,1);
        tmpL(:,t) = op1(:);
        tmpLK(:,t) = op2(:);
        
        
        waveFunction1 = [];
        for c = 1:size(p,3)
            tmp = p(:,:,c);
            tmp(:) = tmp / sum(tmp(:));
            waveFunction1(:,c) = tmp(:).^.5;
        end

        rgb = label2rgb(k);

        
        waveFunction2 = ones(size(waveFunction1,1),1);
        waveFunction2 = FF(:);
        waveFunction2 = waveFunction2 / norm(waveFunction2);
        [m] = qM(OB,waveFunction1,waveFunction2);

        [u edge] = measureUCG(m);

        measureTrace(:,t) = u';

        %{
        figure(h1);
        imshow([oI fusionImage]/255,[]);
        %imshow(rgb,[]);
        hold on
        plotUCG(m,edge)
        hold off
        title(num2str(e))
        drawnow

        figure(h2)
        [strain] = convertToStrain(measureTrace,11,10);
        plot(strain');
        hold on
        plot(TYPEF(e,1:t)*max(strain(:)),'k','LineWidth',3)
        hold off

        figure(h3)
        plot(measureTrace');
        hold on
        plot(TYPEF(e,1:t)*max(measureTrace(:)),'k','LineWidth',3)
        hold off
        %}
    end
    
    L1(:,str:stp) = tmpL;
    L2(:,str:stp) = tmpLK;
    
    str = stp + 1;
    stp = size(rZtmp,4) + str - 1;
    
    figure(h1);
    plot(measureTrace','r');
    hold on
    plot(TYPEF(e,:)*max(measureTrace(:)),'k','LineWidth',3);
    hold off
    figure(h2)
    [tsig] = convertToStrain(measureTrace,11,10);
    plot(tsig','b');
    hold on
    plot(TYPEF(e,:)*max(tsig(:)),'k','LineWidth',3);
    hold off

    drawnow
    SIGA{e} = measureTrace;

    fprintf(['End: ' num2str(e) ':' num2str(size(rZtmp,5)) '\n'])
end
%% THIS WAS USED __ NOW NEXT BLOCK BELOW
close all
uLL = mean(L1,2);
uLL2 = mean(L2,2);
uLLv = reshape(uLL,[136 136]);
uLLv2 = reshape(uLL2,[136 136]);
uLLv = CELLMASK_sub.*uLLv > .6;
uLLv2 = CELLMASK_sub.*uLLv2 > .7;

imshow(uLLv,[]);
figure;
imshow(uLLv2,[]);
drawnow
vidx = find(uLLv);
vidx2 = find(uLLv2);
str = 1;
stp = str + 3*numel(vidx) - 1;
VV = zeros(3*numel(vidx),size(rZ,4));
VV2 = zeros(3*numel(vidx2),size(rZ,4));
cnt = 1;
for e = 1:size(rZ,4)
    tmp = rZ(:,:,:,e);
    tmpS = [];
    tmpS2 = [];
    
    for k = 1:size(rZ,3)
        tempy = tmp(:,:,k);
        tmpS = [tmpS;tempy(vidx)];
        tmpS2 = [tmpS2;tempy(vidx2)];
    end
    VV(:,cnt) = tmpS(:);
    VV2(:,cnt) = tmpS2(:);
    cnt = cnt + 1;
    str = stp + 1;
    stp = str + 3*numel(vidx) - 1;
    e
end
%% NOW USED FOR VIDX
close all
newM = imerode(CELLMASK,strel('disk',27,0));
newM = imresize(CELLMASK,[size(rZ,1) size(rZ,2)]);
vidx = find(newM);
VV = zeros(3*numel(vidx),size(rZ,4));
cnt = 1;
for e = 1:size(rZ,4)
    tmp = rZ(:,:,:,e);
    tmpS = [];
    
    for k = 1:size(rZ,3)
        tempy = tmp(:,:,k);
        tmpS = [tmpS;tempy(vidx)];
    end
    VV(:,cnt) = tmpS(:);
    cnt = cnt + 1;
    e
end
%%
[U_VV,E_VV] = PCA_FIT_FULL_Tws(VV,10);
C_VV = PCA_REPROJ_T(VV,E_VV,U_VV);
%%
[U_VV2,E_VV2] = PCA_FIT_FULL_Tws(VV2,10);
C_VV2 = PCA_REPROJ_T(VV2,E_VV2,U_VV2);
%%
CC_V = reshape(C_VV,[size(C_VV,1) size(rZtmp,4) size(rZtmp,5)]);
CC_V2 = reshape(C_VV2,[size(C_VV2,1) size(rZtmp,4) size(rZtmp,5)]);
CC_zQ = reshape(zC,[size(zC,1) size(rZtmp,4) size(rZtmp,5)]);
%%

%% stack into training style
%% determine if popped NEEDED TO RUN FOR BSIG
close all
F_STACK = [];
bSIG = [];
for e = 1:size(TYPEF,1)
    tsig = CC_V(:,:,e);
    %tsig = SIGA{e};
    %tsig = sort(tsig,2);
    [tsig] = convertToStrain(tsig,11,10);
    F_STACK = cat(4,F_STACK,tsig);
    bSIG = [bSIG TYPEF(e,end)];
    %{
    plot(tsig','r');
    hold on
    plot(TYPEF(e,:),'k')
    axis([0 size(TYPEF,2) -7 7]);
    hold off
    drawnow
    %}
    fprintf(['Done with:' num2str(e) '\n']);
end
%%
kidx = kmeans(C_VV(5,:)',2);
kidx = reshape(kidx,[249 1008]);
%% 
close all
for e = 1%:size(TYPEF,1)
    plot(F_STACK(:,:,1,e)','r');
    hold on
    plot(TYPEF(e,:),'k')
    title(num2str(e))
    hold off
    drawnow
    pause(.3)
end
%%
close all
selE = 3;
simImage = zeros(size(rawImage));
tmp = reshape(E_VV(:,selE),[numel(vidx) 3]);
for k = 1:3
    tempy = simImage(:,:,k);
    tempy(vidx) = tmp(:,k);
    simImage(:,:,k) = tempy;
    simImage(:,:,k) = bindVec(simImage(:,:,k));
end
imshow(simImage,[]);
%% watch a movie
close all
toWatch = 727;
toWatch = fidx(midx);
VV_TMP = PCA_BKPROJ_T(CC_V(:,:,toWatch),E_VV,U_VV);
for e = 1:size(rZtmp,4)
    rawImage = rZtmp(:,:,:,e,toWatch)/255;
    simImage = zeros(size(rawImage));
    tmp = VV_TMP(:,e);
    tmp = reshape(tmp,[numel(vidx) 3]);
    for k = 1:3
        tempy = simImage(:,:,k);
        tempy(vidx) = tmp(:,k);
        simImage(:,:,k) = tempy;
    end
   
    
    
    imshow(cat(2,simImage/255,rawImage),[]);
    title(num2str(e));
    drawnow
    if TYPEF(toWatch,e) == 1
        
        waitforbuttonpress
    end
    

end
%%
lambda = myLDA(WTF,bSIG);
new = WTF*lambda;
nb = fitcnb(new, bSIG);
%%
WTF = [squeeze(mean(F_STACK(1:10,1,1,:),2))' squeeze(mean(F_STACK((end-11:end),1,1,:),2))'];
[Xloadings,Yloadings,Xscores,Yscores, beta,pctVar,mse,stats,Weights] = plsregress(WTF,bSIG',3);
new = [ones(size(WTF,1),1) WTF]*beta;
nb = fitcnb(new, bSIG);
YY = predict(nb,new);
sum(bSIG==YY')/numel(YY)
%%
YY = predict(nb,new);
%%
close all
fidx0 = find(bSIG==0);
fidx1 = find(bSIG==1);
sig0 = squeeze(F_STACK(:,:,1,fidx0));
sig1 = squeeze(F_STACK(:,:,1,fidx1));
%{
plot(sig1,'r');
hold on
plot(sig0,'k','LineWidth',3);
%}
figure;
plot(squeeze(mean(sig1,1)),'r')
hold on
plot(squeeze(mean(sig0,1)),'k')
%%

%%
WTF = squeeze(std(std(F_STACK,1,2),1,1))';
WTF2 = squeeze(mean(F_STACK,2))';
WTF = [WTF' WTF2];
%%
nb = fitcnb(WTF, bSIG);
%%
YY = predict(nb,WTF);
%%
Popped_layer = [imageInputLayer([size(F_STACK,1) size(F_STACK,2) 1]);

          convolution2dLayer([1 40],5,'Padding','same');
          reluLayer();
          maxPooling2dLayer([1 50],'Stride',2);


          fullyConnectedLayer(2);
          softmaxLayer();
          
          classificationLayer()];

options = trainingOptions('sgdm',...
    'MaxEpochs',1000,...
    'InitialLearnRate',.001,...
    'ExecutionEnvironment','parallel',...
    'Plots','training-progress');
CNN_net = trainNetwork(F_STACK,categorical(bSIG'),Popped_layer,options);
%%

%% train the when popped on only the popped
fidx = find(bSIG==1);
%fidx = 1:numel(bSIG);
YYY = [];
GG = [];
GG2 = [];
QQ = [];
SIGGGY  = {};
YYY = {};
YYYT = [];

for e = 1:numel(fidx)
    %{
    % bck to this
    tsig = CC_V(1:10,:,fidx(e));
    tsig2 = CC_V2(1:10,:,fidx(e));
    %tsig = SIGA{e};
    %tsig = sort(tsig,2);
    [tsig] = convertToStrain(tsig,11,10);
    [tsig2] = convertToStrain(tsig2,11,10);
    LLAM = [tsig',tsig2']';
    LLAM = [[tsig',tsig2']';lambda'*[tsig',tsig2']'];
    %}
    
    tsig = CC_V(1:10,:,fidx(e));
    [tsig] = convertToStrain(tsig,11,10);
    %{
    tsig = CC_zQ(1:5,:,fidx(e));
    tsig = convertToStrain(tsig,11,10);
    
    LLAM = tsig;
    %}
    LLAM = [tsig];
    %LLAM = [tsig;lambda'*tsig;lambda2'*tsig];
    %SIGGGY{e} = [tsig;tsig2];
    SIGGGY{e} = LLAM;
    
    
    %{
    newC1 = cluster(GMModel1,LLAM');
    newC2 = cluster(GMModel2,LLAM');
    newC = newC1.*(TYPEF(fidx(e),1:(end-1))' == 0) + (TYPEF(fidx(e),1:(end-1))' ~= 0).*(newC2+3);
    %}
    newC = TYPEF(fidx(e),1:(end-1))' == 1;
    
    
    
    YYY{e} = categorical(newC');
    YYYT = [YYYT;newC'];
    GG = [GG;YYY{e}'];
    tidx = find(double(YYY{e}')==2);
    TT = zeros(size(YYY{e}'));
    try
        TT(tidx(1)) = 1;
    catch
        TT(end) = [];
    end
    TT = imdilate(TT,strel('disk',5));
    %GG = [GG;YYY{e}'];
    GG2 = [GG2;TT];
    %QQ = [QQ;[tsig',tsig2']];
    %QQ = [QQ;[tsig']];
    QQ = [QQ;LLAM'];
end
nX = SIGGGY;
YY = YYYT;
%%
lambda = myLDA(QQ,GG);
lambda2 = myLDA(QQ,GG2);
%%
fidx = find(double(GG)==2);
GMModel = fitgmdist(QQ(fidx,:),2);
%%
GMModel1 = fitgmdist(QQ(double(GG)==1,:),3);
GMModel2 = fitgmdist(QQ(double(GG)==2,:),2);
%%

%%
%{
YYYT = categorical(YYYT);
YY = {};
for e = 1:size(YYYT,1)
    YY{e} = YYYT(e,:);
end
%}
% PLAY TRAIN
 options = trainingOptions('sgdm',...
            'InitialLearnRate',.01,...
            'MaxEpochs',1000,...
            'Shuffle','every-epoch',...
            'Verbose',true,...
            'ExecutionEnvironment','cpu',...
            'Plots','training-progress');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
layers = [ ...
    sequenceInputLayer(size(SIGGGY{1},1))
    lstmLayer(3)
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];

%trainedNet_FM = trainNetwork(SIGGGY,YY,layers,options);

trainedNet_FM = trainNetwork(nX,YY,layers,options);
%%
close all
SET = 50;
for e = 1:249
    img = rZtmp(:,:,:,e,SET)/255;
    imshow(img,[0 1])
    drawnow
end
%% view one of the corrected images
close all
imshow(rZ(:,:,:,15000)/255,[])
%% decompse loaded data - ALL
sz = size(rZ);
nZ = reshape(rZ,[prod(sz(1:3)) prod(sz(4))]);
[~,zC,zU,zE] = PCA_FIT_FULL_T(nZ,25);
%%
FULL_LABEL = TYPEF(:,1:(end-1))';
FULL_LABEL = FULL_LABEL(:);
[zU_0,zE_0] = PCA_FIT_FULL_Tws(nZ,25);
zC = PCA_REPROJ_T(nZ,zE_0,zU_0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% tangle up the game
%% split e-vectors into [T;N]
zUn = zU / norm(zU);
for b = 1:size(zE,2)
    zET(:,b) = zUn*(zE(:,b)'*zUn);
    zEN(:,b) = zE(:,b) - zET(:,b);
    zET(:,b) = zET(:,b) / norm(zET(:,b));
    zEN(:,b) = zEN(:,b) / norm(zEN(:,b));
end
%% view the game
close all
figure;
vU = reshape(zU,sz(1:3));
imshow(vU/255,[]);
vS = std(zC,1,2);
vTAN = [];
for e = 1:5

    sig = (reshape(zET(:,e).^2,sz(1:3)));
    msig = sum(sig,3);
    for k = 1:size(sig,3)
        sig(:,:,k) = bindVec(sig(:,:,k));
    end
    
    vTAN = [vTAN msig];
end
figure;
imshow(vTAN,[]);
title('TAN');
vNOR = [];
for e = 1:5
    sig = (reshape(zEN(:,e).^2,sz(1:3)));
    msig = sum(sig,3);
    for k = 1:size(sig,3)
        sig(:,:,k) = bindVec(sig(:,:,k));
    end
    vNOR = [vNOR msig];
end
figure;
imshow(vNOR,[]);
title('NOR');
vTOT = [];
for e = 1:5
    vTOT = [vTOT sum(reshape(zE(:,e).^2,sz(1:3)),3)];
end
figure;
imshow(vTOT,[]);
title('TOT');


%% try smallest
[smU,smE] = PCA_FIT_FULL_Tws(nZ,3,'smallestreal');

%% make model - start with gathering red cap color
uCell = mean(rZ,4);
cellMask = zeros(size(uCell,1),size(uCell,2));
cellMask(end/2,end/2) = 1;
cellMask = bwdist(cellMask);
cellMask = cellMask < 30;
cellMaskpt = find(cellMask);
redColor = [];
for e = 1:5:5000
    tmpH = [];
    for k = 1:3
        tmp = rZ(:,:,k,e);
        tmpH = [tmpH tmp(cellMaskpt)];
    end
    redColor = [redColor ; tmpH];
    e
end
UredCapColor = mean(redColor);
SredCapColor = cov(redColor);
%% get pvc color
pvcColor = [];
for e = 1:5:5000
    tmp = rZ(:,:,:,e);
    borderMask = tmp(:,:,1) ~= 0;
    sz = size(tmp);
    tmp = reshape(tmp,[prod(sz(1:2)) prod(3)]);
    tmp = mvnpdf(tmp,UredCapColor,SredCapColor);
    tmp = reshape(tmp,sz(1:2));
    tmp = bindVec(tmp);
    tmpMask = tmp > .000000001;
    tmpMask = imfill(tmpMask,'holes');
    tmpMask = bwlarge(tmpMask);


    DIST = bwdist(tmpMask);


    ringMask = DIST < 9 & DIST > 4 & borderMask;
    
    out = flattenMaskOverlay(rZ(:,:,:,e)/255,ringMask);
    ridx = find(ringMask);
    tmpH = [];
    for k = 1:3
        tmpK = rZ(:,:,k,e);
        tmpH = [tmpH tmpK(ridx)];
    end
    pvcColor = [pvcColor ; tmpH];

    R = regionprops(logical(tmpMask),'Centroid','MajorAxis','MinorAxis');
    imshow(out,[]);
drawnow
    e
end

upvcColor = mean(pvcColor);
spvcColor = cov(pvcColor);
%%
obj = fitgmdist([redColor;pvcColor],2);
%%
fullLabel = fitgmdist(fullSegment,2,'Start','plus','Options',options);
%%
close all
for e = 29900:35000
    [idx,model,out] = getCellModel2(rZ(:,:,:,e),fullLabel,fullLabel2,1,false);
    RGB = label2rgb(idx);
    imshow(cat(2,rZ(:,:,:,e)/255,out,double(RGB)/255),[]);
    drawnow
end
%% try next color segment
 tmpSTACK = [];
for e = 1:1000
    baseName = uniqueImageKey{e};
    for fr = 200
        fileName = [FilePath baseName '--f' num2str(fr) '.tif'];
        tmpI = imresize(double(imread(fileName)),[nsz nsz]);
        tmpSTACK = cat(4,tmpSTACK,tmpI);
        %ftmpI = imfilter(tmpI,fspecial('gaussian',[21 21],2));
        ftmpI = tmpI;
        [idx,model,out] = getCellModel2(ftmpI,fullLabel,fullLabel2,1,false);
        RGB = label2rgb(idx);
        imshow(cat(2,tmpI/255,out,double(RGB)/255),[]);
        drawnow
    end
    e
end
%% extract the cap only distribution
for e = 1:10
    tmpI = rZ(:,:,:,e);
    [idx,model,POS,out] = getCellModel2(ftmpI,fullLabel,fullLabel2,1,spaceFunc,false);
end
%% try next color segment
% I THINK HERE IS WHERE I START 
%rZtmp = reshape(rZ,[size(rZ,1) size(rZ,2) 3 249 size(rZ,4)/249]);
TOT = [];
spaceFunc{1} = @(X)(X.*bwlarge(X > .0001));
spaceFunc{2} = @(X)X;
spaceFunc{3} = @(X)X;
close all
h1 = figure;
h2 = figure;
%STORE = [];
disp = true;

%[yi,xi] = ksdensity(TOT(TOT < 100),linspace(0,100,1000));
for e = 1:1008%1:size(rZtmp,5)
    SIG = [];
    L = [];
    for t = 1:1:size(rZtmp,4)
        tmpI = rZtmp(:,:,:,t,e);
        %{
        baseName = uniqueImageKey{e};
        fileName = [FilePath baseName '--f' num2str(t) '.tif'];
        tmpI = imresize(double(imread(fileName)),[nsz nsz]);
        %}
        
        % project into model space
        numC = 5;
        sz = size(tmpI);
        JJ = reshape(tmpI,[prod(sz) 1]);
        JC = PCA_REPROJ_T(JJ,zE(:,1:numC),zU);
        JC = PCA_BKPROJ_T(JC,zE(:,1:numC),zU);
        tmpI = reshape(JC,sz);
        % make disk mask
        zMASK = tmpI(:,:,1) ~= 0;

        %ftmpI = imfilter(tmpI,fspecial('gaussian',[21 21],2));
        ftmpI = tmpI;
        [idx,model,POS,out] = getCellModel2(ftmpI,fullLabel,fullLabel2,1,spaceFunc,false);
        % copy over raw data into model for fusion
        repidx = find(bwlarge(idx==1));
        for k = 1:3
            tmp = rZtmp(:,:,k,t,e);
            tmpTar = tmpI(:,:,k);
            tmpTar(repidx) = tmp(repidx);
            tmpI(:,:,k) = tmpTar;
        end
        [mea,POSRED] = splitEnergyLevels(tmpI,idx,[yi;xi],bwlarge(idx==1));

        
%{
        TOT = [TOT;mea(idx==1)];
%}
        
        
        %{
        SIG(1,t) = sum(idx(:)==1);
        SIG(2,t) = sum(idx(:)==2);
        SIG(3,t) = sum(idx(:)==3);
        
        if t == 1
            last = idx;
        else
            cnt = 1;
            for s1 = 1:3
                for s2 = 1:3
                    if s1 ~= s2
                        %SIG(cnt+3,t) = sum(idx(:) == s1 & last(:) == s2);
                        cnt = cnt + 1;
                    end
                end
           end
        end
       
        %}


        rPOS = POS;
        tL = [];
        for k = 1:3
            tL(k) = norm(rPOS(:,1)-rPOS(:,2));
            %plot(rPOS(2,1:2),rPOS(1,1:2),'k');
            rPOS = circshift(rPOS,[0 1]);
        end
        tL = [tL' ; norm(POS(:,1) - POSRED')];
        L(:,t) = tL;

        %L(4,t) = norm(POS(:,1) - POSRED');
        
         CL = {'r' 'c' 'y'};
        RGB = label2rgb(idx);
        if disp
            figure(h1);
            %imshow(cat(2,tmpI/255,out,double(RGB)/255),[]);
            imshow(tmpI/255,[]);
            hold on
            for k = 1:size(POS,2)
                plot(POS(2,k),POS(1,k),[CL{k} 'o']);
                hold on
            end
            plot(POSRED(1,2),POSRED(1,1),'r*');

            rPOS = POS;
            for k = 1:3
                plot(rPOS(2,1:2),rPOS(1,1:2),'k');
                rPOS = circshift(rPOS,[0 1]);
            end
            plot([POS(2,1) POSRED(2)],[POS(1,1) POSRED(1)],'w');
            plot(size(RGB,2)/2,size(RGB,1)/2,'k.');
            hold off
            title(num2str(e))
            sL = [];
            for s = 1:size(L,1)
                sL(s,:) = imfilter(L(s,:),fspecial('average',[1 10]),'replicate');
            end
            TT = 10;

            initL = mean(sL(:,1:(min(t,TT))),2);
            STR = bsxfun(@minus,sL,initL);
            STR = bsxfun(@minus,STR,initL.^-1);
            STR(end+1,:) = mean(abs(STR),1);
            figure(h2);
            plot(STR');
            hold on
            plot(TYPEF(e,1:t)*5,'k')
            hold off
            drawnow
            e
        end
    end
    %e
    STORE(:,:,e) = L;
    %waitforbuttonpress
end
%%

%% label mean image
[idx,model,out] = getCellModel2(mean(tmpSTACK,4),fullLabel,fullLabel2,1,false);
%% try next color segment - THIS SEGMENT NOT THE RED CAP PLUS PVC
fullSegment = [];
for e = 1:1000
    baseName = uniqueImageKey{e};
    for fr = 200
        fileName = [FilePath baseName '--f' num2str(fr) '.tif'];
        tmpI = imresize(double(imread(fileName)),[nsz nsz]);
        %tmpI = imresize(tmpI,.45);
        tmpSZ = size(tmpI);
        tmpI = reshape(tmpI,[prod(tmpSZ(1:2)) tmpSZ(3)]);
        fullSegment = [fullSegment;tmpI];
    end
    e
end
%% fit distributions
options = statset('Display','Iter');
fullLabel = fitgmdist(fullSegment,2,'Start','plus','Options',options);
capIDX = 1;
[trainIdx] = cluster(fullLabel,fullSegment);
fidx = trainIdx ~= capIDX;
fullLabel2 = fitgmdist(fullSegment(fidx,:),2,'Start','plus','Options',options);
%% gather measure for cap only
JJ = [];
for e = 1:1000
    for fr = 1:50:200
        tmpI = rZtmp(:,:,:,t,e);
        sz = size(tmpI);
        g = reshape(tmpI,[prod(sz(1:2)) sz(3)]);
        tmpI = tmpI(:,:,1);
        [d1 d2] = gradient(tmpI);
        d = (d1.^2 + d2.^2).^5;
        [gidx] = cluster(fullLabel,g);
        gidx = reshape(gidx,sz(1:2));
        gidx = bwlarge(gidx==capIDX);
        gidx = imerode(gidx,strel('disk',21,0));
        d = d.*gidx;
        JJ = [JJ;d(find(gidx))];
    end
    e
end
%%
%ksdensity(JJ)
close all
hJJ = JJ(1:5:end);
[yi,xi] = ksdensity(hJJ(hJJ < 1.8),linspace(0,1.8,30000));
plot(xi,yi)
%% model2
[idx,model,out] = getCellModel2(tmpI,fullLabel,fullLabel2,1,false);
%%
fileName = [FilePath baseName '--f' num2str(fr) '.tif'];
tmpI = imresize(double(imread(fileName)),[nsz nsz]);

[tform] = getCameraPara(FileList{91}{1});
tmpI = 255*getRectifiedImage(FileList{91}{2},tform,[]);
%tmpI = imread(FileList{91}{2});
tmpI = double(tmpI);
tmpI = imfilter(tmpI,fspecial('gaussian',[21 21],7));
[idx,model,out] = getCellModel(tmpI,fullLabel,false);
imshow(idx==1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reshape loaded data
%zCfinal = zC;
SELDIMS = 1:5;
[zCfinal Zmean Zstd] = zscore(zC(SELDIMS,:),1,2);
%zCfinal = [zCfinal;zC];
Hsz = size(zCfinal);
zCfinal = reshape(zCfinal,[Hsz(1) 249 size(zCfinal,2)/249]);
zCfinal = reshape(zCfinal,[size(zCfinal,1) size(zCfinal,2) 1 size(zCfinal,3)]);


rCfinal = reshape(zC(SELDIMS,:),[Hsz(1) 249 size(zC,2)/249]);
rCfinal = reshape(rCfinal,[size(rCfinal,1) size(rCfinal,2) 1 size(rCfinal,3)]);
%% setup rClass
BEST = [];
WINDOW_SZ = 11;
str = 1;
HUM = [];
for tr = 1:size(zCfinal,4)
    stp = str + 249-1;
    tmpS = im2col(zCfinal(:,:,1,tr),[size(zCfinal,1) WINDOW_SZ]);
    tmpS = reshape(tmpS,[size(zCfinal,1) WINDOW_SZ size(tmpS,2)]);
    
    tmpSz = im2col(rCfinal(:,:,1,tr),[size(rCfinal,1) WINDOW_SZ]);
    tmpSz = reshape(tmpSz,[size(rCfinal,1) WINDOW_SZ size(tmpSz,2)]);
    
    Ze = zeros(size(tmpS));
    
    tmpS = reshape(tmpS,[size(tmpS,1) size(tmpS,2) 1 size(tmpS,3)]);
    tmpSz = reshape(tmpSz,[size(tmpSz,1) size(tmpSz,2) 1 size(tmpSz,3)]);
    Ze = reshape(Ze,[size(Ze,1) size(Ze,2) 1 size(Ze,3)]);
    GG = cat(3,tmpS,tmpSz,Ze);
    
    BEST = cat(4,BEST,GG);
    
    tmpBB = TYPEBEST(str:stp);
    tmpBB = tmpBB((WINDOW_SZ-1)/2:(end-(WINDOW_SZ-1)/2-1));
    
    HUM = [HUM;tmpBB];
    tr
    str = stp + 1;
end
%BEST = reshape(BEST,[size(BEST,1) size(BEST,2) 1 size(BEST,3)]);
%% train dynamic
layers = [ ...
    imageInputLayer([size(BEST,1) size(BEST,2) 3])
    convolution2dLayer([1 11],15)
    reluLayer()
    %maxPooling2dLayer([1 3],'Stride',1);
    %convolution2dLayer([1 3],15)
    %reluLayer()
    convolution2dLayer([size(BEST,1) 1],5)
    reluLayer()
    %maxPooling2dLayer([3 1])
    fullyConnectedLayer(2)
    softmaxLayer()
    classificationLayer()];
options = trainingOptions('sgdm','MaxEpochs',8,'InitialLearnRate',0.00005,'ExecutionEnvironment','parallel','Plots','training-progress');
emergenceDynamic = trainNetwork(single(BEST),categorical(HUM),layers,options);
%%

%% my metrics implement - frenet 
K = [];
for e = 1:size(rCfinal,4)
    K(:,:,:,e) = myFrenet(rCfinal(:,:,1,e),[7],@(X)funcG(X));
    e
end
K = abs(K);
ksz = size(K);
K = reshape(K,[ksz(1) prod(ksz(2:4))]);
[nK,Kmu,Ksigma] = zscore(K,1,2);
nK = reshape(nK,ksz);
K = reshape(K,ksz);
%% LSTM
close all
toA = 20;
for e = 1:size(zCfinal,4)
    % stack both QR and zscore QR
    nX{e} = [nK(:,:,:,e)];
end
for e = 1:(size(TYPEF,1))
    YY{e} = categorical(TYPEF(e,2:end));
end

%% LSTM - quantum
close all
toA = 20;
majorStackX = [];
for e = 1:size(STORE,3)
    tmp = STORE(:,:,e);
    tmp = imfilter(tmp,fspecial('average',[1 11]),'replicate');
    init = mean(tmp(:,1:10),2);
    tmp = bsxfun(@minus,tmp,init);
    tmp = bsxfun(@times,tmp,init.^-1);
   
    st = tmp;
    tmp = gradient(tmp);
    tmp = imfilter(tmp,fspecial('average',[1 11]),'replicate');
    tmp = cat(1,tmp,st);
    
    %tmp = st;
    % stack both QR and zscore QR
    nX{e} = tmp;
    majorStackX = cat(3,majorStackX,tmp);
    e
end
majorStackY = [];
for e = 1:(size(TYPEF,1))
    YY{e} = categorical(TYPEF(e,2:end));
    newY = diff(TYPEF(e,1:end));
    newY = imdilate(newY,strel('disk',3));
    majorStackY = [majorStackY;newY];
    e
end
%% linear fun
sigX = [];
sigY = [];
sigIDX = [];
for e = 1:size(majorStackX,3)
    TsigX = im2col(majorStackX(:,:,e),[size(majorStackX,1) 11],'sliding');
    TsigY = im2col(majorStackY(e,:),[1 11]);
    sigX = [sigX;TsigX'];
    sigY = [sigY;TsigY((end-1)/2,:)'];
    sigIDX = [sigIDX;e*ones(size(TsigY,2),1)];
    e
end
%%
close all
u111 = mean(sigX,1);
nsigX = bsxfun(@minus,sigX,u111);

%{
[Xloadings,Yloadings,Xscores,Yscores,beta,pctVar,mse,stats,lambda1] = plsregress(nsigX,sigY,5);
lambda2 = mynLDA(nsigX*lambda,sigY,1,1);
funnySig = nsigX*lambda*lambda2;
%}


fidx = find(sigY==1);
[qS qC qU qE qL qERR qLAM] = PCA_FIT_FULL(sigX(fidx,:),2);
qC = PCA_REPROJ(sigX,qE,qU);
qA = PCA_BKPROJ(qC,qE,qU);
qERR = sum((qA - sigX).^2,2).^.5;
hope = [qC qERR];
funnySig = hope;
%{
grp = kmeans(nsigX(fidx,:),3);
GG = sigY;
GG(fidx) = grp;
%}
%GG = sigY;
%Mdl = fitcnb(funnySig,GG);
%sigYPRE = predict(Mdl,funnySig);

nsigX = funnySig;



u111 = mean(sigX,1);
nsigX = bsxfun(@minus,sigX,u111);


[Xloadings,Yloadings,Xscores,Yscores,beta,pctVar,mse,stats,lambda1] = plsregress(nsigX,sigY,3);
lambda2 = mynLDA(nsigX*lambda1,sigY,1,3);
funnySig = nsigX*lambda1*lambda2;

Mdl = fitcnb(funnySig,sigY);
sigYPRE = predict(Mdl,funnySig);





UQ = unique(sigIDX);
mag = .1*max(funnySig(:));

h1 = figure;
h2 = figure;
for u = 1:numel(UQ)
    fidx = sigIDX == u;
    TsigX = funnySig(fidx,:);
    TsigY = sigY(fidx,:);
    TsigYPRE = sigYPRE(fidx,:);
    
    
    figure(h1);
    plot(TsigX);
    hold on
    plot(TsigY*.2*mag);
    plot(TsigYPRE*.1*mag,'g');
    hold off
    axis([0 numel(TsigY) -mag mag]);
    drawnow
    pause(.1)
    
    
    figure(h2);
    plot(sigX(fidx,:));
    hold on
    plot(TsigY*.2);
    hold off
    waitforbuttonpress
end

%% view signal
close all
mag = max(nK(:));
mag = 0;
for e = 1:numel(nX)
    tmpSig = nX{e}';
    plot(tmpSig);
    hold on
   
    mag = max([mag;tmpSig(:)]);
     axis([0 250 -mag mag])
    plot(mag*(double(YY{e})-1),'r')
    hold off
    drawnow
    pause(.2)
    
end
%% i think this will train local - i dont think i need this
%emergenceDynamic = trainNetwork(nX,YY,layers,options);

%% I dont need this - i think - now on GPU
%{
optimVars = [
    optimizableVariable('lstmStates',[2 20],'Type','integer')
    optimizableVariable('InitialLearnRate',[1e-3 5e-2],'Transform','log')
    optimizableVariable('Momentum',[0.8 0.95])
    optimizableVariable('L2Regularization',[1e-10 1e-2],'Transform','log')];
maxE = 10;
toSlow = 1;
exeEnvironment = {'cpu',true};


optVarsTest.lstmStates = 2;
optVarsTest.InitialLearnRate = 1e-3;
optVarsTest.Momentum = .8;
optVarsTest.L2Regularization = 1e-10;



SLAP = 750;
[valError,cons] = makeObjFcn_emerge(...
                                    optVarsTest,...
                                    nX(1:SLAP),...
                                    YY(1:SLAP),...
                                    nX((SLAP+1):end),...
                                    YY((SLAP+1):end),...
                                    toSlow,...
                                    'none',...
                                    maxE,...
                                    'cpu');

maxE = 400;
toSlow = 1;
exeEnvironment = {'gpu',false};
[BayesObject] = hyperPdeploy_emerge(optimVars,nX(1:SLAP),...
                           YY(1:SLAP),...
                            nX((SLAP+1):end),...
                            YY((SLAP+1):end),...
                            exeEnvironment,maxE,toSlow);
%}
%% remote GPU - train on GPU via condor and hyper parameters

% check inputs
rm = [];
for e = 1:numel(nX)
    if ~(size(nX{e},2) == size(YY{e},2))
        nX{e} = nX{e}(:,1:size(YY{e},2));
        e
    end
    if (numel(unique(YY{e})) ~= 2)
        rm = [rm e];
        %newY = (double(YY{e})-1) ~= 0;
        %YY{e} = categorical(newY);
        e
    end
end
nX(rm) = [];
YY(rm) = [];
SLAP = 750;
SLAP = 600;
func = cFlow('hyperPdeploy_emerge');
func.setMCRversion('v930');
func.setMemory('8000');
func.setGPU(1);
% max function evaluations
maxE = 200;
% slow down training to get stable error
toSlow = 1;
% 
maxEval = 10000;
% max time for training hyper parameters
maxTime = 2*60*60;
% execution environment gpu
exeEnvironment = {'gpu',false};
% setup the variables to optimize
optimVars = [
    optimizableVariable('lstmStates',[2 10],'Type','integer')
    optimizableVariable('InitialLearnRate',[1e-3 5e-2],'Transform','log')
    optimizableVariable('Momentum',[0.8 0.95])
    optimizableVariable('L2Regularization',[1e-10 1e-2],'Transform','log')];
% call the func and store the result in beta0 = b0
b_full2 = func(optimVars,nX(1:SLAP),...
                           YY(1:SLAP),...
                            nX((SLAP+1):end),...
                            YY((SLAP+1):end),...
                            exeEnvironment,maxE,toSlow,maxE,maxTime,size(nX{1},1));

auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
func.submitDag(auth,50,50);
%% GPU emerge
%{
SLAP = 750;
SLAP = 600;
func = cFlow('hyperPdeploy_emerge');
func.setMCRversion('v930');
func.setMemory('8000');
func.setGPU(1);
% max function evaluations
maxE = 200;
% slow down training to get stable error
toSlow = 1;
% 
maxEval = 10000;
% max time for training hyper parameters
maxTime = 2*60*60;
% execution environment gpu
exeEnvironment = {'gpu',false};
% setup the variables to optimize
optimVars = [
    optimizableVariable('lstmStates',[2 10],'Type','integer')
    optimizableVariable('InitialLearnRate',[1e-3 5e-2],'Transform','log')
    optimizableVariable('Momentum',[0.8 0.95])
    optimizableVariable('L2Regularization',[1e-10 1e-2],'Transform','log')];
% call the func and store the result in beta0 = b0
b_fullEM2 = func(optimVars,eX(1:SLAP),...
                           eY(1:SLAP),...
                            eX((SLAP+1):end),...
                            eY((SLAP+1):end),...
                            exeEnvironment,maxE,toSlow,maxE,maxTime,size(eX{1},1),'last');

auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
func.submitDag(auth,50,50);
%}
%% LOCAL EMERGE
hParaE2 = cFlowLoader(b_fullEM2);
%hPara.XAtMinObjective = hPara.XAtMinObjective*100;
[valErrorE3,consE3,emergeNET3] = makeObjFcn_emerge(...
                                    hParaE2.XAtMinEstimatedObjective,...
                                    eX,...
                                    eY,...
                                    '',...
                                    '',...
                                    0,...
                                    'training-progress',...
                                    2000,...
                                    'cpu',...
                                    size(eX{1},1),...
                                    'last');
%% regrade
FilePath = '/mnt/tetra/nate/RED/';
redMATFileList = {};
FileExt = {'mat'};
redMATFileList = gdig(FilePath,redMATFileList,FileExt,1);
redTHRESH = linspace(.51,.8,20);
JfinalScore = {};

redTHRESH = linspace(.51,.8,10);
JfinalScore = {};
oPath = '/mnt/tetra/nate/forJEFFredCapps4/'
mkdir(oPath)
for t = 1:numel(redTHRESH)
    tPath = [oPath 'threshold' num2str(t) filesep];
    mkdir(tPath)
    parfor e = 1:numel(redMATFileList)
        tic
        [pth,nm,ext] = fileparts(redMATFileList{e});
        tName = [tPath nm '.csv'];
        [JfinalScore{t,e},~,~] = gradeRedCaps(redMATFileList{e},emergeNET3,trainedNet_FM,redTHRESH(t),.5);
        csvwrite(tName,JfinalScore{t,e});
        toc
    end
end
%% query camera
qS = ['20170131_Camera4'];
qS = ['20170317_Camera2'];
qS = ['20170705_Camera1'];
qS = ['20171013_Rack1_Camera3'];
for e = 1:numel(redMATFileList)
    if ~isempty(strfind(redMATFileList{e},qS))
        fnd = e;
    end
end
[JfinalScore{t,e},~,~] = gradeRedCaps(redMATFileList{fnd},emergeNET3,trainedNet_FM,redTHRESH(t),.5);
%% test with guosheng data
matFile = '/mnt/tetra/nate/phytoMorph_Image_Phenomics_Tool_Kit_analysis1-2018-02-01-16-19-20.3/output/20180104Camera3_dataPackage.mat';
matFile = '/mnt/tetra/nate/phytoMorph_Image_Phenomics_Tool_Kit_analysis1-2018-02-08-18-43-20.2/output/20180104Camera3_dataPackage.mat';
gfs = gradeRedCaps(matFile,emergeNET3,trainedNet_FM,.5,.5);
csvwrite('/mnt/tetra/nate/phytoMorph_Image_Phenomics_Tool_Kit_analysis1-2018-02-01-16-19-20.3/output/20180104Camera3_new.csv',gfs);
%% attempt the training session on the local cluster - same as above but local
exeEnvironment = {'cpu',true};
optimVars = [
    optimizableVariable('lstmStates',[2 20],'Type','integer')
    optimizableVariable('InitialLearnRate',[1e-3 5e-2],'Transform','log')
    optimizableVariable('Momentum',[0.8 0.95])
    optimizableVariable('L2Regularization',[1e-10 1e-2],'Transform','log')];

bO_local = hyperPdeploy_emerge(optimVars,nX(1:SLAP),...
                           YY(1:SLAP),...
                            nX((SLAP+1):end),...
                            YY((SLAP+1):end),...
                            exeEnvironment,maxE,toSlow,maxE,maxTime,size(nX{1},1));

%% perform local training with hyper parameters from condor
%{
    YYYT = categorical(YYYT);
    YY = {};
    for e = 1:size(YYYT,1)
        YY{e} = YYYT(e,:);
    end
%}
% check inputs
rm = [];
for e = 1:numel(nX)
    if ~(size(nX{e},2) == size(YY{e},2))
        nX{e} = nX{e}(:,1:size(YY{e},2));
        e
    end
    if (numel(unique(YY{e})) ~= 2)
        rm = [rm e];
        %newY = (double(YY{e})-1) ~= 0;
        %YY{e} = categorical(newY);
        e
    end
end
nX(rm) = [];
YY(rm) = [];
hPara = cFlowLoader(b_full2);
%hPara.XAtMinObjective = hPara.XAtMinObjective*100;
[valError,cons,trainedNet_FM] = makeObjFcn_emerge(...
                                    hPara.XAtMinEstimatedObjective,...
                                    nX,...
                                    YY,...
                                    '',...
                                    '',...
                                    0,...
                                    'training-progress',...
                                    2000,...
                                    'cpu',...
                                    size(nX{1},1),...
                                    'sequence');
                                %%
%hPara.XAtMinObjective = hPara.XAtMinObjective*100;
[valErrorER,consER,emergeNET] = makeObjFcn_emerge(...
                                    hPara.XAtMinEstimatedObjective,...
                                    eX,...
                                    eY,...
                                    '',...
                                    '',...
                                    0,...
                                    'training-progress',...
                                    2000,...
                                    'cpu',...
                                    size(eX{1},1),...
                                    'last');
                                %%
                                trainedNet_FM;
                                trainedNet_FM = trainNetwork(trainedNet_FM,nX,YY);
                                %%
                                disp = 0;
                                close all
                                pY = {};
                                delta = [];
                                for e = 1:numel(nX)
                                    pY{e} = trainedNet_FM.predict(nX{e});
                                    trainedNet_FM.resetState();
                                    
                                    [~,ng] = max(double(pY{e}),[],1);
                                    
                                    
                                    [stepSig,frame] = extractFrameFromProb(pY{e}(2,:),.4,5,21);
                                    
                                    %{
                                    ng = (ng - 1) == 1;
                                    
                                    
                                    ng = [zeros(size(ng));ng;zeros(size(ng))];
                                    v = ng(2,end);
                                    ng(2,end) = 0;
                                    ng = imclearborder(ng);
                                    ng = ng(2,:);
                                    ng(end) = v;
                                    ng = bwareaopen(ng,5);
                                    ng = imclose(ng,strel('disk',40));
                                    %ng = ng > 3;
                                    %ng = double(pY{e}(2,:)) > .8;
                                    %ng = processProb(pY{e}(2,:),double(YY{e}) == 2);
                                    %{
                                    if sum(ng) > 30
                                        ng = bwareaopen(ng,5);
                                        ng = imclose(ng,strel('disk',5));
                                        tng = ng;
                                        ng = [zeros(size(ng));ng;zeros(size(ng))];
                                        
                                        ng = imcomplement(ng);
                                        
                                        ng = imfill(logical(ng),[2 size(ng,2)]);
                                        ng = ng(2,:);
                                        ng = tng == ng;
                                        %ng = imcomplement(ng(2,:));
                                    end
                                    %}
                                    %{
                                   %}
                                    %}
                                    %ng = ng -1;
                                    %ng = double(pY{e}(2,:)) > .5;
                                    
                                    
                                    if disp
                                        plot((double(YY{e}) > 1)+.1,'r');
                                        hold on;
                                        plot(stepSig==1,'g')
                                        hold off 
                                        drawnow
                                    end
                                   % pause(.1)
                                    
                                    
                                    pidx = find(stepSig==1);
                                    tidx = find(double(YY{e}) > 1);
                                    
                                    
                                    if ~isempty(pidx) & ~isempty(tidx)
                                        delta(e) = (pidx(1) - tidx(1));
                                        if delta(e) > 5
                                            if disp
                                                peak = zeros(size(stepSig));
                                                peak(pidx(1)) = 1;
                                                hold on
                                                plot(pY{e}(2,:),'m');
                                                plot(peak,'b');
                                                hold off
                                                waitforbuttonpress
                                            end
                                        end
                                    end
                                    mean(abs(delta))
                                    mean(delta)
                                    %e
                                end
                                %% view problems
                                [v,midx] = max(abs(delta))
%% this is a test script that I dont think i need any more
%{
func = cFlow('myGPUtest');
func.setMCRversion('v930');
func.setGPU(1);

layers = [ ...
    sequenceInputLayer(5)
    lstmLayer(5)
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];

net = func(nX,YY,layers,'gpu',500);


auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
func.submitDag(auth,50,50);
%}
%%

%% stack for hmm

EXTRA = predict(emergenceDynamic,single(BEST));
hmmBEST = BEST;
hmmBEST(:,:,3,:) = [];
hmmBEST = permute(hmmBEST,[1 3 2 4]);
Bsz = size(hmmBEST);

hmmBEST = reshape(hmmBEST,[prod(Bsz(1:2))  Bsz(3) Bsz(4)]);
fidx0 = find(HUM==0);
fidx1 = find(HUM==1);
[hmm] = makeChainRepeat_foEmergence(U0,C0,U1,C1,D,HARDLINE_HOLD,n,nCOMP,gmmNUM);

%% condor apply - local
for e = [13:numel(FileList)]
    try
       [finalScoreFM{e}] = detectEmergence(FileList{e},'',trainedNet_FM,151,zU,zE,zU,zE,11,Zmean,Zstd,CELLMASK,SELDIMS,toMatch,BOX,Kmu,Ksigma,Tscore); 
        
    catch ME
        
    end
end
%% raw test
        
%%
%{
%% backgroud adjustment
parfor e = 1:numel(FileList{8})
    I{e} = imread(FileList{8}{e});
    e
end
%% crop out small wood background
[sub BOX] = imcrop(I{1});
%%
close all
subI = [];
for e = 1:numel(I)
    subI(:,:,:,e) = imcrop(I{e},BOX);
    imshow(subI(:,:,:,e)/255,[]);
    title(num2str(e))
    drawnow
end
%%
close all
imshow(I{2},[])
%%
[optimizer, metric] = imregconfig('monomodal');
parfor e = 3:size(subI,4)
    tform{e} = imregtform(rgb2gray(subI(:,:,:,e)/255), rgb2gray(subI(:,:,:,2)/255), 'translation', optimizer, metric);
    e
end
%%
close all
N = 200;
outI = [];
sz = size(subI);

for k = 1:3
    outI(:,:,k) = imwarp(subI(:,:,k,N),tform{N},'OutputView',imref2d(sz(1:2)));
end
imshow(cat(3,subI(:,:,1:2,2),outI(:,:,3))/255)
%%
close all
sz = size(I{2});
for N = 3:numel(I)
    for k = 1:3
        nI{N}(:,:,k) = imwarp(I{N}(:,:,k),tform{N},'OutputView',imref2d(sz(1:2)));
    end
    N
end
imshow(cat(3,I{2}(:,:,1:2),nI{200}(:,:,3)))
%%
close all
imshow(cat(3,subI(:,:,1:2,2),subI(:,:,3,end)))
%}
%% deploy on DE
cornPopper = @(X)detectEmergence(X,'',trainedNet,151,zU,zE,zU,zE,11,Zmean,Zstd,CELLMASK,SELDIMS,toMatch,BOX,Kmu,Ksigma,Tscore);
%pF = partialFunction(cornPopper,'cornPopperNNapp');
%pF.publish();
% run local =- cornPopper(FileList{91})
%%
%cornPopper = @(X)detectEmergence_ver2(X,CNN_net,trainedNet_FM,nsz,zU,zE,CELLMASK,BOX,vidx,Tscore);
%cornPopper = @(X)detectEmergence_ver2(X,emergeNET3,trainedNet_FM,nsz,U_VV,E_VV,CELLMASK,BOX,vidx,Tscore);

%pF = partialFunction(cornPopper,'cornPopperNNapp');
%pF.publish();
% run local =- cornPopper(FileList{9});
%
%TestPath = '/iplant/home/jgustin/maizeData/coleoptileEmergence/20180130_Rack1_Camera1/';
%% download local for fixes
websave(['./cornPopperNNapp.mat'],'https://de.cyverse.org/dl/d/6EB4ADC9-1D8F-424E-95D5-03E2C715B142/cornPopperNNapp.mat');               
%%
load(['./cornPopperNNapp.mat']);
%% sept 12 upgrade
cornPopper = @(X)detectEmergence_ver2(X,emergeNET3,trainedNet_FM,nsz,U_VV,E_VV,CELLMASK,BOX,vidx,Tscore,GMModel);
pF = partialFunction(cornPopper,'cornPopperNNapp');
pF.publish();
%% Nov 21, 2019
websave(['./cornPopperNNapp.mat'],'https://de.cyverse.org/dl/d/6EB4ADC9-1D8F-424E-95D5-03E2C715B142/cornPopperNNapp.mat');     
load(['./cornPopperNNapp.mat']);
%% sept 19 upgrade make another
FilePath = '/mnt/tetra/nate/fixPOP/';
tFileList = {};
FileExt = {'tiff','TIF'};
tFileList = sdig(FilePath,tFileList,FileExt,1);
perD = .03;
MS = [];
toS = round(3*2.5*10000);
toI =  10;
str = 1;
TOT = toI*toS*numel(tFileList);
MS = zeros(TOT,3);
for s = 1:numel(tFileList)
    rnd = randi(numel(tFileList{s}),toI,1);
    for e = 1:numel(rnd)
        stp = str + toS - 1;
        I = double(imread(tFileList{s}{rnd(e)}))/255;
        %I = imresize(I,perD);
        sz = size(I);
        I = reshape(I,[prod(sz(1:2)) sz(3)]);
        fidx = randi(size(I,1),[toS 1]);
        MS(str:stp,:) = I(fidx,:);
        str = stp + 1;
        size(MS)
        e
    end
    s
end
%%
options = statset('Display','iter');
GMModel2 = fitgmdist(MS,4,'Options',options);
%cornPopperTEST = @(X)detectEmergence_ver2(X,emergeNET3,trainedNet_FM,nsz,U_VV,E_VV,CELLMASK,BOX,vidx,Tscore,GMModel2);
%%
S = functions(obj.func);
F = S.workspace{1}.fields;
emergeNET3 = S.workspace{1}.emergeNET3;
trainedNet_FM = S.workspace{1}.trainedNet_FM;
nsz = S.workspace{1}.nsz;
U_VV = S.workspace{1}.U_VV;
E_VV = S.workspace{1}.E_VV;
CELLMASK = S.workspace{1}.CELLMASK;
BOX = S.workspace{1}.BOX;
vidx = S.workspace{1}.vidx;
Tscore = S.workspace{1}.Tscore;

%cornPopperTEST = @(X)detectEmergence_ver2(X,emergeNET3,trainedNet_FM,nsz,U_VV,E_VV,CELLMASK,BOX,vidx,Tscore,GMModel2);

%%
FilePath = '/mnt/tetra/nate/junk/20180130_Rack1_Camera3/';
FilePath = '/mnt/tetra/nate/fixPOP/20180302_Rack2_Camera2/';
FilePath = '/mnt/tetra/nate/fixPOP/next3/20181101_Rack2_Camera6/';
FilePath = '/mnt/tetra/nate/fixPOP/next4/';

%FilePath = '/mnt/tetra/nate/fixPOP/next6/';
%FilePath = '/mnt/tetra/nate/fixPOP/next5/20181101_Rack2_Camera6/';
%FilePath = '/mnt/tetra/nate/fixPOP/20180521_Rack2_Camera3/';
tFileList = {};
FileExt = {'tiff','TIF'};
tFileList = sdig(FilePath,tFileList,FileExt,1);
%cornPopper(tFileList)
%%
for e = 1:numel(tFileList)
    I = imread(tFileList{e}{5});
    Lab = rgb2lab(I);
    tmp = Lab(:,:,2);
    fidx = find(tmp > 0);
    tmp(fidx) = tmp(fidx)/max(tmp(fidx));
    [py,px] = ksdensity(tmp(fidx));
    lidx = find(imdilate(py,strel('disk',10)) == py);
    [~,sidx] = sort(py(lidx),'descend');
    lidx = lidx(sidx(1:2));
    z = ones(size(py));
    z(lidx(1):lidx(2)) = py(lidx(1):lidx(2));
    [~,threshI(e)] = min(z);
    Y{e} = [px;py;z];
    thresh(e)=  px(threshI(e));
    e
    imshow(tmp > thresh(e),[]);
    drawnow
    
end
%% 
HI = [yi;xi];
fullLABELERS{1} = fullLabel;
fullLABELERS{2} = fullLabel2;
Tscore = [];
TcornPopper = @(X)detectEmergence(X,'',trainedNet_FM,151,zU,zE,zU,zE,11,Zmean,Zstd,CELLMASK,SELDIMS,toMatch,BOX,Kmu,Ksigma,HI,fullLABELERS,Tscore);
TcornPopper(FileList{1});
%%
close all
tough = zscore(rawX{f}(:,:,1),2);
plot(tough');
%% 
%Xdata = {};
for f = 1:numel(FileList)
    try
        %if isempty(Xdata{f})
            [cornScore2{f},Xdata{f},rawX{f}] = cornPopper(FileList{f}(1:200));
        %end
    catch ME
        CORNDEAD{f} = ME;
        getReport(ME)
        
    end
end
%% backup nX and YY
%% add to training data
for e = 1:numel(rawX)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(rawX{e})
        Tnum = e;
        [Tpth] = fileparts(FileList{Tnum}{1});
        fidx = strfind(Tpth,filesep);
        % make labels
        LL = {'A' 'B' 'C' 'D' 'E' 'F' 'G'};
        NL = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12'};
        HL = {'d' 'p'};
        LABELS = {};
        for e1 = 1:numel(HL)
            for e2 = 1:numel(LL)
                for e3 = 1:numel(NL)
                    LABELS{end+1} = [HL{e1} '-' LL{e2} NL{e3}];
                end
            end
        end
        
        
         TkeyT = lower([Tpth((fidx(end)+1):end) '-' LABELS{1}]);
         sc = handData_t2.get(TkeyT);
        if ~isempty(sc)
            LABELS = lower(LABELS);
            TESTValues = [];
            for l = 1:numel(LABELS)
                Tkey{l} = lower([Tpth((fidx(end)+1):end) '-' LABELS{l}]);
                sc = handData_t2.get(Tkey{l});

                sc = sc(1:min(250,numel(FileList{e})));
                TESTValues = [TESTValues ; sc];
            end

            TESTValues = reshape(TESTValues,[min(250,numel(FileList{e})) 168]);
            Tscore = [];
            for t = 1:size(TESTValues,2)
                fidx = find(TESTValues(:,t));
                if ~isempty(fidx)
                    Tscore(t) = fidx(1);
                else
                    Tscore(t) = 0;
                end
            end
            
            GGG{e} = Tscore;
            DELTA{e} = Tscore - cornScore{e};
            
            
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if ~isempty(rawX{e})
        for tr = 1:size(rawX{e},3)
            if Tscore(tr) ~= 0 
                newSig = convertToStrain(rawX{e}(:,:,tr),11,10);
                nX{end+1} = newSig;
                YY{end+1} = categorical(TESTValues(1:(end-1),tr)');
            end
        end
    end
end
%%

%% generate large EMERGE SET
eX = {};
eY = [];
for e = 1:numel(rawX)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(rawX{e})
        Tnum = e;
        [Tpth] = fileparts(FileList{Tnum}{1});
        fidx = strfind(Tpth,filesep);
        % make labels
        LL = {'A' 'B' 'C' 'D' 'E' 'F' 'G'};
        NL = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12'};
        HL = {'d' 'p'};
        LABELS = {};
        for e1 = 1:numel(HL)
            for e2 = 1:numel(LL)
                for e3 = 1:numel(NL)
                    LABELS{end+1} = [HL{e1} '-' LL{e2} NL{e3}];
                end
            end
        end
        
        
         TkeyT = lower([Tpth((fidx(end)+1):end) '-' LABELS{1}]);
         sc = handData_t2.get(TkeyT);
        if ~isempty(sc)
            LABELS = lower(LABELS);
            TESTValues = [];
            for l = 1:numel(LABELS)
                Tkey{l} = lower([Tpth((fidx(end)+1):end) '-' LABELS{l}]);
                sc = handData_t2.get(Tkey{l});

                sc = sc(1:min(250,numel(FileList{e})));
                TESTValues = [TESTValues ; sc];
            end

            TESTValues = reshape(TESTValues,[min(250,numel(FileList{e})) 168]);
            Tscore = [];
            for t = 1:size(TESTValues,2)
                fidx = find(TESTValues(:,t));
                if ~isempty(fidx)
                    Tscore(t) = fidx(1);
                else
                    Tscore(t) = 0;
                end
            end
            
            GGG{e} = Tscore;
            %DELTA{e} = Tscore - cornScore{e};
            
            
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if ~isempty(rawX{e})
        for tr = 1:size(rawX{e},3)
            %newSig = convertToStrain(rawX{e}(:,:,tr),11,10);
            newSig = convertToDoesEmergeSignal(rawX{e}(:,:,tr),11,10);

            %newSig = rawX{e}(:,:,tr);
            %newSig = bsxfun(@minus,newSig,mean(newSig,2));
            eX{end+1} = newSig;
            EP = (TESTValues(size(newSig,2),tr)') > 0;
            eY = [eY ; categorical(EP)];
        end
    end
end
%% backup the eX and eY
eXBK = eX;
eYBK = eY;
%% train emergeNET

options = trainingOptions('sgdm',...
            'InitialLearnRate',.01,...
            'MaxEpochs',1000,...
            'Shuffle','every-epoch',...
            'Verbose',true,...
            'ExecutionEnvironment','cpu',...
            'Plots','training-progress');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EmergeNetlayers = [ ...
    sequenceInputLayer(size(eX{1},1))
    lstmLayer(7,'OutputMode','last')
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];

%trainedNet_FM = trainNetwork(SIGGGY,YY,layers,options);

emergeNET = trainNetwork(eX,eY,EmergeNetlayers,options);
%% watch high res
I = imread(FileList{f}{2});
[tmpJ,tmpBOX] = imcrop(I);
cnt = 1;
for e = 2:2:numel(FileList{f})
    I = imread(FileList{f}{e});
    JJ(:,:,:,cnt) = imcrop(I,tmpBOX);
    cnt = cnt + 1;
end
%%
close all
I1 = rgb2gray(JJ(:,:,:,1));
Iend = rgb2gray(JJ(:,:,:,end));
for e = 1:size(JJ,4)
    Im = rgb2gray(JJ(:,:,:,e));
    KKM = cat(3,I1,Im,Iend);
    imshow(KKM,[])
    title(num2str(e));
    drawnow
end
%% regrade
reSig = convertToStrain(rawX{f},11,10);
for e = 1:size(reSig,3)
    tmpP = trainedNet_FM.predict(reSig(:,:,e));
    [~,tmpP] = max(tmpP,[],1);
    tmpP = tmpP == 2;
    reGrade
end
%% FINAL GRADE
DELTA = {};
for e = 1:13%numel(finalScore2)
    if ~isempty(finalScore2{e})
        
        Tnum = e;
        [Tpth] = fileparts(FileList{Tnum}{1});
        fidx = strfind(Tpth,filesep);
        % make labels
        LL = {'A' 'B' 'C' 'D' 'E' 'F' 'G'};
        NL = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12'};
        HL = {'d' 'p'};
        LABELS = {};
        for e1 = 1:numel(HL)
            for e2 = 1:numel(LL)
                for e3 = 1:numel(NL)
                    LABELS{end+1} = [HL{e1} '-' LL{e2} NL{e3}];
                end
            end
        end
        
        
         TkeyT = lower([Tpth((fidx(end)+1):end) '-' LABELS{1}]);
         sc = handData_t2.get(TkeyT);
        if ~isempty(sc)
            LABELS = lower(LABELS);
            TESTValues = [];
            for l = 1:numel(LABELS)
                Tkey{l} = lower([Tpth((fidx(end)+1):end) '-' LABELS{l}]);
                sc = handData_t2.get(Tkey{l});

                sc = sc(1:min(250,numel(FileList{e})));
                TESTValues = [TESTValues ; sc];
            end

            TESTValues = reshape(TESTValues,[min(250,numel(FileList{e})) 168]);
            Tscore = [];
            for t = 1:size(TESTValues,2)
                fidx = find(TESTValues(:,t));
                if ~isempty(fidx)
                    Tscore(t) = fidx(1);
                else
                    Tscore(t) = 0;
                end
            end
            GGG{e} = Tscore;
            DELTA{e} = Tscore - cornScore{e};
        end
        
    end
end
%% stats on DELTA
MT = [];
P = [];
close all

for e = 1:numel(DELTA)
    if ~isempty(DELTA{e})
        tmp = DELTA{e};
        tmp(finalScore2{e}==0) = 0;
        tmp2 = GGG{e};
        tmp2(finalScore2{e}==0) = 0;
        tmp3 = finalScore2{e};
        tmp3(finalScore2{e}==0) = 0;
        MT = [MT;tmp(:)];
        P = [P;[tmp3(:) tmp2(:) e*ones(size(tmp(:)))]]; 
    end
end
plot3(P(:,1),P(:,2),P(:,3),'.')
title([num2str(mean(abs(MT))) '--' num2str(corr(P(:,1),P(:,2)))]);
kidx = P(:,1) ~= 0 & P(:,2) ~= 0;
mean(abs(MT(kidx)))
hold on
plot3(P(kidx,1),P(kidx,2),P(kidx,3),'ro')
%% apply BEST
%[score] = detectEmergence(FileList{5},emergenceNet,BESTEmergence,151,zU,zE,11,[],[]);
%[score1,score2,raw1,raw2] = detectEmergence(FileList{1},emergenceNet3,emergenceDynamic,151,zU,zE,zU,zE,11,Zmean,Zstd,CELLMASK,SELDIMS,toMatch,BOX,Tscore);
[finalScore] = detectEmergence(FileList{2},'',emergenceDynamic,151,zU,zE,zU,zE,11,Zmean,Zstd,CELLMASK,SELDIMS,toMatch,BOX,Tscore);
%[score1,score2,raw1,raw2] = detectEmergence(FileList{1},'',emergenceDynamic,151,zU,zE,11,Zmean,Zstd,CELLMASK,SELDIMS,toMatch,BOX,Tscore);
%%
func = @(X)detectEmergence(X,emergenceNet3,emergenceDynamic,151,zU,zE,zU,zE,11,Zmean,Zstd,CELLMASK,SELDIMS,toMatch,BOX,Tscore);
pf = partialFunction(func,'cornPopperNNapp');
pf.publish();
%%
layers = [imageInputLayer([size(rZ,1) size(rZ,2) 3]);
          convolution2dLayer([12 12],20);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          convolution2dLayer([5 5],20);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(2);
          softmaxLayer();
          classificationLayer()];
layers = [imageInputLayer([size(rZ,1) size(rZ,2) 3]);
          convolution2dLayer([21 21],7);
          reluLayer();
          maxPooling2dLayer(8,'Stride',2);
          %convolution2dLayer([5 5],20);
          %reluLayer();
          %maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(2);
          softmaxLayer();
          classificationLayer()];
layers = [imageInputLayer([size(rZ,1) size(rZ,2) 3]);
          convolution2dLayer([11 11],9);
          reluLayer();
          maxPooling2dLayer([8 8],'Stride',6);
          convolution2dLayer([3 3],4);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(2);
          softmaxLayer();
          classificationLayer()];
%selIdx = randperm(numel(TYPE));
%selIdx = 1:249:numel(TYPE);
%subX = imgSIM(:,:,:,1:4:end);
%subY = TYPE(1:4:end);
options = trainingOptions('sgdm','MaxEpochs',5,'InitialLearnRate',0.0001,'ExecutionEnvironment','parallel');
emergenceNet3 = trainNetwork(single(rZ),categorical(TYPE),layers,options);
%emergenceNet4 = trainNetwork(subX,categorical(subY),layers,options);
%%
TESTN = predict(emergenceNet4,imgSIM(:,:,:,1:3*249));
%%
[score] = detectEmergence(FileList{3},emergenceNet,151);
%% regress 
layers = [ ...
    imageInputLayer([size(zCfinal,1) size(zCfinal,2) 1])
    convolution2dLayer([2 15],10)
    reluLayer
    maxPooling2dLayer([2 2],'Stride',1);
    %convolution2dLayer([10 1],5)
    %reluLayer
    fullyConnectedLayer(1)
    regressionLayer];
options = trainingOptions('sgdm','MaxEpochs',700,'InitialLearnRate',0.00001,'ExecutionEnvironment','parallel');
nVALUE = VALUE;
nVALUE(find(isinf(VALUE))) = 0;
UnVALUE = mean(nVALUE);
nVALUE = nVALUE - UnVALUE;
RtrainedNet = trainNetwork(single(zCfinal),nVALUE,layers,options);

%% TEST BEST
close all
[YBEST PROB] = classify(BESTEmergence,single(BEST));
plot(YBEST)
%% test prediction for NON regression
TEST = predict(emergenceNet,single(rZ));
%% test regression
close all
RTEST = predict(RtrainedNet,single(zCfinal)) + UnVALUE;
plot(VALUE,RTEST,'.')
%% reshape
nsz = 101;
rZ = zeros([nsz nsz 3 numel(Z)]);
for e = 1:numel(Z)
    rZ(:,:,:,e) = Z{e};
    e
end
%% check a test
Tnum = 2;
[Tpth] = fileparts(FileList{Tnum}{1});
fidx = strfind(Tpth,filesep);
% make labels
LL = {'A' 'B' 'C' 'D' 'E' 'F' 'G'};
NL = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12'};
HL = {'d' 'p'};
LABELS = {};
for e1 = 1:numel(HL)
    for e2 = 1:numel(LL)
        for e3 = 1:numel(NL)
            LABELS{end+1} = [HL{e1} '-' LL{e2} NL{e3}];
        end
    end
end
LABELS = lower(LABELS);
TESTValues = [];
for l = 1:numel(LABELS)
    Tkey{l} = lower([Tpth((fidx(end)+1):end) '-' LABELS{l}]);
    sc = handData_t2.get(Tkey{l});
    TESTValues = [TESTValues ; sc];
end
TESTValues = reshape(TESTValues,[250 168]);
Tscore = [];
for t = 1:size(TESTValues,2)
    fidx = find(TESTValues(:,t));
    if ~isempty(fidx)
        Tscore(t) = fidx(1);
    else
        Tscore(t) = 0;
    end
end
%%
DELTA = abs(double(score2) - Tscore);
[~,sidx] = sort(DELTA);
%%

%%
Rscore = [];
mag = 5;
for t = 1:size(score,1)
    fidx = find(squeeze(score(t,2,:) > .5));
    if ~isempty(fidx)
        Rscore(t) = mag*fidx(1);
    else
        Rscore(t) = -1;
    end
end
close all
DELTA = abs(Tscore - Rscore);
[J,sidx] = sort(DELTA);
figure;
plot(Tscore,Rscore,'.')
hold on
plot(Tscore(sidx(end)),Rscore(sidx(end)),'ro');
figure;
plot(squeeze(score(sidx(end),2,:)))
fidx = find(Rscore == 10);
figure;
plot(squeeze(score(fidx(1),2,:)))

%%
        
        
        SZ = size(STACK);
        STACK = reshape(STACK,[prod(SZ(1:2)) SZ(3:5)]);
        STACK = permute(STACK,[1 3 2 4]);
        
        
        
        close all
        for cb = 1:size(STACK,5)
            for t = 1:size(STACK,4)
                imshow(STACK(:,:,:,t,cb),[]);
                drawnow
            end
        end
        
        

%%
close all
for e = 1:size(CheckerBoard,1)
    [cameraParameters,imagePoints,worldPoints] = getCameraPara(CheckerBoard{e,1}{1});
    for s = 1:100
        tic
       
        [rI{e}(:,:,s)] = getRectifiedImage(CheckerBoard{e,1}{s},cameraParameters,imagePoints,worldPoints);
        toc
    end
    
    %{
    SZ = size(subI{e});
    tform = fitgeotrans(imagePoints,worldPoints,'affine');
    
    
    
    cI = imwarp(subI{e},tform);
    
    
    
    E = edge(rgb2gray(subI{e}),'Canny');
    [H,T,R] = hough(E','theta',linspace(-15,15,200));
    P  = houghpeaks(H,40,'NHoodSize',[101 11]);
    lines = houghlines(E',T,R,P,'FillGap',600,'MinLength',400);
    
    [H2,T2,R2] = hough(E,'theta',linspace(-20,20,200));
    P2  = houghpeaks(H2,40,'NHoodSize',[101 11]);
    lines2 = houghlines(E,T2,R2,P2,'FillGap',600,'MinLength',400);
    
    close all
    figure, imshow(subI{e}), hold on
    max_len = 0;
    for k = 1:length(lines)
       xy = [lines(k).point1; lines(k).point2];
       xy = flipdim(xy,2);
       
       plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

       % Plot beginnings and ends of lines
       plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
       plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

       % Determine the endpoints of the longest line segment
       len = norm(lines(k).point1 - lines(k).point2);
       if ( len > max_len)
          max_len = len;
          xy_long = xy;
       end
    end
    
    max_len = 0;
    for k = 1:length(lines2)
       xy = [lines2(k).point1; lines2(k).point2];
       
       plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','blue');

       % Plot beginnings and ends of lines
       plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
       plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

       % Determine the endpoints of the longest line segment
       len = norm(lines2(k).point1 - lines2(k).point2);
       if ( len > max_len)
          max_len = len;
          xy_long = xy;
       end
    end
    
    
    plot(imagePoints(:,1),imagePoints(:,2),'k*');
    
    %{
    imshow(subI{e},[]);
    drawnow
    %}
    drawnow
    %}
end

%% find frist images
for e = 1:numel(FileList)
    [pth,nm,ext] = fileparts(FileList{e});
    nm = str2num(nm);
    if nm == 1
        kp(e) = true;
    else
        kp(e) = false;
    end
end
FileList = FileList(kp);
%%
dlP = '/home/nate/Downloads/JEFF/';
mkdir(dlP);
cnt = 1;
for e = 30:numel(FileList)
    CMD = ['iget ' FileList{e} ' ' dlP num2str(cnt) '.tif'];
    cnt = cnt + 1;
    system(CMD)
    CMD
end
%%
FilePath = 'W:\';
FilePath = '/home/nate/Downloads/JEFF/';
FileList = {};
FileExt = {'tif','TIF'};
FileList = gdig(FilePath,FileList,FileExt,1)
%%
I = {};
for e = 1:numel(FileList)
    try
        I{end+1} = imread(FileList{e});
        imshow(I{end},[]);
        drawnow
    catch
    end
    
end
%%
un = CheckerBoard.checkerBoard;
for e = 1:size(un,1)
    CB{e} = un(e,:);
end
CheckerBoard.checkerBoard = CB'
%%
options = trainingOptions('sgdm','MaxEpochs', 10,'ExecutionEnvironment','auto');%,'InitialLearnRate', 1e-6);
%layers = data.layers;
layers = [imageInputLayer([75 75 3])
          convolution2dLayer([15,15],10)
          reluLayer
          maxPooling2dLayer([3 3],'Stride',2)
          convolution2dLayer([7,7],3)
          reluLayer
          maxPooling2dLayer([3 3],'Stride',2)
          fullyConnectedLayer(2)
          softmaxLayer()
          classificationLayer()];


checkerDetector = trainFasterRCNNObjectDetector(CheckerBoard, layers, options,'SmallestImageDimension',400,'NumStrongestRegions',10);


%%
img = imread(CheckerBoard{1,1}{1});
[bbox, score, label] = detect(checkerDetector, img,'SelectStrongest',false);

detectedImg = insertShape(img, 'Rectangle', bbox);
figure
imshow(detectedImg)


%% get values for a