%% mount max folder local
remoteBase = '/iplant/home/mbraud/UIUC_2019_Tassel_Images/';
localBase = '/home/nate/mntCyverse/tasselMax/';
mountCyVerse(remoteBase,localBase)
%sudo ./mountCyverse mbraud/UIUC_2019_Tassel_Images/ /home/nate/mntCyverse/tasselMax/
%% pull data here - this is a mess but works
FilePath = '/iplant/home/mbraud/UIUC_2019_Tassel_Images/';
FileList = {};
FileExt = {'tiff','tif','jpg','TIF','TIFF','PNG'};
FileList = idig(FilePath,FileList,FileExt);
targetD = '/mnt/spaldingdata/nate/octerineDataStore/tassel/UIUC_2019_Tassel_images/';
xfer_get_p(FileList,targetD);
%% dig for image files - dig at mount point
FilePath = '/home/nate/mntCyverse/tasselMax/';
FileList = {};
FileExt = {'jpg'};
FileList = fdig(FilePath,FileList,FileExt,1);
%% dig for files in pile - dig over local pile
FilePath = '/mnt/spaldingdata/nate/octerineDataStore/tassel/UIUC_2019_Tassel_images/';
FileList = {};
FileExt = {'jpg'};
rawFileList = fdig(FilePath,FileList,FileExt,1);
%% attach local pile data set
dataPath = '/mnt/spaldingdata/nate/octerineDataStore/tassel/UIUC_2019_Tassel_images/';
project.attachDataCollection('localPile_test',dataPath);
%% attach function of background testing
func = @(in,out)getTasselMask(in,10,out);
project.attachFunction('tasselBackground',func);
%% run the background test
project.runCollection('localPile_test',Inf,'tasselBackground');
%% run the background test
project.runCollection('localPile_test',1200);
%% run the background test
project.runCollection('localPile_test',inf,'main_ver2');
%% profile the main
project.profile('localPile_test',3);
%% get function
func = project.getFunction('main');
func.operate(rawFileList{3},makeTempLocation())
%% load from the background results and see if we an train a AI system
FilePath = '/mnt/spaldingdata/nate/octerineDataStore/results/{projectName_TIPS}/{functionName_tasselBackground}/{createDate_2020-04-28-12-16-55}/';
FileList = {};
FileExt = {'tif'};
outFileList = fdig(FilePath,FileList,FileExt);
samN = 8000;
sam0 = [];
sam1 = [];
sam2 = [];
sam3 = [];
for e = 1:numel(outFileList)
    if contains(outFileList{e},'overlay')

        [p,n,ex] = fileparts(outFileList{e});
        maskName = strrep(outFileList{e},'overlay','mask');
        fidx = strfind(p,filesep);
        in = p((fidx(end)+1):end);
        rawName = find(contains(rawFileList,in));
        
        
        
        
        
        I = imread(rawFileList{rawName});
        oI = I;
        M = imread(maskName);
        O = imread(outFileList{e});
        DT = bwdist(M);
        N = (DT < 30) & ~ M;
        NN = ~M & ~N;
        %{
        % bester
        sample = sampleColorFeatures1(I);
        sample = sampleNhoodFeatures(I,sample);
        sample = sampleGrayGradFeatures(I,sample);
        %}
        
        idx0 = find(M==0);
        idx1 = find(M==1);
        
        idx2 = find(N==1);
        idx3 = find(NN==1);
        
        
        idx0 = idx0(randperm(numel(idx0)));
        idx1 = idx1(randperm(numel(idx1)));
        idx2 = idx2(randperm(numel(idx2)));
        idx3 = idx3(randperm(numel(idx3)));


        szI = size(I);
        I = reshape(I,[prod(szI(1:2)) szI(3)]);

    
        if exist('Mdl')
            y = Mdl.predict(double(I));
            y = reshape(y,szI(1:2));
            
            y2 = net(double(I)');
            y2o = reshape(y2(2,:),szI(1:2));
            
            y3 = tree.predict(double(I));
            y3 = reshape(y3,szI(1:2));
            
            out = flattenMaskOverlay(oI,logical(y));
            out2 = flattenMaskOverlay(oI,logical(y2o>.5));
            out3 = flattenMaskOverlay(oI,logical(y3>.5));
            out = [O,out,out2,out3];
        end
        
        
        

        tmpN0 = min(numel(idx0),samN);
        tmpN1 = min(numel(idx1),samN);
        tmpN2 = min(numel(idx2),samN);
        tmpN3 = min(numel(idx3),samN);
        
        
        idx0 = idx0(1:tmpN0);
        idx1 = idx1(1:tmpN1);
        idx2 = idx2(1:tmpN2);
        idx3 = idx3(1:tmpN3);

        sam0 = [sam0;I(idx0,:)];
        sam1 = [sam1;I(idx1,:)];
        sam2 = [sam2;I(idx2,:)];
        sam3 = [sam3;I(idx3,:)];
        size(sam0)
        size(sam1)
        e
    end
end
%% train bayes
grp1 = 1*ones(size(sam1,1),1);
grp0 = 0*ones(size(sam0,1),1);
Y = [grp0;grp1];
X = [sam0;sam1];
Mdl = fitcnb(double(X),Y);
%% train ANN
Y = [grp0;grp1];
X = double([sam0;sam1]);
Y = full(ind2vec(Y'+1))';
net = patternnet(5);
net = train(net,X',Y','useParallel','yes');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% near the tassel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% train bayes
grp2 = 0*ones(size(sam2,1),1);
grp3 = 0*ones(size(sam3,1),1);
grp1 = 1*ones(size(sam1,1),1);
Y = [grp2;grp3;grp1];
X = [sam2;sam3;sam1];
Mdl = fitcnb(double(X),Y);
%% train ANN
Y = [grp2;grp3;grp1];
X = [sam2;sam3;sam1];
Y = full(ind2vec(Y'+1))';
net = patternnet(5);
net = train(net,double(X)',Y','useParallel','yes');
%% train dTree
Y = [grp2;grp3;grp1];
X = [sam2;sam3;sam1];
% Y = sim(net,double(X)');Y = Y(2,:)' > .5;
tree = fitctree(double(X),Y);
%tree = fitctree(double(X),Y,'OptimizeHyperparameters','auto');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make list of AI methods - functionStack
aiMethods(1) = formalFunc('naiveBayes',@(x)Mdl.predict(x));
aiMethods(2) = formalFunc('ANN',@(x)sim(net,x'));
aiMethods(3) = formalFunc('cTree',@(x)tree.predict(x));
%% 
for e = 1%:numel(failList)
    fileName = failList(e).fileName;
    [tasselMask,I] = getTasselMask_ver2(fileName,10,aiMethods);
    %montage(failList(1).maskName,'size',[2 2]);
end
%% make TIPS networks
genFunction(net,'TIPSnet')
%% 
in = tiFileList{10};
oPath = './JUNK/';
out = TIPS_noBackground_verCondor(in,oPath,tree,Mdl);
%% try onto grid
FilePath = '/iplant/home/mbraud/UIUC_2019_Tassel_Images/';
iFileList = {};
FileExt = {'tiff','tif','jpg','TIF','TIFF','PNG'};
iFileList = idig(FilePath,iFileList,FileExt);

tiFileList = issueBulkTicket(iFileList);
tipsFunc = cFlow('TIPS_noBackground_verCondor');
tipsFunc.setMCRversion('v930');
tipsFunc.addDirectoryMap('output>/mnt/spaldingdata/nate/joeTest/');
for e = 1:numel(tiFileList)
    tipsFunc(tiFileList{e},'./output/',tree,Mdl);
end
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
tipsFunc.submitDag(auth,50,50);
%% loop and find landscape images
for e = 1:numel(rawFileList)
    imgInfo = imfinfo(rawFileList{e});
    sz = [imgInfo.Height imgInfo.Width];
    if sz(1) < sz(2)
        rawFileList{e}
    end
end
%% attach main function
bgFunc = @(x,out)getTasselMask_ver2(x,10,aiMethods,out);
func = @(in,out)TIPS_noBackground_ver2(in, out, bgFunc);
%%
project.attachFunction('main_ver2',func);
%% dig at CyVerse
FilePath = '/iplant/home/mbraud/UIUC_2019_Tassel_Images/';
iFileList = {};
FileExt = {'tiff','tif','jpg','TIF','TIFF','PNG'};
iFileList = idig(FilePath,iFileList,FileExt);
tiFileList = issueBulkTicket(iFileList);
%%
tipsFunc = cFlow('TIPS_noBackground_ver2');
tipsFunc.setMCRversion('v930');
tipsFunc.addDirectoryMap('output>/mnt/spaldingdata/nate/joeTest2/');
tipsFunc(tiFileList{20},'./output',bgFunc);
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
tipsFunc.submitDag(auth,50,50);
%% try anonyous onto grid
funky = @(x,y)letsAdd(x,y);
funky(1,3)
funkyFlow = cFlow('applyThisFunky');
funkyFlow.setMCRversion('v930');
%tipsFunc.addDirectoryMap('output>/mnt/spaldingdata/nate/joeTest/');
funkyFlow(1,3,funkyFlow);
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
funkyFlow.submitDag(auth,50,50);
%% run main with 5 images
project.runCollection('localPile_test',5);
%% this will run TIPS no background
TIPS_noBackground(FileList{4},'./outputT/');
%%
basePath = './outputT/';
parfor e = 1:5
    fileName = FileList{e};
    [p,n,ext] = fileparts(fileName);
    oPath = [basePath n filesep];
    getTasselMask(fileName,10,oPath);
end
%% read of some data
% rough context numbers at a pixel level
% such as the of the n-hood-std-nhood
close all
N = 30;
pidx = randperm(numel(FileList));
pidx = pidx(1:N);
reSZ = .25;
samN = 50000;
colorSamples = [];

for e = 1:numel(pidx)
    
    I = imread(FileList{pidx(e)});
    I = double(I)/255;
    
    uI = imfilter(I,fspecial('disk',31),'replicate');
    
    
    if exist('refImage')
        I = imhistmatch(I,refImage);
    end
    %Lab = rgb2lab(I);
    
    if exist('gm')
        %I = imresize(I,reSZ);
        Lab = rgb2lab(I);
        
        szI = size(I);
        J = reshape(I,[prod(szI(1:2)) szI(3)]);
        Ju = reshape(uI,[prod(szI(1:2)) szI(3)]);
        JLAB = reshape(Lab,[prod(szI(1:2)) szI(3)]);
        data = [J,Ju,JLAB];
        
        pre = gm.cluster(data);
        pre = reshape(pre,szI(1:2));
        preRGB = label2rgb(pre);
        imshow(preRGB,[]);
        waitforbuttonpress
        
        
        
        for c = 1:cl
            tmpPRE = pre;
            fidx = (pre == c);
            subData = data(fidx,:);
            pre_sub = gm_sub{c}.cluster(subData);
            pre_sub = pre_sub + cl;
            tmpPRE(fidx) = pre_sub;
            preRGB = label2rgb(tmpPRE);
            imshow(preRGB,[]);
            waitforbuttonpress
        end
        
        %{
        fidx = (pre == 3);
        subData = data(fidx,:);
        
        ridx = randperm(size(subData,1));
        subSample = subData(ridx(1:samN),:);
        sampleData_2 = [sampleData_2,subSample];
        %}
    end
    
    I = imresize(I,reSZ);
    uI = imresize(uI,reSZ);
    Lab = rgb2lab(I);
    
    
    szI(e,:) = size(I);
    
    samI = reshape(I,[prod(szI(e,1:2)) szI(e,3)]);
    samuI = reshape(uI,[prod(szI(e,1:2)) szI(e,3)]);
    samLab = reshape(Lab,[prod(szI(e,1:2)) szI(e,3)]);
    
    
    
    ridx = randperm(size(samI,1));
    
    
   
    imageSample = pidx(e)*ones(samN,1);
    pixelLocation = ridx(1:samN)';
    
    
    colorSample = samI(ridx(1:samN),:);
    colorSampleU = samuI(ridx(1:samN),:);
    colorSampleLAB = samLab(ridx(1:samN),:);
    
   
    
    
    
    sampleData = [imageSample,pixelLocation,colorSample,colorSampleU,colorSampleLAB];
    
    
    colorSamples = [colorSamples;sampleData];
    
    imshow(I,[]);
    drawnow
end
%% fit distribution to pull out the background
options = statset('Display','iter');
data = double(colorSamples(:,3:end));
cl = 3;
gm = fitgmdist(data,cl,'Options',options);
clData = gm.cluster(data);
for c = 1:cl
    fidx = (clData == c);
    subData = data(fidx,:);
    gm_sub{c} = fitgmdist(subData,3,'Options',options);
end
%% make reference image from samples
refImage = [];
for k = 1:3
    refImage = cat(3,refImage,colorSamples(:,k+2));
end
%% make some image masks
oPath = '/mnt/tetra/nate/tasselMasks/';
mmkdir(oPath);
parfor e = 1:numel(FileList)
    makeImageMask(FileList{e},oPath,10);
end
%% get a mask file list of the local data
FilePath = '/mnt/tetra/nate/tasselMasks/';
maskFileList = {};
FileExt = {'tif'};
maskFileList = fdig(FilePath,maskFileList,FileExt,1);
%% loop over local data - gather mask features for gmm below
% mask features are the first layer of clustering
% the maskfeature clusters help choose which type of skeleton point
% is being sampled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

domainW = 50;
domainW_num = 101;
domainL = 50;
domainL_num = 101;

% make sample domain
%%%%%%%%%%%%%%%%%%%%%%%%%
[x1,x2] = ndgrid(linspace(-domainW,domainW,domainW_num),linspace(-domainL,domainL,domainL_num));
szX = size(x1);
x = [x1(:) x2(:) ones(size(x1(:)))];



fT = @(P,gD)simpleAffine(P,gD);
cnt = 1;
fm = [];

skelSamSZ = 2000;

% for each mask image
for e = 1:numel(maskFileList)
    
    % read themask
    M = imread(maskFileList{e});
    % get the skeleton
    skeletonI = bwmorph(M,'skeleton',inf);
    % find the skeleton points
    skel = [];
    [skel(:,2),skel(:,1)] = find(skeletonI);
    % randperms the points
    ridx = randperm(size(skel,1));
    skel = skel(ridx,:);
    skel = skel(1:skelSamSZ,:);
    
    
    M = double(M);
    
    maskBasis = basisF(M);
    skeletonBasis = basisF(skeletonI);
    
    tic
    tmpF = [];
    % loop over the skeleton points
    parfor s = 1:size(skel,1)
        try
            boxW = 100;
            tic;
            samplePoint = skel(s,:);
            %samplePoint = basisT.affine2D(samplePoint');
            samplePoint = basisT(samplePoint');
            box = samplePoint.point2box([boxW boxW]);
            sampleMask = maskBasis*box;
            sampleSkeleton = skeletonBasis*box;


            tmpM = sampleMask.E;

            tmpM = tmpM > .5;
            hsz = hsize(tmpM);
            tmpM = imfill(~tmpM,hsz+.5,8) & tmpM;

            tmpR = regionprops(logical(tmpM),'MajorAxis','MinorAxis','Area','Perimeter','ConvexArea','Eccentricity');
            tmpF(s,:) = [tmpR.Area,tmpR.MajorAxisLength,tmpR.MinorAxisLength,tmpR.Eccentricity,tmpR.ConvexArea,tmpR.Perimeter];
            toc;
        catch ME
            ME
        end
    end
   rm = find(all(tmpF == 0,2));
   tmpF(rm,:) = [];
    
    toc
    e
    
    fm = [fm;tmpF];
end
%% fit gmm to mask features
%%%%%%%%%%%%%%%%%%%%%%%%
grps = 3;
zfm = zscore(fm);
options = statset('Display','iter','MaxIter',300);
gm = fitgmdist(zfm,grps,'Options',options);
%% loop over the data again and cluster the skeleton points
% make sample domain
[x1,x2] = ndgrid(linspace(-domainW,domainW,domainW_num),linspace(-domainL,domainL,domainL_num));
szX = size(x1);
x = [x1(:) x2(:) ones(size(x1(:)))];

fT = @(P,gD)simpleAffine(P,gD);

masterPointList = [];
for e = 1:10%numel(maskFileList)
    fprintf(['Starting image:' num2str(e) '\n']);tic;
    
    % read the mask image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M = imread(maskFileList{e});
    % make the skeleton
    skeletonI = bwmorph(M,'skeleton',inf);
    % find the skeleton points
    skel = [];[skel(:,2),skel(:,1)] = find(skeletonI);
    
    
    
    
    % sample the skeleton points for binary mask features
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M = double(M);
    tic
    tmpF = [];
    % for each skeleton point
    parfor s = 1:size(skel,1)
        % get the skeleton point
        samplePoint = skel(s,:);
        
        
        
        samplePoint = basisT(samplePoint');
        box = samplePoint.point2box([boxW boxW]);
        sampleMask = maskBasis*box;
            
            
        % sample the point
        sampleP = funcP(samplePoint,fT);
        sampleP.getT([]);
    
        tmpM = sampleP.getF(M,x,szX);
        
        
        tmpM = tmpM > .5;
        hsz = hsize(tmpM);
        tmpM = imfill(~tmpM,hsz+.5,8) & tmpM;
        
        tmpR = regionprops(logical(tmpM),'MajorAxis','MinorAxis','Area','Perimeter','ConvexArea','Eccentricity');
        tmpF(s,:) = [tmpR.Area,tmpR.MajorAxisLength,tmpR.MinorAxisLength,tmpR.Eccentricity,tmpR.ConvexArea,tmpR.Perimeter];
    end
    toc
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kidx = gm.cluster(zscore(tmpF));
    imshow(M,[]);hold on
    for k = 1:grps
        plot(skel(kidx==k,1),skel(kidx==k,2),CL{k})
    end
    
    tmpPointList = [e*ones(size(skel,1),1) skel kidx];
    
    
    masterPointList = [masterPointList;tmpPointList];
end
%% sub cluster group-3 based on internal distance matrix
% source=image3 target=image4
grp = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = 3;
sourceIdx = find(masterPointList(:,1) == s & masterPointList(:,4) == grp);
source_PixelList = masterPointList(sourceIdx,2:3);
maskPath = '/mnt/tetra/nate/tasselMasks/';
source_rgbImage = strrep(FileList{s},'/home/nate/mntCyverse/tasselMax/','/iplant/home/mbraud/UIUC_2019_Tassel_Images/');
[~,nm] = fileparts(source_rgbImage);
source_maskImage = [maskPath nm '.tif'];
% setup the source
sourceData.rgb_fileName = source_rgbImage;
sourceData.mask_fileName = source_maskImage;
sourceData.pixelList = source_PixelList;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 4;
targetIdx = find(masterPointList(:,1) == t & masterPointList(:,4) == grp);
target_PixelList = masterPointList(targetIdx,2:3);
maskPath = '/mnt/tetra/nate/tasselMasks/';
target_rgbImage = strrep(FileList{t},'/home/nate/mntCyverse/tasselMax/','/iplant/home/mbraud/UIUC_2019_Tassel_Images/');
[~,nm] = fileparts(target_rgbImage);
target_maskImage = [maskPath nm '.tif'];
% setup the target
targetData.rgb_fileName = target_rgbImage;
targetData.mask_fileName = target_maskImage;
targetData.pixelList = target_PixelList;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
domainW = 50;
domainL = 50;
domainW_num = round(101);
domainL_num = round(101);
domainData.xData = [domainW domainW_num];
domainData.yData = [domainL domainL_num];
h = distanceJob(sourceData,targetData,domainData);
%%
test_rgbFile = FileList{2};
[~,nm]= fileparts(test_rgbFile);
test_maskFile = [oPath nm '.tif'];
M = imread(test_maskFile);
I = imread(test_rgbFile);
out = flattenMaskOverlay(I,M);

%%
sourceData.rgb_fileName;
sourceData.mask_fileName;
sourceData.pixelList;


targetData.rgb_fileName;
targetData.mask_fileName;
targetData.pixelList;


%%
samW = 50;
samN = samW*2 + 1;
% make sample square
[sam1,sam2] = ndgrid(linspace(-samW,samW,samN),linspace(-samW,samW,samN));
samSZ = size(sam1);
SAM = [sam1(:),sam2(:),ones(size(sam1(:)))];

% for each file
for e = 1:numel(FileList)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the file name
    [~,nm] = fileparts(FileList{e});
    % make the mask name
    maskName = [oPath nm '.tif'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read the image
    I = imread(FileList{e});
    I = double(I)/255;
    % read the mask
    M = imread(maskName);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % process the mask - get crop box
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make logical and get crop box
    M = logical(M);
    M = bwlarge(M);
    R = regionprops(M,'BoundingBox');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % crop the image and the mask
    I = imcrop(I,R(1).BoundingBox);
    M = imcrop(M,R(1).BoundingBox);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % make overlay
    out = flattenMaskOverlay(I,M);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make skeleton, endpoints, and branch points list
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make skeleton
    skeleton = bwmorph(M,'skeleton',inf);
    % find skeleton points
    skel = [];
    [skel(:,2),skel(:,1)] = find(skeleton);
    
    % find end points
    ep = [];
    endpoints = bwmorph(skeleton,'endpoints');
    [ep(:,2),ep(:,1)] = find(endpoints);
    [epi] = findIndexInSpace(skel,ep);
    
    % find branch points
    bp = [];
    branchpoints = bwmorph(skeleton,'branchpoints');
    [bp(:,2),bp(:,1)] = find(branchpoints);
    [bpi] = findIndexInSpace(skel,bp);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make distance transform
    DT = double(bwdist(~M));
    DT = imfilter(DT,fspecial('gaussian',[31 31],7),'replicate');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    [gDT1,gDT2] = gradient(DT);
    affine = zeros(1,3);
    affine(3) = 1;
    
    
    %M = cat(3,double(M),double(M),double(M));
    
    M = double(M);
    tmpI = zeros([samSZ size(I,3) size(skel,1)]);
    tmpM = zeros([samSZ size(M,3) size(skel,1)]);
   
    for p = 1:size(skel,1)
        n1 = gDT1(skel(p,2),skel(p,1));
        n2 = gDT2(skel(p,2),skel(p,1));
        t1 = -n2;
        t2 = n1;
        
        NOR = [n1 n2];
        NOR = NOR / norm(NOR);
        TAN = [t1 t2];
        TAN = TAN / norm(TAN);
        
        DIS = [skel(p,1) skel(p,2)];
        T = [[TAN',NOR',DIS'];affine];
        tSAM = mtimesx(T,SAM,'T')';
        
       
        tmp = ba_interp2(I,tSAM(:,1),tSAM(:,2));
        tmp = reshape(tmp,[samSZ size(I,3)]); 
        tmpI(:,:,:,p) = tmp;
        
        
       
        tmp = ba_interp2(M,tSAM(:,1),tSAM(:,2));
        tmp = reshape(tmp,[samSZ size(M,3)]);
        tmpM(:,:,:,p) = tmp;
        
        p
        size(skel,1)
    end
    
    
    %compute features
    NPP = 100;
    pi = zeros(size(tmpM,4),NPP);
    pSpace = linspace(0,norm(2*[samW samW]),NPP);
    parfor p = 1:size(tmpM,4)
        m = tmpM(:,:,:,p);
        i = tmpI(:,:,:,p);
        mn = m(:) / sum(m(:));
        mx = [];
        [mx(:,2),mx(:,1)] = find(m);
        Ex(p,:) = mtimesx(mn,'T',SAM);
        cp = Ex(p,:);
        cp = cp(1:2);
        m = bsxfun(@minus,mx,cp);
        dm = sum(m.*m,2).^.5;
        pi(p,:) = ksdensity(dm,pSpace);
        p
    end
    
    
    pi = bsxfun(@times,pi,sum(pi,2).^-1);
    en = -sum(log(pi+eps).*pi,2);
    [iS,iC,iU,iE,iL,iERR,iLAM] = PCA_FIT_FULL(pi,2);
    w = sweepPCA(iC,iE,iU,3*iLAM(1).^.5,1,5);
    [N] = hist3(iC,[1000 1000]);
    N = imfilter(N,fspecial('gaussian',[31 31],11),'replicate');
    plot(dX,iC(:,1),'.');
    
    
    
    dX = sum(Ex(:,1:2).*Ex(:,1:2),2).^.5;
    [sdX,sidx] = sort(dX);
    
    ng = 3;
    %tmpX = [dX iC(:,1)];
    tmpX = [dX en];
    tmpX = zscore(tmpX);
    kidx = kmeans(tmpX,ng);
    DC = [];
    CL = {'r.','b.','g.','k.','m.','c.','y.'};
    figure;
    for k = 1:ng
        fidx = find(kidx==k);
        subX = tmpX(fidx,:);
        plot(subX(:,1),subX(:,2),CL{k});
        hold on
        tmpU = mean(subX,1);
        
        MIN = min(subX,[],1);
        MAX = max(subX,[],1);
        tmpU = .5*[MAX - MIN] + MIN;
        
        subX = bsxfun(@minus,subX,tmpU);
        subX = sum(subX.*subX,2).^.5;
        [~,midx] = min(subX);
        
        DC = [DC , tmpI(:,:,:,fidx(midx))];
        
    end
    
    figure;
    imshow(DC,[]);
    hold on
    pause(1);
    drawnow
    xp = 10:size(tmpI,2):numel(CL)*size(tmpI,2);
    yp = 50*ones(size(xp));
    for p = 1:numel(CL)
        plot(xp(p),yp(p),CL{p});
    end
    
    % show clustered skeleton
    figure;
    imshow(I,[]);
    hold on
    for k = 1:ng
        fidx = find(kidx==k);
        plot(skel(fidx,1),skel(fidx,2),CL{k})
    end
    
    
    for p = 1:100:numel(sidx)
        imshow(tmpI(:,:,:,sidx(p)),[]);
        drawnow
    end
    
    
    close all
    for p = 1:100:size(tmpV,4)
        imshow(tmpV(:,:,:,p),[]);
        drawnow
    end
    
    
    %{
    % make adj matrix
    d = pdist(skel);
    d = squareform(d);
    dm = d <= 2^.5;
    d = double(dm).*d;
    
    % make graph from d
    g = graph(d);
    
  
    disp = true;
    if disp
        imshow(out,[]);
        hold on;
        drawnow
    end
    traceD = [];
    for i = 1:numel(epi)
        for j = 1:numel(bpi)
            [traceP{i,j},traceD(i,j)] = g.shortestpath(epi(i),bpi(j));
        end
        [~,sidx] = min(traceD(i,:));
        if disp
            
            plot(ep(i,1),ep(i,2),'g.');
            plot(bp(sidx,1),bp(sidx,2),'r.');
            tidx = traceP{i,sidx};
            %plot(skel(:,1),skel(:,2),'b.')
            plot(skel(tidx,1),skel(tidx,2),'k');
            drawnow
        end
    end
    %}
    
    
    
    
    
    
    out = flattenMaskOverlay(I,M);
    imshow(out,[]);
    hold on
    plot(skel(:,1),skel(:,2),'b.');
    hold off
    drawnow
end