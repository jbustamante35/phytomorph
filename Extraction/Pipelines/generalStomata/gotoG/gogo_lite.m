%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gather all maizeData
dataPath = ['/iplant/home/leakey_cyverse/maizeData/stomataTopoData%'];
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
cnt = 1;
for e = 1:numel(r)
    [~,~,ext] = fileparts(r(e).DATA_NAME);
    if strcmp(ext,'.nms')
        maizeFileList{cnt} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
        cnt = cnt + 1;
    end
end
%% test on new image
disp = false;
close all
for i = 2
    [res] = gogoCount(nptSZ,maizeFileList{i},npt,G,tE,tU,trainedNet2,tree,GLM,STEP,beta,nb,pnet,FDA,fIDX,fIDX_step,disp,'');
    STORE_Maize(:,:,i) = res;
end
%%
%close all
vS = reshape(STORE_Maize,[nptSZ(1:2) size(res,2),size(STORE_Maize,3)]);
disp = 1;
THL = [];
CNT_60 = [];
for i = 1:size(vS,4)
  close all
        if disp
            I = imread(maizeFileList{i});
            for e = 1:4
                I(1:initSize,:) = [];
                I = imrotate(I,90);
            end
        end




        sig = min(vS(:,:,1,i),[],3);
        
        sig = imfilter(sig,fspecial('gaussian',[5]),'replicate');
        mm = sig > .21111;
        mm = bwareaopen(mm,40);
        mmOLD = bwareaopen(mm,70);
        
        
        
        
        R = regionprops(logical(mm),'Centroid','Area','Perimeter','PixelIdxList');


        mm = zeros(size(mm));
        mm2 = zeros(size(mm));


        fidx = count([R.Area])==1 & count([R.Perimeter]+1)==1;
        %fidx = count([R.Area])==1;
        fidx2 = count([R.Area])==2 & count([R.Perimeter]+1)==2;

        fidx = count([R.Area])==1;
        fidx2 = count([R.Area])==2;


        fidx = find(fidx == 1);
        for f = 1:numel(fidx)
            mm(R(fidx(f)).PixelIdxList) = 1;
        end



        fidx2 = find(fidx2 == 1);
        for f = 1:numel(fidx2)
            mm2(R(fidx2(f)).PixelIdxList) = 1;
        end

        mm = logical(mm);
        mm2 = logical(mm2);

        %CNT(i) = size(R,1);
        CNT_60(i) = numel(fidx) + 2*numel(fidx2);






        if disp
            %[pth,nm,ext] = fileparts(QClist{i});
            %imwrite(imresize(I,size(mm)),[oPath nm '.tif']);
            out = flattenMaskOverlay(bindVec(imresize(I,size(mm))/255),mm,.55,'r');
            %out = flattenMaskOverlay(out,mmOLD,.55,'b');
            out = flattenMaskOverlay(out,mm2,.55,'g');
            %figure;imshow(sig,[]);
            figure;imshow(out,[]);
            title([num2str(CNT_60(i)) '--' num2str(GT{i,2}) '---->' num2str(i)]);

            drawnow
            waitforbuttonpress



    
            close all
            
            
        end
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create disk grid
G = [];
g = [];
N = [100 60];
R = [10 40 60];
for r = 1:numel(R)
    [g(:,:,1),g(:,:,2)] = ndgrid(linspace(-pi,pi,N(1)),linspace(0,R(r),N(2)));
    G(:,:,:,r) = cat(3,g(:,:,2).*sin(g(:,:,1)),g(:,:,2).*cos(g(:,:,1)));
end
%% create bug disk
G = [];
g = [];

NL = [5 5];
RL = [40];
GL = [];
g = [];
for l = 1:numel(RL)
    [g(:,:,1),g(:,:,2)] = ndgrid(linspace(-pi,pi,NL(1)),linspace(0,RL(l),NL(2)));
    GL(:,:,:,l) = cat(3,g(:,:,2).*sin(g(:,:,1)),g(:,:,2).*cos(g(:,:,1)));
end
NU = [100 20];
RU = [20];
GU = [];
g = [];
for u = 1:numel(RU)
    [g(:,:,1),g(:,:,2)] = ndgrid(linspace(-pi,pi,NU(1)),linspace(0,RU(u),NU(2)));
    GU(:,:,:,u) = cat(3,g(:,:,2).*sin(g(:,:,1)),g(:,:,2).*cos(g(:,:,1)));
    %[g(:,:,1),g(:,:,2)] = ndgrid(linspace(-RU(u),RU(u),NU(1)),linspace(-RU(u),RU(u),NU(1)));
    %GU(:,:,:,u) = cat(3,g(:,:,1),g(:,:,2));
    
end
G = zeros(size(GU,1),size(GU,2),size(GL,1),size(GL,2),size(GL,4),size(GU,4),2);
zoomCNT = 1;
for i1 = 1:size(GL,1)
    for i2 = 1:size(GL,2)
        for zoomLower = 1:size(GL,4)
            for zoomUpper = 1:size(GU,4)
                G(:,:,i1,i2,zoomLower,zoomUpper,1) = GL(i1,i2,1,zoomLower)+GU(:,:,1,zoomUpper);
                G(:,:,i1,i2,zoomLower,zoomUpper,2) = GL(i1,i2,2,zoomLower)+GU(:,:,2,zoomUpper);
                zoomCNT = zoomCNT + 1;
            end
        end
    end
    i1
end
%% create square grid
G = [];
g = [];
N = [50 60];
R = [20 40 60];
for r = 1:numel(R)
    [g(:,:,1),g(:,:,2)] = ndgrid(linspace(-R(r),R(r),N(1)),linspace(-R(r),R(r),N(1)));
    G(:,:,:,r) = cat(3,g(:,:,1),g(:,:,2));
end
%% create sample grid
imageFraction = .35;
wholeSize = 512;
initSize = abs(min(G(:)));
pt = [];
[pt(:,:,1),pt(:,:,2)] = ndgrid(initSize:round(wholeSize*imageFraction),initSize:round(wholeSize*imageFraction));
ptSZ = size(pt);
pt = reshape(pt,[prod(ptSZ(1:2)) ptSZ(3)]);
%% gather all sorghumData
dataPath = ['/iplant/home/leakey_cyverse/sorghumData/stomataTopoData%'];
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
cnt = 1;
for e = 1:numel(r)
    [~,~,ext] = fileparts(r(e).DATA_NAME);
    if strcmp(ext,'.nms')
        sorFileList{cnt} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
        cnt = cnt + 1;
    end
end
%% exclude from a folder
FilePath = '/mnt/tetra/nate/STOMATA/counts_March_29_2019/';
exFileList = {};
FileExt = {'jpg'};
exFileList = fdig(FilePath,exFileList,FileExt,1);
% remove contains data
for e = 1:numel(exFileList)
    toKP(e) = ~contains(exFileList{e},'raw');
end
exFileList = exFileList(toKP);
for e = 1:numel(exFileList)
    [~,exNM{e}] = fileparts(exFileList{e});
end
for e = 1:numel(sorFileList)
    [~,soNM{e}] = fileparts(sorFileList{e});
end
[C,iC] = setdiff(soNM,exNM);
sorFileList = sorFileList(iC);
%% generate tickets
TsorFileList = issueBulkTicket(sorFileList);
%% generate list for training
NI = 50;
trainFileList = {};
for e = 1:NI
    trainFileList{e} = sorFileList{e};
end
for e = 1:10
    trainFileList{end+1} = QClist{takeIDX(e)};
end
%% gather data for AI

%sData = zeros(size(pt,1),NI,30,20,25);
sData = zeros(size(pt,1),NI,size(G,1),size(G,2),size(G,4));
disp = false;
close all
for i = 1:NI
    I = imread(sorFileList{i});
    tmpSUB = zeros(size(pt,1),size(G,1),size(G,2),size(G,4));
    %tmpSUB = zeros(size(pt,1),30,20,25);
    parfor p = 1:size(pt,1)
        if disp
            imshow(I,[]);
            hold on
            plot(pt(p,2),pt(p,1),'r*')
        end
        
        
        
        % normal work
        tmpD = cat(3,G(:,:,1,:)+pt(p,1),G(:,:,2,:)+pt(p,2));
        subI = ba_interp2(I,squeeze(tmpD(:,:,2,:)),squeeze(tmpD(:,:,1,:)));
        %subI = bsxfun(@minus,subI,mean(subI,1));
        %subI = abs(fft(subI,[],1));
        
        %{
        % complex work
        tmpD = cat(7,G(:,:,:,:,:,:,1)+pt(p,1),G(:,:,:,:,:,:,2)+pt(p,2));
        subI = ba_interp2(I,tmpD(:,:,:,:,:,:,2),tmpD(:,:,:,:,:,:,1));
        %for r = 1:size(subI,3)
        %    imshow(squeeze(subI(:,:,r,20)),[]);
        %    drawnow
        %end
        %subI = bsxfun(@minus,subI,mean(subI,1));
        subI = abs(fft(subI,[],1));
        subI = subI(1:30,:,:,:);
        tmpSZ = size(subI);
        subI = reshape(subI,[tmpSZ(1:2) prod(tmpSZ(3:4))]);
        %}
        
        
        
        
        tmpSUB(p,:,:,:) = subI;
        
        if disp
            for m = 1:size(subI,3)
                tmp = subI(:,:,m);
                tmp = bindVec(tmp);
                subI(:,:,m) = tmp;
            end

            imshow(subI(:,:,1),[]);
            hold off
            drawnow
        end
        
        
        fprintf([num2str(p) ':' num2str(size(pt,1)) ':' num2str(i) ':' num2str(NI) '\n']);
        
    end
    sData(:,i,:,:,:) = tmpSUB;
end
%% sample for clicks
for i = 1:numel(trainFileList)
    I = imread(trainFileList{i});
    sam = ba_interp2(I,pt(:,2),pt(:,1));
    sam = reshape(sam,ptSZ(1:2));
    imshow(sam,[]);
    title(num2str(i));
    drawnow
    samI(:,:,i) = sam;
end
%% collect clicks
for i = 51:numel(trainFileList)
    [col{i},row{i},~] = impixel(samI(:,:,i),[]);
end
%% make masks
MSK = [];
close all
for i = 1:NI
    msk = zeros(size(samI,1),size(samI,2));
    for p = 1:numel(col{i})
        msk(row{i}(p),col{i}(p)) = 1;
    end
    msk = imdilate(msk,strel('disk',5,0));
    out = flattenMaskOverlay(samI(:,:,i),logical(msk));
    imshow(out,[]);
    drawnow
    MSK = [MSK;msk(:)];
end
%% compress samples
X = sData;
X = X(:,:,1:40,:,:);
szX = size(X);
%%
szX = size(sData);
X = reshape(sData,[prod(szX(1:2)) szX(3:5)]);
X = permute(X,[2 3 4 1]);
%% from X other place
szXX = size(X);
tData = reshape(X,[prod(szXX(1:3)) prod(szXX(4))])';
[tU,tE,tL] = PCA_FIT_FULLws(tData(find(MSK==1),:),10);
%[tU,tE,tL] = PCA_FIT_FULLws(tData,20);
[tC,err] = PCA_REPROJ(tData,tE,tU);
tC = [tC err];
%[tC] = PCA_REPROJ(tData,tE,tU);
%% show compressed data
%tCTMP = reshape(tC,[szX(1:2) 20]);
iC = reshape(tC,[ptSZ(1:2) NI size(tC,2)]);
for v = 1:size(iC,3)
    tmp = squeeze(iC(:,:,v,:));
    vw = [1 2 10];
    for k = 1:3
        tmp(:,:,k) = bindVec(tmp(:,:,vw(k)));
    end
    imshow(tmp(:,:,1:3),[]);
    waitforbuttonpress
end
%% rank features
fIDX = [];
for r = 1:100
    IDX1 = find(MSK==1);
    IDX0 = find(MSK==0);
    IDX0 = IDX0(randperm(numel(IDX0)));
    IDX = [IDX0(1:3*numel(IDX1));IDX1];
    [fIDX(:,r), Z] = rankfeatures(tC(IDX,:)', MSK(IDX)');
end
fNC = 4;
fIDX = fIDX(1:fNC,:);
UQ = unique(fIDX);
S = [];
for u = 1:numel(UQ)
    S(u) = sum(fIDX(:) == UQ(u));
end
[~,sidx] = sort(S,'descend');
fIDX = UQ(sidx(1:fNC))
fIDX = 1:size(tC,2);
%% train patten net on compressed data
IDX1 = find(MSK==1);
IDX0 = find(MSK==0);
IDX0 = IDX0(randperm(numel(IDX0)));
IDX = [IDX0(1:4*numel(IDX1));IDX1];
tT = tC;
pnet = patternnet(15);
pnet = train(pnet,tT(IDX,fIDX)',[logical(MSK(IDX)');~logical(MSK(IDX)')]);
%% train pls on compressed data
trainX = tC(IDX,fIDX);
trainY = MSK(IDX);
fprintf(['Starting training of PLS.\n']);
for l = 1:10
    [~,~,~,~, beta,~,~,~,~] = plsregress(trainX,trainY,l);
    yPre = [ones(numel(IDX),1) trainX]*beta;
    CV(l) = corr(yPre,trainY);
    fprintf(['Done fitting model (' num2str(l) ').\n'])
end
fprintf(['End training of PLS.\n']);
[~,~,~,~, beta,~,~,~,~] = plsregress(trainX,trainY,6);
%% tree
trainX = tC;
trainY = MSK;
IDX1 = find(MSK==1);
IDX0 = find(MSK==0);
IDX0 = IDX0(randperm(numel(IDX0)));
IDX = [IDX0(1:(2*numel(IDX1)));IDX1];
fprintf(['Starting training of TREE.\n']);
tree = fitctree(trainX(IDX,fIDX),trainY(IDX),'ScoreTransform','logit','OptimizeHyperparameters','auto');
fprintf(['End training of TREE.\n']);
%% nb
trainX = tC;
trainY = MSK;
IDX1 = find(MSK==1);
IDX0 = find(MSK==0);
IDX0 = IDX0(randperm(numel(IDX0)));
IDX = [IDX0(1:2*numel(IDX1));IDX1];
nb = fitcnb(trainX(IDX,fIDX), trainY(IDX));
%% glm
trainX = tC;
trainY = MSK;
IDX1 = find(MSK==1);
IDX0 = find(MSK==0);
IDX0 = IDX0(randperm(numel(IDX0)));
IDX = [IDX0(1:2*numel(IDX1));IDX1];
fprintf(['Starting training of GLM.\n']);
GLM = fitglm(trainX(IDX,fIDX),trainY(IDX),'quadratic','Distribution','binomial','link','logit','BinomialSize',1);
fprintf(['End training of GLM.\n']);
%% step
fIDX_step = [2 3 11];
fIDX_step = [1 3 8];
fIDX_step = [1 2 3];
fIDX_step = [2 8 10];
fIDX_step = [2 5 17 19];
fIDX_step = [1 3 5 11];
fIDX_step = [1 3 6 11];
fIDX_step = [1 3 4 11];
fIDX_step = [4 5 7 9];
fIDX_step = [7 8 9 11];
trainX = tC;
trainY = MSK;
IDX1 = find(MSK==1);
IDX0 = find(MSK==0);
IDX0 = IDX0(randperm(numel(IDX0)));
IDX = [IDX0(1:2*numel(IDX1));IDX1];
fprintf(['Starting training of STEP.\n']);
STEP = stepwiseglm(trainX(IDX,fIDX_step),trainY(IDX),'quadratic','Distribution','binomial','link','logit');
fprintf(['End training of STEP.\n']);
%% FDA
trainX = tC;
trainY = MSK;
IDX1 = find(MSK==1);
IDX0 = find(MSK==0);
IDX0 = IDX0(randperm(numel(IDX0)));
IDX = [IDX0(1:4*numel(IDX1));IDX1];
FDA = fitcdiscr(trainX(IDX,fIDX),trainY(IDX),'DiscrimType','quadratic');
%% downsample to create CNN data
szT = size(sData);
X = reshape(sData,[prod(szT(1:2)) szT(3:5)]);
X = permute(X,[2 3 4 1]);
X = X(1:40,:,:,:);
%% new gogo
szX = size(sData);
X = reshape(sData,[prod(szX(1:2)) szX(3:5)]);
X = permute(X,[2 3 4 1]);
%% create layers for CNN
layers = [
    imageInputLayer([size(X,1),size(X,2),size(X,3)])

    convolution2dLayer(5,15,'Padding','same')
    batchNormalizationLayer
    reluLayer
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,5,'Padding',1)
    batchNormalizationLayer
    reluLayer
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,5,'Padding',1)
    batchNormalizationLayer
    reluLayer
    maxPooling2dLayer(2,'Stride',2)

    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];

%{
layers = [
    imageInputLayer([size(X,1),size(X,2),size(X,3)])

    convolution2dLayer(5,5,'Padding','same')
    batchNormalizationLayer
    reluLayer
    maxPooling2dLayer(6,'Stride',4)
    

    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];
    %}
options = trainingOptions('sgdm','ExecutionEnvironment','parallel','MaxEpochs',20,'InitialLearnRate',0.001,'Verbose',false,'Plots','training-progress');
%% train CNN
IDX1 = find(MSK==1);
IDX0 = find(MSK==0);
IDX0 = IDX0(randperm(numel(IDX0)));
IDX = [IDX0(1:2*numel(IDX1));IDX1];
trainedNet_60 = trainNetwork(X(:,:,:,IDX),categorical(MSK(IDX)),layers,options);
%%
options = trainingOptions('sgdm','ExecutionEnvironment','cpu','MaxEpochs',20,'InitialLearnRate',0.001,'Verbose',false,'Plots','training-progress');
trainedNet2 = trainNetwork(X(:,:,:,IDX),categorical(MSK(IDX)),trainedNet.Layers,options);
%%
trainedNet = trainNetwork(X,categorical(MSK),layers,options);
%% simulate pattern network
close all
for i = 1:NI
    Ypre = sim(pnet,tT(:,fIDX)');
    Ypre = Ypre(1,:);
    Ypre = reshape(Ypre,szT(1:2));
    Ypre = reshape(Ypre,[ptSZ(1:2) NI]);
    imshow(Ypre(:,:,i),[]);
    drawnow
    figure;
    imshow(samI(:,:,i),[]);
    waitforbuttonpress
    close all
end
%% sim tree
close all
for i = 1:NI
    [~,Ypre] = tree.predict(tC(:,fIDX));
    Ypre = reshape(Ypre(:,2),szT(1:2));
    Ypre = reshape(Ypre,[ptSZ(1:2) NI]);
    imshow(Ypre(:,:,i),[]);
    drawnow
    figure;
    imshow(samI(:,:,i),[]);
    waitforbuttonpress
    close all
end
%% NB 
close all
for i = 1:NI
    [~,Ypre] = nb.predict(tC(:,fIDX));
    Ypre = reshape(Ypre(:,2),szT(1:2));
    Ypre = reshape(Ypre,[ptSZ(1:2) NI]);
    imshow(Ypre(:,:,i),[]);
    drawnow
    figure;
    imshow(samI(:,:,i),[]);
    waitforbuttonpress
    close all
end

%% sim FDA
close all
for i = 1:NI
    Ypre = FDA.predict(tC(:,fIDX));
    Ypre = reshape(Ypre,szT(1:2));
    Ypre = reshape(Ypre,[ptSZ(1:2) NI]);
    imshow(Ypre(:,:,i),[]);
    drawnow
    figure;
    imshow(samI(:,:,i),[]);
    waitforbuttonpress
    close all
end
%% simulate CNN
Ypre = trainedNet_60.predict(X,'MiniBatch',4096);
Ypre = Ypre(:,2);
Ypre = reshape(Ypre,szT(1:2));
Ypre = reshape(Ypre,[ptSZ(1:2) NI]);
%%
close all
for i = 1:NI
imshow(Ypre(:,:,i),[]);
drawnow
figure;
imshow(samI(:,:,i),[]);


mm = Ypre(:,:,i) > .5;
mm = bwareaopen(mm,50);
out = flattenMaskOverlay(bindVec(imresize(samI(:,:,i),size(mm))/255),mm);
imshow(out,[]);
waitforbuttonpress
close all
end
%%
imageFraction = 1;
wholeSize = 512;
initSize = abs(min(G(:)));
npt = [];
[npt(:,:,1),npt(:,:,2)] = ndgrid(initSize:(wholeSize-initSize),initSize:(wholeSize-initSize));
nptSZ = size(npt);
npt = reshape(npt,[prod(nptSZ(1:2)) nptSZ(3)]);
%% test on new image
disp = false;
close all
for i = 156:numel(QClist)
    [res] = gogoCount(QClist{i},npt,G,tE,tU,trainedNet_60,tree,GLM,STEP,beta,nb,pnet,FDA,fIDX,fIDX_step,disp);
    STORE_60(:,:,i) = res;
end
%%
QClist = issueBulkTicket(QClist);
%% here runs condor on all ********

rPath = '/iplant/home/phytomorphuser/counts_March_18_2019/';
[rPath,iticket] = issueTicket(rPath(1:end-1),10*numel(TsorFileList),'write');


ggC = cFlow('gogoCount');
ggC.setMCRversion('v930');
ggC.setMemory('8000');
res = {};
for i = 1:10%numel(TsorFileList)
   res{i} = ggC(nptSZ,TsorFileList{i},npt,G,tE,tU,trainedNet,tree,GLM,STEP,beta,nb,pnet,FDA,fIDX,fIDX_step,disp,rPath);
    i
end
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
ggC.submitDag(auth,150,150);
%%
ff = @(X,Y)gogoCount(nptSZ,X,npt,G,tE,tU,trainedNet,tree,GLM,STEP,beta,nb,pnet,FDA,fIDX,fIDX_step,disp,Y);
pf = partialFunction(ff,'March_20_2019_sorghumStomata');
pf.publish();
%%
myPFwrapper('March_20_2019_sorghumStomata',TsorFileList{i},rPath);
%% try patrial function on CHTC
rPath = '/iplant/home/phytomorphuser/counts_March_29_2019/';
system(['imkdir ' rPath]);
[rPath,iticket] = issueTicket(rPath(1:end-1),10*numel(TsorFileList),'write');
ggC = cFlow('myPFwrapper');
ggC.setMCRversion('v930');
ggC.setMemory('8000');
res = {};
for i = 1:numel(TsorFileList)
   res{i} = ggC('March_20_2019_sorghumStomata',TsorFileList{i},rPath);
   fprintf(['Done building job:' num2str(i) '\n'])
end
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
ggC.submitDag(auth,150,150);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read new ground truths
GT = readtext('~/Ground truth 393x393 images .csv');
GT(1,:) = [];
CHECK = cell2mat(GT(:,2));
for e = 1:size(GT,1)
    [pth,nm,ext] = fileparts(QClist{e});
    toTest = GT{e,1}(3:end);
    if ~strcmp(toTest,[nm])
        nm
        toTest
        e
    end
end
%%
%close all
vS = reshape(STORE_60,[nptSZ(1:2) size(res,2),size(STORE_60,3)]);
vS = reshape(STOREs,[nptSZ(1:2) size(res,2),size(STOREs,3)]);
%vS = reshape(min(cat(4,STORE_60,STOREs),[],4),[nptSZ(1:2) size(res,2),size(STORE_60,3)]);
disp = 0;
MM = 3:21;
%MM = linspace(.1,.5,10);
TH = linspace(.1,.6,10);
THL = [];
CNT_60 = [];
%for m = 1:numel(MM)
%  for l = 1:numel(TH)
%oPath = ['/mnt/tetra/nate/forJohn/'];
CNT_60 = [];
for i = 1:size(vS,4)
  close all
        if disp
            I = imread(QClist{i});
            for e = 1:4
                I(1:initSize,:) = [];
                I = imrotate(I,90);
            end
        end


        %sig = mean(vS(:,:,:,i),3);
        %sig = .5*(max(vS(:,:,:,i),[],3) + min(vS(:,:,:,i),[],3));
        %sig = mean(vS(:,:,:,i),3);
        %{
        G = sort(vS(:,:,:,i),3);
        sig = prod(vS(:,:,:,i),3);
        G = reshape(G,[size(G,1)*size(G,2) size(G,3)]);
        [gS,gC,gU,gE,gL,gERR,gLAM] = PCA_FIT_FULL(G,3);
        gC = PCA_REPROJ(G,gE,gU);
        gC = reshape(gC,[size(mm) size(gC,2)]);
        sig = gC;
        %}



        sig = min(vS(:,:,1,i),[],3);
        
        %sig = imfilter(sig,fspecial('gaussian',[71 71]),'replicate');
        %sig = imfilter(sig,fspecial('gaussian',[MM(m)]),'replicate');
        
        %sig = imfilter(sig,fspecial('gaussian',[3]),'replicate');
        sig = imfilter(sig,fspecial('gaussian',[5]),'replicate');
        %mm = sig > TH(l);
        %mm = sig > 0.53;
        mm = sig > .21111;
        %mm = sig > 0.1556;
        mm = bwareaopen(mm,40);
        mmOLD = bwareaopen(mm,70);
        
        
        
        
        %{
        sig = min(vS(:,:,1,i),[],3);
        sig = imfilter(sig,fspecial('gaussian',[21]),'replicate');
        LEVEL = MM(m);
        mm = sig - LEVEL;
        RECON = imreconstruct(mm,sig);
        BO = (mm - RECON) > - LEVEL;
        BO = BO & mean(vS(:,:,1,i),3) > TH(l);
        %BO = imfill(BO,'holes');
        %BO = bwareaopen(BO,80);
        mm = BO;
         mm = bwareaopen(mm,40);
        %}
        
        
        
        
        
        R = regionprops(logical(mm),'Centroid','Area','Perimeter','PixelIdxList');


        mm = zeros(size(mm));
        mm2 = zeros(size(mm));


        fidx = count([R.Area])==1 & count([R.Perimeter]+1)==1;
        %fidx = count([R.Area])==1;
        fidx2 = count([R.Area])==2 & count([R.Perimeter]+1)==2;

        fidx = count([R.Area])==1;
        fidx2 = count([R.Area])==2;


        fidx = find(fidx == 1);
        for f = 1:numel(fidx)
            mm(R(fidx(f)).PixelIdxList) = 1;
        end



        fidx2 = find(fidx2 == 1);
        for f = 1:numel(fidx2)
            mm2(R(fidx2(f)).PixelIdxList) = 1;
        end

        mm = logical(mm);
        mm2 = logical(mm2);
        %mm = bwareaopen(mm,70);
        %{
        mm = imdilate(sig,strel('disk',11,0)) == sig & sig > .50;
        mm = imdilate(mm,strel('disk',5,0));
        %}

        %{
        sig = mean(vS(:,:,7,i),3);
        mm = bwareaopen(sig > .5,80);
        mmBIG = logical(bwareaopen(mm,500));
        mm = logical(mm - mmBIG);
        %}
        %mm = imdilate(sig,strel('disk',21,0)) == sig & sig > .50;
        %mm = imdilate(mm,strel('disk',5,0));


        %{
        LEVEL = 1;
        mm = sig - LEVEL;
        RECON = imreconstruct(mm,sig);
        BO = (mm - RECON) > - LEVEL;
        BO = BO & mean(vS(:,:,:,i),3) > .25;
        BO = imfill(BO,'holes');
        %BO = bwareaopen(BO,80);
        mm = BO;
        %}
        %mm = mean(sig,3) > .5;
        %mm = bwareaopen(mm,40);


        %R = regionprops(logical(mm),'Centroid');


        %{

        tmpM = zeros(size(mm));
        for e = 1:numel(R)
            LOC(e,:) = R(e).Centroid;
            tmpM(round(LOC(e,2)),round(LOC(e,1))) = 1;
        end
        DIST = bwdist(tmpM);

        DIS = squareform(pdist(LOC));
        MM = [];
        for e = 1:size(LOC,1)
            tmp = DIS(e,:);
            tmp(e) = inf;
            [MM(e)] = min(tmp);
        end
        lambda = 1/mean(MM);
        PROB = 1- exp(-lambda*DIST);
        imshow(mm,[]);

        %}


        %CNT(i) = size(R,1);
        CNT_60(i) = numel(fidx) + 2*numel(fidx2);






        if disp
            %[pth,nm,ext] = fileparts(QClist{i});
            %imwrite(imresize(I,size(mm)),[oPath nm '.tif']);
            out = flattenMaskOverlay(bindVec(imresize(I,size(mm))/255),mm,.55,'r');
            %out = flattenMaskOverlay(out,mmOLD,.55,'b');
            out = flattenMaskOverlay(out,mm2,.55,'g');
            %figure;imshow(sig,[]);
            figure;imshow(out,[]);
            title([num2str(CNT_60(i)) '--' num2str(GT{i,2}) '---->' num2str(i)]);

            drawnow
            waitforbuttonpress



    
            prompt = {'Enter deltaCount:'};
            dlg_title = 'Input';
            num_lines = 1;
            defaultans = {'0'};
            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
            deltaC(i) = str2num(answer{1});

           % corr(CNT_60(1:i)',deltaC(1:i)'+CHECK(1:i))

            close all
            
            
        end
    end
    
    
 %{
    %THL(l,m) = corr(CNT_60',CHECK(1:numel(deltaC))+deltaC');
    THL(l,m) = corr(CNT_60',CHECK2);
    max(THL)
    if l > 2 & m > 2
        mesh(THL);
        drawnow
        hold on
    end
    
   
    end
end
    %}
%%
csvwrite('/mnt/tetra/nate/forJohn_march7.csv',[CHECK CHECK2 CNT_60'])
%% look for the worst example(s)
[JJ,takeIDX] = sort(abs(CNT_60'-CHECK),'descend');
%% clip the image from test
%res = reshape(res(:,2),nptSZ(1:2));
res = reshape(res,[nptSZ(1:2) size(res,2)]);
for e = 1:4
    I(1:initSize,:) = [];
    I = imrotate(I,90);
end
%% smooth
for e = 1:size(res,3)
    sres(:,:,e) = imfilter(res(:,:,e),fspecial('gaussian',[31 31]),'replicate');
end
%%
I = imread(QClist{i});
    for e = 1:4
        I(1:initSize,:) = [];
        I = imrotate(I,90);
    end
%% show the test image and results
close all

sig = min(res(:,:,1),[],3);
%sig = imfilter(sig,fspecial('gaussian',[71 71]),'replicate');
sig = imfilter(sig,fspecial('gaussian',[21 21]),'replicate');
mm = imdilate(sig,strel('disk',11,0)) == sig & sig > .50;
mm = imdilate(mm,strel('disk',5,0));

R = regionprops(logical(mm),'Centroid');
    
figure;imshow(mean(res(:,:,1),3),[]);
figure;imshow(min(res(:,:,1),[],3),[]);
figure;imshow(I,[]);
mm = mean(res(:,:,1),3) > .5;
mm = bwareaopen(mm,80);
out = flattenMaskOverlay(bindVec(imresize(I,size(mm))/255),mm);
figure;imshow(out,[]);
